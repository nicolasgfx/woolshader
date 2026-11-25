# SPDX-License-Identifier: GPL-3.0-or-later
"""Convert a photon_path_tracer package (scene.json, PLY meshes, .pptstrands fibres) into a .blend.

    blender -b --factory-startup --python blender/ppt_package_to_blend.py -- <package> <out.blend>

The fibres become native Curves with the package's points and radii, and a part's fibre shape (the wool
fibres' thin ends) goes into the radii as the renderer applies it at compile. Cycles draws them as round
linear segments ('THICK_LINEAR'), which is the renderer's OPTIX_PRIMITIVE_TYPE_ROUND_LINEAR.
Everything is baked from metres Y-up to metres Z-up: (x, y, z) -> (x, -z, y), a rotation.

The viewer's Export step runs this script for its Blender format (src/export/format_blend.cpp). It
handles what a viewer export writes: hair, fabric (with its translucency and its highlight), plastic, coated_diffuse
and diffuse materials, identity transforms, an analytical environment. Anything else stops the conversion with a message, so nothing
is dropped silently.

Fabric has no exact Cycles equivalent. The renderer's BRDF is
    f = rho/pi (1 - c E(mo)) (1 - c E(mi)) + c C D V        (Charlie sheen, fabric_bsdf.h)
and a node cannot see mi. Under a constant dome only the directional albedo counts, so the group
below reproduces that: the diffuse term replaces (1 - c E(mi)) by its cosine-weighted mean K, and the
sheen term is Cycles' Microfiber lobe times g(mo), which fits it to the Charlie albedo.

A fabric may carry the plastic's highlight under its sheen (the renderer's R1, photon_path_tracer dd8870fda), layered
from outside in: sheen, highlight, diffuse body. The group draws the highlight as a Principled BSDF without a base
colour behind the sheen factor, crossed in and out, and gives the body and its translucency what a mirror of the
highlight's F0 does not reflect: the renderer's body sees the highlight's true directional albedo, so this is an
approximation, close for the yarn's F0 of about 0.03.
"""
import hashlib
import json
import math
import os
import sys

import bpy
import numpy as np
from mathutils import Matrix, Vector

# Fibres deep in the mass see less of the dome, mostly near their normal, where 1 - E is larger: the
# analytic K is multiplied by this, fitted in situ against a render of heart_knit by the renderer.
DIFFUSE_IN_SITU_GAIN = 1.0521
# Microfiber at the renderer's roughness + 0.2 has the shape closest to Charlie (furnace test at 0.5,
# rms 0.011).
MICROFIBER_ROUGHNESS_OFFSET = 0.2
# g(mu) = Charlie albedo / Microfiber albedo from the same furnace test, times a gain of 1.055 from
# the same in-situ fit (clipped to 1 at mu = 0). Measured at sheen roughness 0.5, the yarn's default;
# other roughnesses keep this table, so they are approximate.
SHEEN_FIT = [(0.0, 1.0), (0.1, 0.950), (0.2, 0.899), (0.35, 0.886), (0.5, 0.908),
             (0.65, 0.932), (0.8, 0.932), (0.95, 0.880), (1.0, 0.844)]

# The renderer's default display curve (Narkowicz ACES fit, display_kernels.cu), then the sRGB OETF.
TONE_A, TONE_B, TONE_C, TONE_D, TONE_E = 2.51, 0.03, 2.43, 0.59, 0.14

STRAND_MIN_RADIUS = 1e-6  # metres: the renderer's kStrandMinRadius, document/asset_data.h

DIELECTRIC_F0 = 0.04  # F0 at ior 1.5: the renderer's and Blender's specular 0.5 (F0 = 2 x 0.04 x specular)
# The highlight's Principled BSDF sockets, fed per material (connect_highlight) into the fabric group.
HIGHLIGHT_SOCKETS = ("Roughness", "Specular IOR Level", "Anisotropic", "Anisotropic Rotation", "Tangent")
# The specular block's keys the conversion models (the plastic's and the fabric's, document/material.h).
HIGHLIGHT_KEYS = {"specular", "roughness", "roughness_y", "aniso_rotation", "roughness_tex", "reflectance_tex",
                  "color", "specular_tint", "coverage"}


def to_z_up(points):
    """(n, 3) Y-up array -> Z-up array."""
    return np.stack([points[:, 0], -points[:, 2], points[:, 1]], axis=1)


def read_ply(path):
    """Binary PLY as written by export_scene.cpp. Returns (vertex array, triangles, property names)."""
    with open(path, "rb") as fh:
        data = fh.read()
    end = data.find(b"end_header\n")
    if end < 0:
        sys.exit(f"{path}: no PLY header")
    lines = data[:end].decode("ascii").split("\n")
    if "format binary_little_endian 1.0" not in lines:
        sys.exit(f"{path}: not little-endian binary")
    vertex_count = next(int(l.split()[2]) for l in lines if l.startswith("element vertex"))
    face_count = next(int(l.split()[2]) for l in lines if l.startswith("element face"))
    names = [l.split()[2] for l in lines if l.startswith("property float")]
    start = end + len(b"end_header\n")
    vertices = np.frombuffer(data, "<f4", vertex_count * len(names), start).reshape(-1, len(names))
    faces = np.frombuffer(data, np.dtype([("n", "u1"), ("v", "<u4", 3)]), face_count,
                          start + vertices.nbytes)
    if faces.size and (faces["n"] != 3).any():
        sys.exit(f"{path}: faces other than triangles")
    if start + vertices.nbytes + faces.nbytes != len(data):
        sys.exit(f"{path}: unexpected file size")
    return vertices, faces["v"], names


def read_strands(path):
    """.pptstrands v1. Returns (positions (V, 3), radii (V), first (N + 1)), Y-up metres."""
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:8] != b"PPTSTRND":
        sys.exit(f"{path}: not a strands file")
    version, basis = np.frombuffer(data, "<u4", 2, 8)
    vertex_count, strand_count, _ = (int(x) for x in np.frombuffer(data, "<u8", 3, 16))
    if version != 1 or basis != 0:
        sys.exit(f"{path}: version {version} basis {basis}, expected 1 and 0 (round linear)")
    if len(data) != 40 + 16 * vertex_count + 4 * (strand_count + 1):
        sys.exit(f"{path}: unexpected file size")
    positions = np.frombuffer(data, "<f4", 3 * vertex_count, 40).reshape(-1, 3)
    radii = np.frombuffer(data, "<f4", vertex_count, 40 + 12 * vertex_count)
    first = np.frombuffer(data, "<u4", strand_count + 1, 40 + 16 * vertex_count)
    return positions, radii, first


def build_mesh(name, path, material):
    vertices, triangles, names = read_ply(path)
    mesh = bpy.data.meshes.new(name)
    mesh.vertices.add(len(vertices))
    mesh.loops.add(3 * len(triangles))
    mesh.polygons.add(len(triangles))
    mesh.vertices.foreach_set("co", to_z_up(vertices[:, :3]).ravel())
    mesh.polygons.foreach_set("loop_start", np.arange(len(triangles), dtype=np.int32) * 3)
    mesh.loops.foreach_set("vertex_index", triangles.ravel().astype(np.int32))
    mesh.update(calc_edges=True)
    mesh.shade_smooth()
    # The stored normals are the yarn's, not what the triangles would give: keep them.
    mesh.normals_split_custom_set_from_vertices(to_z_up(vertices[:, 3:6]))
    if "u" in names:
        uv = vertices[:, [names.index("u"), names.index("v")]]
        layer = mesh.uv_layers.new(name="UVMap")
        layer.uv.foreach_set("vector", uv[triangles.ravel()].ravel())
    mesh.materials.append(material)
    return bpy.data.objects.new(name, mesh)


def shaped_radii(name, positions, radii, first, shape):
    """The radii the renderer draws (document/fibre_shape.cpp): r diameter_scale (1 + t (tip_scale - 1)),
    t the normalised arc length along the strand."""
    if shape["thickness"] != "imported" or shape["diameter_variation"] != 0.0:
        sys.exit(f"{name}: part_fibre_shapes other than imported radii without variation are not handled")
    steps = np.linalg.norm(np.diff(positions.astype(np.float64), axis=0), axis=1)
    steps[first[1:-1] - 1] = 0.0  # from one strand's last point to the next strand's first
    along = np.concatenate(([0.0], np.cumsum(steps)))
    strand = np.repeat(np.arange(len(first) - 1), np.diff(first))
    start = along[first[:-1]][strand]
    length = along[first[1:] - 1][strand] - start
    t = np.divide(along - start, length, out=np.zeros_like(along), where=length > 0)
    shaped = radii * shape["diameter_scale"] * (1.0 + t * (shape["tip_scale"] - 1.0))
    return np.maximum(shaped, STRAND_MIN_RADIUS).astype(np.float32)


def build_curves(name, path, material, shape):
    positions, radii, first = read_strands(path)
    if shape is not None:
        radii = shaped_radii(name, positions, radii, first, shape)
    curves = bpy.data.hair_curves.new(name)
    curves.add_curves(np.diff(first).tolist())
    curves.set_types(type='POLY')
    curves.position_data.foreach_set("vector", to_z_up(positions).ravel())
    radius = curves.attributes.new("radius", 'FLOAT', 'POINT')
    radius.data.foreach_set("value", radii)
    curves.materials.append(material)
    return bpy.data.objects.new(name, curves)


def clamped_output(nodes, links, value):
    node = nodes.new("ShaderNodeClamp")
    links.new(value, node.inputs["Value"])
    return node.outputs[0]


def sheen_energy_coefficients(roughness):
    """A, B, P of E(mu) = clamp(A (1 - sqrt(mu))^P + B), the fit of the renderer's lobes.h."""
    a = max(roughness, 0.2)
    return (((-1.378665 * a + 3.385171) * a - 3.029331) * a + 1.534094,
            ((0.103124 * a - 0.339355) * a + 0.456232) * a - 0.088798,
            ((-1.267677 * a + 2.937879) * a - 2.261111) * a + 1.526667)


def sheen_diffuse_scale(roughness, coverage):
    """K = 2 * integral of (1 - c E(mu)) mu dmu over [0, 1] by the midpoint rule, with the in-situ gain."""
    a, b, p = sheen_energy_coefficients(roughness)
    mu = (np.arange(4096) + 0.5) / 4096
    energy = np.clip(a * (1.0 - np.sqrt(mu)) ** p + b, 0.0, 1.0)
    return float(2.0 * np.mean((1.0 - coverage * energy) * mu)) * DIFFUSE_IN_SITU_GAIN


def build_fabric_group(roughness, coverage, transmit_share, highlight):
    """The fabric shader of one sheen roughness, coverage and transmit share, with or without the highlight.
    Inputs: Albedo, Sheen Color, Normal; with the highlight also the Principled BSDF's sockets HIGHLIGHT_SOCKETS
    (fed by connect_highlight).

    The renderer's translucency (TranslucencyParams, document/material.h): W of the body's albedo A
    leaves through the far side as Lambert transmission and A (1 - W) is reflected, one budget; the
    sheen stays outside, crossed where the light enters and where it leaves. Here the body's diffuse
    term keeps (1 - W) and a Translucent BSDF carries W, both behind the same sheen factor. On a closed
    yarn tube the light crosses the near wall inward and the far wall outward, as the renderer's bundle
    proxy. W = 0 builds the group as before, node for node.

    The highlight (when there is one): the renderer's fabric reflects with the dielectric highlight what reaches it
    through the sheen, (1 - c E(mu)) on the way in and again on the way out, and passes the rest on to the body.
    Here a Principled BSDF with a black base colour carries the highlight times (1 - c E(mu))^2, and the body's
    diffuse and translucent terms keep (1 - F0), F0 = 2 F0(1.5) x Specular IOR Level as Blender and the renderer
    both define it. Without a highlight the group is built as before, node for node."""
    name = f"PPT Fabric r{roughness:g} c{coverage:g}" + (f" t{transmit_share:g}" if transmit_share > 0 else "") + \
        (" h" if highlight else "")
    tree = bpy.data.node_groups.new(name, "ShaderNodeTree")
    tree.interface.new_socket("Albedo", in_out='INPUT', socket_type='NodeSocketColor')
    tree.interface.new_socket("Sheen Color", in_out='INPUT', socket_type='NodeSocketColor')
    tree.interface.new_socket("Normal", in_out='INPUT', socket_type='NodeSocketVector')
    if highlight:
        for socket in HIGHLIGHT_SOCKETS:
            tree.interface.new_socket(socket, in_out='INPUT',
                                      socket_type='NodeSocketVector' if socket == "Tangent" else 'NodeSocketFloat')
    tree.interface.new_socket("Shader", in_out='OUTPUT', socket_type='NodeSocketShader')
    nodes, links = tree.nodes, tree.links
    group_in = nodes.new("NodeGroupInput")
    group_out = nodes.new("NodeGroupOutput")

    def math_node(operation, a, b=None):
        node = nodes.new("ShaderNodeMath")
        node.operation = operation
        for socket, value in zip(node.inputs, (a, b)):
            if isinstance(value, bpy.types.NodeSocket):
                links.new(value, socket)
            elif value is not None:
                socket.default_value = value
        return node.outputs[0]

    def scale_color(color, factor):
        node = nodes.new("ShaderNodeVectorMath")
        node.operation = 'SCALE'
        links.new(color, node.inputs[0])
        links.new(factor, node.inputs["Scale"])
        return node.outputs[0]

    geometry = nodes.new("ShaderNodeNewGeometry")
    dot = nodes.new("ShaderNodeVectorMath")
    dot.operation = 'DOT_PRODUCT'
    links.new(group_in.outputs["Normal"], dot.inputs[0])
    links.new(geometry.outputs["Incoming"], dot.inputs[1])
    mu = clamped_output(nodes, links, dot.outputs["Value"])

    # (1 - c E(mu)) K
    energy_a, energy_b, energy_p = sheen_energy_coefficients(roughness)
    root = math_node('SUBTRACT', 1.0, math_node('SQRT', mu))
    energy = math_node('ADD', math_node('MULTIPLY', math_node('POWER', root, energy_p), energy_a), energy_b)
    energy = math_node('MULTIPLY', clamped_output(nodes, links, energy), coverage)
    through_sheen = math_node('SUBTRACT', 1.0, energy)
    keep = math_node('MULTIPLY', through_sheen, sheen_diffuse_scale(roughness, coverage))
    if highlight:  # the body keeps what the highlight's mirror does not reflect: 1 - F0
        f0 = math_node('MULTIPLY', group_in.outputs["Specular IOR Level"], 2.0 * DIELECTRIC_F0)
        keep = math_node('MULTIPLY', keep, math_node('SUBTRACT', 1.0, clamped_output(nodes, links, f0)))

    fit = nodes.new("ShaderNodeFloatCurve")
    curve = fit.mapping.curves[0]
    while len(curve.points) < len(SHEEN_FIT):
        curve.points.new(0.5, 0.5)
    for point, (x, y) in zip(curve.points, SHEEN_FIT):
        point.location = (x, y)
    fit.mapping.update()
    links.new(mu, fit.inputs["Value"])

    diffuse = nodes.new("ShaderNodeBsdfDiffuse")
    reflected = keep if transmit_share <= 0 else math_node('MULTIPLY', keep, 1.0 - transmit_share)
    links.new(scale_color(group_in.outputs["Albedo"], reflected), diffuse.inputs["Color"])
    links.new(group_in.outputs["Normal"], diffuse.inputs["Normal"])
    sheen = nodes.new("ShaderNodeBsdfSheen")
    sheen.distribution = 'MICROFIBER'
    sheen.inputs["Roughness"].default_value = min(roughness + MICROFIBER_ROUGHNESS_OFFSET, 1.0)
    fitted = math_node('MULTIPLY', fit.outputs[0], coverage)
    links.new(scale_color(group_in.outputs["Sheen Color"], fitted), sheen.inputs["Color"])
    links.new(group_in.outputs["Normal"], sheen.inputs["Normal"])
    add = nodes.new("ShaderNodeAddShader")
    links.new(diffuse.outputs[0], add.inputs[0])
    links.new(sheen.outputs[0], add.inputs[1])
    if transmit_share > 0:
        translucent = nodes.new("ShaderNodeBsdfTranslucent")
        links.new(scale_color(group_in.outputs["Albedo"], math_node('MULTIPLY', keep, transmit_share)), translucent.inputs["Color"])
        links.new(group_in.outputs["Normal"], translucent.inputs["Normal"])
        with_transmission = nodes.new("ShaderNodeAddShader")
        links.new(add.outputs[0], with_transmission.inputs[0])
        links.new(translucent.outputs[0], with_transmission.inputs[1])
        add = with_transmission
    if highlight:
        principled = nodes.new("ShaderNodeBsdfPrincipled")
        principled.inputs["Base Color"].default_value = (0.0, 0.0, 0.0, 1.0)  # the highlight alone: the body is above
        principled.inputs["IOR"].default_value = 1.5
        for socket in HIGHLIGHT_SOCKETS:
            links.new(group_in.outputs[socket], principled.inputs[socket])
        links.new(group_in.outputs["Normal"], principled.inputs["Normal"])
        # Through the sheen in and out; a Mix Shader with an empty first input scales the closure by its factor.
        behind_sheen = nodes.new("ShaderNodeMixShader")
        links.new(math_node('MULTIPLY', through_sheen, through_sheen), behind_sheen.inputs["Fac"])
        links.new(principled.outputs[0], behind_sheen.inputs[2])
        with_highlight = nodes.new("ShaderNodeAddShader")
        links.new(add.outputs[0], with_highlight.inputs[0])
        links.new(behind_sheen.outputs[0], with_highlight.inputs[1])
        add = with_highlight
    links.new(add.outputs[0], group_out.inputs["Shader"])
    return tree


def load_packed_image(package, relative_path, colorspace):
    """The package's image, packed: the package is temporary, the .blend must hold the pixels.
    Materials share one image: colourwork has a yarn material per colour, all with the same ply map."""
    stored_path = "//" + relative_path
    for image in bpy.data.images:
        if image.filepath_raw == stored_path:
            return image
    image = bpy.data.images.load(os.path.join(package, relative_path))
    image.colorspace_settings.name = colorspace
    image.pack()
    image.filepath_raw = stored_path
    return image


def albedo_output(spec, package, nodes, links):
    """The diffuse colour: the albedo, times its texture if there is one (sRGB, as the renderer reads it)."""
    albedo = spec["diffuse"]
    constant = nodes.new("ShaderNodeRGB")
    constant.outputs[0].default_value = (*albedo["albedo"], 1.0)
    if "albedo_tex" not in albedo:
        return constant.outputs[0]
    texture = nodes.new("ShaderNodeTexImage")
    texture.image = load_packed_image(package, albedo["albedo_tex"], 'sRGB')
    product = nodes.new("ShaderNodeVectorMath")
    product.operation = 'MULTIPLY'
    links.new(constant.outputs[0], product.inputs[0])
    links.new(texture.outputs["Color"], product.inputs[1])
    return product.outputs[0]


def normal_output(spec, package, nodes, links, on_strands):
    """The shading normal. On strands the renderer's fibre frame: the view direction made perpendicular
    to the axis, so a fibre has no terminator; Hair Info computes the same vector. Maps do not apply there."""
    if on_strands:
        return nodes.new("ShaderNodeHairInfo").outputs["Tangent Normal"]
    if "normal_tex" not in spec:
        return nodes.new("ShaderNodeNewGeometry").outputs["Normal"]
    texture = nodes.new("ShaderNodeTexImage")
    texture.image = load_packed_image(package, spec["normal_tex"], 'Non-Color')
    normal_map = nodes.new("ShaderNodeNormalMap")
    normal_map.uv_map = "UVMap"
    normal_map.inputs["Strength"].default_value = spec.get("normalmap_scale", 1)
    links.new(texture.outputs["Color"], normal_map.inputs["Color"])
    return normal_map.outputs[0]


def build_hair_shader(spec, nodes):
    """The renderer's hair lobe is Chiang 2016, as Cycles': the same colour mapping and roughness terms."""
    fibre = spec["fibre"]
    shader = nodes.new("ShaderNodeBsdfHairPrincipled")
    shader.model = 'CHIANG'
    shader.parametrization = 'COLOR'
    shader.inputs["Color"].default_value = (*fibre["color"], 1.0)
    shader.inputs["Roughness"].default_value = fibre["longitudinal_roughness"]
    shader.inputs["Radial Roughness"].default_value = fibre["azimuthal_roughness"]
    shader.inputs["Coat"].default_value = 0.0
    shader.inputs["IOR"].default_value = fibre["ior"]
    shader.inputs["Offset"].default_value = math.radians(fibre["cuticle_tilt_deg"])
    return shader


def build_fabric_shader(spec, package, nodes, links, normal, fabric_groups):
    sheen = spec["sheen"]
    highlight = has_highlight(spec)
    group_key = (sheen["roughness"], sheen.get("coverage", 1.0), spec.get("translucency", {}).get("transmit_share", 0.0),
                 highlight)
    if group_key not in fabric_groups:
        fabric_groups[group_key] = build_fabric_group(*group_key)
    shader = nodes.new("ShaderNodeGroup")
    shader.node_tree = fabric_groups[group_key]
    if highlight:
        connect_highlight(spec["specular"], package, nodes, links, shader)
    albedo = albedo_output(spec, package, nodes, links)
    links.new(albedo, shader.inputs["Albedo"])
    # C = color + tint (albedo - color). The colour sockets of a Mix node are inputs 6, 7 and output 2.
    sheen_color = nodes.new("ShaderNodeMix")
    sheen_color.data_type = 'RGBA'
    sheen_color.inputs[0].default_value = sheen["tint"]
    sheen_color.inputs[6].default_value = (*sheen["color"], 1.0)
    links.new(albedo, sheen_color.inputs[7])
    links.new(sheen_color.outputs[2], shader.inputs["Sheen Color"])
    links.new(normal, shader.inputs["Normal"])
    return shader


def build_coated_shader(spec, package, nodes, links, normal):
    coat = spec["coat"]
    shader = nodes.new("ShaderNodeBsdfPrincipled")
    links.new(albedo_output(spec, package, nodes, links), shader.inputs["Base Color"])
    shader.inputs["Specular IOR Level"].default_value = 0.0  # the diffuse base has no specular of its own
    shader.inputs["Coat Weight"].default_value = coat.get("coverage", 1.0)
    shader.inputs["Coat Roughness"].default_value = coat["roughness"]
    shader.inputs["Coat IOR"].default_value = coat["ior"]
    links.new(normal, shader.inputs["Normal"])
    links.new(normal, shader.inputs["Coat Normal"])
    return shader


def build_diffuse_shader(spec, package, nodes, links, normal):
    shader = nodes.new("ShaderNodeBsdfDiffuse")
    links.new(albedo_output(spec, package, nodes, links), shader.inputs["Color"])
    links.new(normal, shader.inputs["Normal"])
    return shader


def scaled_map(package, nodes, links, scalar, texture_path):
    """A scalar, times its map where there is one, clamped to [0, 1]: how the renderer folds a roughness or
    reflectance map into the authored value (Non-Color, as it reads them)."""
    value = nodes.new("ShaderNodeValue")
    value.outputs[0].default_value = scalar
    if not texture_path:
        return value.outputs[0]
    texture = nodes.new("ShaderNodeTexImage")
    texture.image = load_packed_image(package, texture_path, 'Non-Color')
    product = nodes.new("ShaderNodeMath")
    product.operation = 'MULTIPLY'
    product.use_clamp = True
    links.new(value.outputs[0], product.inputs[0])
    links.new(texture.outputs["Color"], product.inputs[1])
    return product.outputs[0]


def has_highlight(spec):
    """A specular block that is switched on (coverage 0 keeps the settings but draws nothing)."""
    return "specular" in spec and spec["specular"].get("coverage", 1.0) > 0.0


def check_highlight(spec):
    """Stops on a highlight the conversion does not model, so nothing is dropped silently. `color` is only the hue
    `specular_tint` moves towards, never a factor, so it matters only with a tint."""
    specular = spec["specular"]
    if set(specular) - HIGHLIGHT_KEYS:
        sys.exit(f"material {spec['id']}: unsupported highlight keys {sorted(set(specular) - HIGHLIGHT_KEYS)}")
    if specular.get("specular_tint", 0.0) != 0.0:
        sys.exit(f"material {spec['id']}: a tinted highlight (specular_tint) is not handled")
    if specular.get("coverage", 1.0) not in (0.0, 1.0):
        sys.exit(f"material {spec['id']}: a highlight coverage other than 0 or 1 is not handled")


def connect_highlight(specular, package, nodes, links, shader):
    """Feeds the highlight into the shader's sockets HIGHLIGHT_SOCKETS (a Principled BSDF or the fabric group): the
    roughness and reflectance maps scale the authored values (the renderer's roughness_tex scales both roughness
    axes; Specular IOR Level 0.5 is the plain dielectric F0 in Blender as in the renderer). The anisotropy is an
    approximation: the renderer gives the roughness along its tangent and across it (roughness_y), Blender one
    roughness and an anisotropy stretching the highlight along the tangent, so the tangent is turned a quarter turn
    from the renderer's rotation (half-turns; Blender's unit is a full turn). Both measure it from the u direction
    of the uv map, which the Tangent node gives (without it Cycles takes a radial tangent around the object's Z)."""
    links.new(scaled_map(package, nodes, links, specular["roughness"], specular.get("roughness_tex")), shader.inputs["Roughness"])
    links.new(scaled_map(package, nodes, links, specular["specular"], specular.get("reflectance_tex")), shader.inputs["Specular IOR Level"])
    roughness_y = specular.get("roughness_y", specular["roughness"])
    anisotropic, rotation = 0.0, 0.0
    if roughness_y > specular["roughness"]:
        anisotropic = 1.0 - specular["roughness"] / roughness_y
        rotation = (specular.get("aniso_rotation", 0.0) / 2.0 + 0.25) % 1.0
    shader.inputs["Anisotropic"].default_value = anisotropic
    shader.inputs["Anisotropic Rotation"].default_value = rotation
    tangent = nodes.new("ShaderNodeTangent")
    tangent.direction_type = 'UV_MAP'
    tangent.uv_map = "UVMap"
    links.new(tangent.outputs["Tangent"], shader.inputs["Tangent"])


def build_plastic_shader(spec, package, nodes, links, normal):
    """A yarn of fibres baked into maps (the package's plastic material): a Principled BSDF with the cavity map
    under the albedo, the normal map, and the highlight (connect_highlight). With coverage 0 it is the body alone."""
    if "specular" not in spec:
        sys.exit(f"material {spec['id']}: a plastic without its specular block")
    shader = nodes.new("ShaderNodeBsdfPrincipled")
    links.new(albedo_output(spec, package, nodes, links), shader.inputs["Base Color"])
    if has_highlight(spec):
        connect_highlight(spec["specular"], package, nodes, links, shader)
    else:
        shader.inputs["Specular IOR Level"].default_value = 0.0
    links.new(normal, shader.inputs["Normal"])
    return shader


def build_material(spec, package, fabric_groups, on_strands):
    known = {"id", "name", "scope", "schema_version", "provenance", "type", "diffuse", "sheen", "coat",
             "fibre", "normal_tex", "normalmap_scale", "translucency", "specular"}
    if set(spec) - known:
        sys.exit(f"material {spec['id']}: unsupported keys {sorted(set(spec) - known)}")
    if "translucency" in spec and (spec["type"] != "fabric" or set(spec["translucency"]) != {"transmit_share"}):
        sys.exit(f"material {spec['id']}: translucency is a fabric's transmit_share only")
    if spec["type"] not in ("hair", "fabric", "coated_diffuse", "diffuse", "plastic"):
        sys.exit(f"material {spec['id']}: unsupported type {spec['type']}")
    if "specular" in spec and spec["type"] not in ("plastic", "fabric"):
        sys.exit(f"material {spec['id']}: a specular lobe on a {spec['type']} material")
    if "specular" in spec:
        check_highlight(spec)
        if on_strands and has_highlight(spec):
            sys.exit(f"material {spec['id']}: a highlight on strands (they have no uv for its tangent)")
    if on_strands and "albedo_tex" in spec.get("diffuse", {}):
        sys.exit(f"material {spec['id']}: a texture on strands (they have no uv)")
    if spec["type"] == "hair" and not on_strands:
        sys.exit(f"material {spec['id']}: a hair material on a mesh")

    material = bpy.data.materials.new(spec["id"])
    material.use_nodes = True
    nodes, links = material.node_tree.nodes, material.node_tree.links
    nodes.clear()
    output = nodes.new("ShaderNodeOutputMaterial")
    if spec["type"] == "hair":
        shader = build_hair_shader(spec, nodes)
    else:
        normal = normal_output(spec, package, nodes, links, on_strands)
        if spec["type"] == "fabric":
            shader = build_fabric_shader(spec, package, nodes, links, normal, fabric_groups)
        elif spec["type"] == "coated_diffuse":
            shader = build_coated_shader(spec, package, nodes, links, normal)
        elif spec["type"] == "plastic":
            shader = build_plastic_shader(spec, package, nodes, links, normal)
        else:
            shader = build_diffuse_shader(spec, package, nodes, links, normal)
    links.new(shader.outputs[0], output.inputs["Surface"])
    return material


def build_collections(scene, scene_json):
    """One collection per package group, the project's as the root (parts without a group go there).
    An inner group repeating its parent's name holds the yarn."""
    root = bpy.data.collections.new(scene_json["project_name"])
    scene.collection.children.link(root)
    collections = {None: root}
    for group in scene_json.get("groups", []):
        if "parent_id" not in group:
            collections[group["id"]] = root
            continue
        parent = collections[group["parent_id"]]
        name = "yarn" if group["name"] == parent.name else group["name"]
        collection = bpy.data.collections.new(name)
        parent.children.link(collection)
        collections[group["id"]] = collection
    return collections


def build_objects(scene_json, package, collections):
    strands = scene_json.get("strands", [])
    assets = scene_json["meshes"] + strands
    strand_paths = {s["path"] for s in strands}
    specs = {m["id"]: m for m in scene_json["materials"]}
    fabric_groups = {}
    materials = {}
    seen_meshes = {}
    for instance in scene_json["instances"]:
        transform = instance["transform"]
        if (transform["translate"], transform["rotate_euler_deg"], transform["scale"]) != \
                ([0, 0, 0], [0, 0, 0], [1, 1, 1]):
            sys.exit(f"instance {instance['id']}: transform is not identity")
        part = str(instance["asset_mesh_index"])
        shapes = instance.get("part_fibre_shapes", {})
        if set(shapes) - {part}:
            sys.exit(f"instance {instance['id']}: part_fibre_shapes for a part it does not place")
        asset = assets[instance["asset_mesh_index"]]
        material_id = instance["material"]["id"]
        path = os.path.join(package, asset["path"])
        on_strands = asset["path"] in strand_paths
        key = (material_id, on_strands)
        if key not in materials:
            materials[key] = build_material(specs[material_id], package, fabric_groups, on_strands)

        if on_strands:
            obj = build_curves(asset["name"], path, materials[key], shapes.get(part))
        else:
            # A hair-only import's body is one mesh written once per chart: coincident copies add nothing.
            with open(path, "rb") as fh:
                digest = (hashlib.sha256(fh.read()).hexdigest(), material_id)
            if digest in seen_meshes:
                seen_meshes[digest].name = material_id
                seen_meshes[digest].data.name = material_id
                print(f"skipped {asset['name']}: identical to an earlier mesh")
                continue
            obj = build_mesh(asset["name"], path, materials[key])
            seen_meshes[digest] = obj
        collections[instance.get("group_id")].objects.link(obj)


def build_camera(scene, spec):
    """One package camera. Returns its object."""
    data = bpy.data.cameras.new(spec["name"])
    data.sensor_fit = 'VERTICAL'
    data.sensor_height = 24.0
    data.angle_y = math.radians(spec["fov_deg"])  # the renderer's fov_deg is vertical
    data.clip_start = 0.01
    data.clip_end = 100.0
    data.dof.use_dof = False
    z_up = lambda p: Vector((p[0], -p[2], p[1]))
    eye, target, up = z_up(spec["position"]), z_up(spec["look_at"]), z_up(spec["up"])
    forward = (target - eye).normalized()
    right = forward.cross(up).normalized()
    camera = bpy.data.objects.new(spec["name"], data)
    rotation = Matrix((right, right.cross(forward), -forward)).transposed()
    camera.matrix_world = Matrix.Translation(eye) @ rotation.to_4x4()
    scene.collection.objects.link(camera)
    return camera


def build_world(scene, env):
    if env["mode"] != "analytical":
        sys.exit(f"environment mode {env['mode']} is not handled")
    world = bpy.data.worlds.new("ppt_env")
    world.use_nodes = True
    background = world.node_tree.nodes["Background"]
    background.inputs["Color"].default_value = (*env["sky_color"], 1.0)
    background.inputs["Strength"].default_value = env["sky_intensity"]
    scene.world = world


def build_tone_curve(scene):
    """x (A x + B) / (x (C x + D) + E) per channel, clamped. The Standard view adds the sRGB OETF."""
    tree = bpy.data.node_groups.new("PPT tone curve", "CompositorNodeTree")
    tree.interface.new_socket("Image", in_out='OUTPUT', socket_type='NodeSocketColor')
    nodes, links = tree.nodes, tree.links

    def math_node(operation, a, b, clamp=False):
        node = nodes.new("ShaderNodeMath")
        node.operation = operation
        node.use_clamp = clamp
        for socket, value in zip(node.inputs, (a, b)):
            if isinstance(value, bpy.types.NodeSocket):
                links.new(value, socket)
            else:
                socket.default_value = value
        return node.outputs[0]

    layers = nodes.new("CompositorNodeRLayers")
    split = nodes.new("CompositorNodeSeparateColor")
    merge = nodes.new("CompositorNodeCombineColor")
    output = nodes.new("NodeGroupOutput")
    links.new(layers.outputs["Image"], split.inputs[0])
    for index in range(3):
        x = split.outputs[index]
        top = math_node('MULTIPLY', x, math_node('ADD', math_node('MULTIPLY', x, TONE_A), TONE_B))
        bottom = math_node('ADD', math_node('MULTIPLY', x, math_node(
            'ADD', math_node('MULTIPLY', x, TONE_C), TONE_D)), TONE_E)
        links.new(math_node('DIVIDE', top, bottom, clamp=True), merge.inputs[index])
    links.new(split.outputs[3], merge.inputs[3])
    links.new(merge.outputs[0], output.inputs["Image"])
    scene.compositing_node_group = tree
    scene.view_settings.view_transform = 'Standard'
    scene.view_settings.look = 'None'


def apply_render_settings(scene, spec):
    scene.unit_settings.system = 'METRIC'
    scene.unit_settings.scale_length = 1.0
    scene.render.engine = 'CYCLES'
    scene.render.resolution_x = spec["width"]
    scene.render.resolution_y = spec["height"]
    scene.render.resolution_percentage = 100
    scene.render.use_compositing = True
    scene.cycles_curves.shape = 'THICK_LINEAR'
    cycles = scene.cycles
    cycles.device = 'GPU'
    cycles.samples = spec["spp"]
    cycles.preview_samples = spec["spp"]
    cycles.use_adaptive_sampling = False
    cycles.use_preview_adaptive_sampling = False
    cycles.use_denoising = False
    cycles.use_preview_denoising = False
    cycles.pixel_filter_type = 'BOX'
    cycles.filter_width = 1.0
    cycles.max_bounces = spec["max_bounces"]
    cycles.diffuse_bounces = spec["max_bounces"]
    cycles.glossy_bounces = spec["max_bounces"]
    cycles.transmission_bounces = spec["max_bounces"]
    cycles.transparent_max_bounces = spec["max_bounces"]
    cycles.volume_bounces = 0
    cycles.sample_clamp_direct = 0.0
    cycles.sample_clamp_indirect = 0.0
    cycles.blur_glossy = 0.0


def prepare_viewports():
    """Open on the camera, and let a rendered viewport run the tone curve like F12 does."""
    for screen in bpy.data.screens:
        for area in screen.areas:
            for space in area.spaces:
                if space.type == 'VIEW_3D':
                    space.shading.use_compositor = 'ALWAYS'
                    space.region_3d.view_perspective = 'CAMERA'
                    space.clip_start = 0.01


def main():
    if bpy.app.version < (5, 1, 0):
        sys.exit(f"needs Blender 5.1 or later, this is {bpy.app.version_string}")
    args = sys.argv[sys.argv.index("--") + 1:]
    package, out_path = os.path.abspath(args[0]), os.path.abspath(args[1])
    with open(os.path.join(package, "scene.json"), encoding="utf-8") as fh:
        scene_json = json.load(fh)
    if scene_json.get("lights"):
        sys.exit("the package has area lights, which the Blender conversion does not handle")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    build_objects(scene_json, package, build_collections(scene, scene_json))
    cameras = {spec["id"]: build_camera(scene, spec) for spec in scene_json["cameras"]}
    default_id = scene_json["default_camera_id"]
    scene.camera = cameras[default_id]
    build_world(scene, scene_json["env"])
    apply_render_settings(scene, next(c for c in scene_json["cameras"] if c["id"] == default_id))
    build_tone_curve(scene)
    prepare_viewports()
    # A preset's credit: the .blend may travel without the file beside it.
    credits = os.path.join(package, "CREDITS.txt")
    if os.path.exists(credits):
        with open(credits, encoding="utf-8") as fh:
            block = bpy.data.texts.new("CREDITS.txt")
            block.write(fh.read())
            block.use_fake_user = True
    bpy.ops.wm.save_as_mainfile(filepath=out_path, compress=True)


main()
