# Blender file export script (GPL-3.0-or-later)

`ppt_package_to_blend.py` turns an exported package (scene.json, PLY meshes, `.pptstrands` fibres) into a `.blend`. The
viewer's Export step runs it for the Blender format, in a Blender 5.1 or later that the user installed:

    blender -b --factory-startup --python blender/ppt_package_to_blend.py -- <package> <out.blend>

The script uses Blender's Python API (`bpy`), which Blender publishes under the GNU General Public License, so it is
released under the GNU General Public License, version 3 or later (`LICENSE` in this folder). It is the only part of Wool
Shader under that licence. The program itself is MIT (`../LICENSE`): it starts Blender as a separate program and passes
it this script, the package folder and the output path. Blender itself is not compiled into it or shipped with it;
the only files from Blender's source in Wool Shader are the MPL-2.0 SLIM files in `../third_party/slim/src`.
Without Blender the Blender format is greyed out, and everything else works.

Blender is a registered trademark of the Blender Foundation. Wool Shader is not made, endorsed or supported by the
Blender Foundation.
