"""Which files a header change recompiles, per component (src/ARCHITECTURE.md).

Reads the component declarations (build/components.json, written by cmake/components.cmake) and the
header dependencies Ninja recorded in the last build (`ninja -t deps`), and prints:
  1. per component X: the translation units outside X that include a header of X, by including
     component, split into X's public and private headers (private must be 0 outside X's own tests);
  2. the touch check: for X's public and private header with the most dependents (and every --header),
     the translation units Ninja recompiles when that header changes;
  3. the static graph check: direct includes of a component that is not a declared dependency, P edges
     no public header uses, and "../" or absolute include spellings in src/.
Without components.json (a tree before the components) every folder of src/ counts as a component and
every header as public.

Usage: python -I cmake/component_report.py <build dir> [--header src/<c>/<path>.h ...] [--verify]
  --verify  builds, then touches each --header, builds again and compares what Ninja compiled with the
            prediction. Run it where the compiler is set up (the MSVC developer environment).
Needs ninja on PATH.
"""
import collections
import json
import os
import re
import subprocess
import sys

INCLUDE = re.compile(r'^\s*#\s*include\s+"([^"]+)"')


def norm(path):
    return path.replace("\\", "/").lower()


def load_components(build, source):
    path = os.path.join(build, "components.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        return {c["name"]: c for c in data["components"]}
    names = sorted(d for d in os.listdir(os.path.join(source, "src")) if os.path.isdir(os.path.join(source, "src", d)))
    return {n: {"name": n, "public_deps": [], "private_deps": [], "public_headers": None, "private_headers": []} for n in names}


def read_deps(build):
    text = subprocess.run(["ninja", "-C", build, "-t", "deps"], capture_output=True, text=True, check=True).stdout
    deps, current = {}, None
    for line in text.splitlines():
        if line and not line.startswith(" ") and ": #deps" in line:
            current = line.split(": #deps")[0]
            deps[current] = []
        elif line.strip() and current:
            deps[current].append(os.path.normpath(os.path.join(build, line.strip())))
    return deps


class Classifier:
    def __init__(self, components, source, build):
        self.components = components
        self.src = norm(os.path.join(source, "src")) + "/"
        self.tests = norm(os.path.join(source, "tests")) + "/"

    def header(self, path):
        """(component, 'public'|'private', path in the component) for a header under src/, else None."""
        p = norm(path)
        if not p.startswith(self.src):
            return None
        parts = p[len(self.src):].split("/", 1)
        if len(parts) < 2 or parts[0] not in self.components:
            return None
        c = self.components[parts[0]]
        if c["public_headers"] is None:
            return parts[0], "public", parts[1]
        return parts[0], ("private" if parts[1] in [h.lower() for h in c["private_headers"]] else "public"), parts[1]

    def owner(self, obj):
        """The component (or 'tests') a translation unit belongs to, from its object path."""
        o = norm(obj)
        m = re.search(r"cmakefiles/([^/]+)\.dir/(.*)$", o)
        if not m:
            return None
        target, rest = m.groups()
        for name, c in self.components.items():
            if target == c.get("target", "knit_" + name).lower() or target == "knit_" + name:
                return name
        if target.startswith("knit_tests") or "tests/" in o.split("cmakefiles")[0]:
            return "tests"
        first = rest.split("/")[0]
        return first if first in self.components else None


def static_check(components, source):
    problems = []
    for name, c in components.items():
        public = set(c["public_headers"] or [])
        declared = set(c["public_deps"]) | set(c["private_deps"]) | {name}
        used_by_public = set()
        for root, _, files in os.walk(os.path.join(source, "src", name)):
            for f in files:
                if not f.endswith((".h", ".cuh", ".cpp", ".cu")):
                    continue
                full = os.path.join(root, f)
                rel = norm(os.path.relpath(full, os.path.join(source, "src", name)))
                with open(full, encoding="utf-8", errors="replace") as fh:
                    for number, line in enumerate(fh, 1):
                        m = INCLUDE.match(line)
                        if not m:
                            continue
                        spelled = m.group(1)
                        if spelled.startswith("..") or ":" in spelled or spelled.startswith("/"):
                            problems.append(f"{name}/{rel}:{number}: include \"{spelled}\" is relative or absolute")
                            continue
                        target = spelled.split("/")[0]
                        if target not in components or "/" not in spelled:
                            continue
                        if c["public_headers"] is not None and target not in declared:
                            problems.append(f"{name}/{rel}:{number}: includes {target}, not a declared dependency")
                        if rel in [h.lower() for h in public]:
                            used_by_public.add(target)
        for dep in c["public_deps"]:
            if dep not in used_by_public:
                problems.append(f"{name}: PUBLIC dependency {dep} is included by no public header (make it PRIVATE)")
    return problems


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    build = args[0]
    extra = [args[i + 1] for i, a in enumerate(args) if a == "--header"]
    verify = "--verify" in args
    with open(os.path.join(build, "CMakeCache.txt"), encoding="utf-8") as f:
        source = re.search(r"CMAKE_HOME_DIRECTORY:INTERNAL=(.*)", f.read()).group(1).strip()
    components = load_components(build, source)
    cls = Classifier(components, source, build)
    deps = read_deps(build)

    users = collections.defaultdict(set)  # (component, vis, path) -> objects including it
    owners = {}
    for obj, files in deps.items():
        owner = cls.owner(obj)
        if owner is None:
            continue
        owners[obj] = owner
        for f in files:
            h = cls.header(f)
            if h:
                users[h].add(obj)

    order = list(components) + ["tests"]
    print("1. Translation units outside X that include a header of X (public / private), by including component")
    print("X".ljust(10) + "".join(o[:8].rjust(10) for o in order))
    for x in components:
        row = x.ljust(10)
        for o in order:
            pub = {obj for (c, v, _), objs in users.items() if c == x and v == "public" for obj in objs if owners[obj] == o}
            priv = {obj for (c, v, _), objs in users.items() if c == x and v == "private" for obj in objs if owners[obj] == o}
            row += ("-" if o == x else f"{len(pub)}/{len(priv)}").rjust(10)
        print(row)
    total = collections.Counter(owners.values())
    print("TUs".ljust(10) + "".join(str(total[o]).rjust(10) for o in order))

    print("\n2. Touch check: TUs recompiled when the header changes, by component")
    chosen = []
    for x in components:
        for vis in ("public", "private"):
            candidates = [(len(objs), h) for h, objs in users.items() if h[0] == x and h[1] == vis]
            if candidates:
                chosen.append(max(candidates)[1])
    for e in extra:
        h = cls.header(os.path.join(source, e))
        if h:
            chosen.append(h)
    for h in chosen:
        by = collections.Counter(owners[o] for o in users[h])
        print(f"  {h[0]}/{h[2]} ({h[1]}): {sum(by.values())} TUs: " + ", ".join(f"{k} {v}" for k, v in sorted(by.items())))

    print("\n3. Static graph check")
    problems = static_check(components, source)
    print("\n".join("  " + p for p in problems) if problems else "  none")

    if verify:
        # A real build: `ninja -n` stops at the regeneration check the CONFIGURE_DEPENDS glob adds to every run.
        print("\n4. Verify with a build")
        subprocess.run(["ninja", "-C", build], capture_output=True, text=True)  # up to date first
        for e in extra:
            path = os.path.join(source, e)
            h = cls.header(path)
            predicted = {norm(o) for o in users.get(h, set())}
            os.utime(path, None)
            out = subprocess.run(["ninja", "-C", build], capture_output=True, text=True).stdout
            rebuilt = {norm(m) for m in re.findall(r"Building \w+ object (\S+)", out)}
            ok = rebuilt == predicted
            print(f"  {e}: ninja rebuilds {len(rebuilt)}, predicted {len(predicted)}: {'same' if ok else 'DIFFERENT'}")
            for o in sorted(rebuilt ^ predicted):
                print(f"    {'only ninja' if o in rebuilt else 'only predicted'}: {o}")


if __name__ == "__main__":
    main()
