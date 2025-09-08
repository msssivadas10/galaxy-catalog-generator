from setuptools import setup
from setuptools.command.build_ext import build_ext
import os, os.path, sys, subprocess, re, shutil

PACKAGE_NAME = "haloutils"
SRC_DIR      = "src"
BUILD_DIR    = "build_fortran"

def parse_fortran_file(filepath):
    # Extract the module name and used modules from a Fortran 90 file.

    # Regex patterns (case-insensitive)
    module_re  = re.compile(r'^\s*module\s+([a-z_][a-z0-9_]*)', re.IGNORECASE)
    use_re     = re.compile(r'^\s*use\s+([a-z_][a-z0-9_]*)'   , re.IGNORECASE)
    modproc_re = re.compile(r'^\s*module\s+procedure'         , re.IGNORECASE)

    module_name  = None
    used_modules = []
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.split("!")[0].strip() # Strip inline comments
            if not line: continue
            if modproc_re.match(line): continue # Skip "module procedure"

            m = module_re.match(line) # module denfinition
            if m and module_name is None:
                module_name = m.group(1).lower()
                continue

            u = use_re.match(line) # use statement
            if u: used_modules.append(u.group(1).lower())

    return module_name, used_modules

def scan_fortran_sources(root_dir, extensions=(".f90", ".F90")):
    # Return a dict mapping file paths to (module_name, used_modules).
    result = {}
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.endswith(extensions):
                fpath     = os.path.join(dirpath, fname)
                mod, uses = parse_fortran_file(fpath)
                if mod: result[fpath] = (mod, uses)
    return result

def build_dependency_graph(root_dir):
    # Build a dependency graph mapping file -> (module, [dependent files]).
    mapping = scan_fortran_sources(root_dir)
    module_to_file = {mod: path for path, (mod, _) in mapping.items()} # Reverse map: module -> file

    # Build resolved dependency graph
    graph = {}
    for path, (mod, uses) in mapping.items():
        resolved = []
        for u in uses:
            if u in module_to_file:
                resolved.append(module_to_file[u])  # link to file
        graph[path] = (mod, resolved)

    return graph

def topo_sort(graph):
    # Perform topological sort on a dependency graph.
    visited = {}   # file -> "temporary" or "permanent"
    order   = []
    
    def dfs(node):
        if node in visited:
            if visited[node] == "temp":
                raise ValueError(f"Cycle detected involving {node}")
            return
        visited[node] = "temp"
        _, deps = graph[node]
        for dep in deps:
            if dep.startswith("<unresolved:"): continue # ignore unresolved modules
            if dep not in graph: continue  # external file, ignore
            dfs(dep)
        visited[node] = "perm"
        order.append(node)

    for node in graph:
        if node not in visited:
            dfs(node)

    return order

class BuildShared(build_ext):

    def compile_as_sharedlib(self, files, lib_filename):

        if not os.path.exists(BUILD_DIR): os.mkdir(BUILD_DIR)

        # Detect platform-specific shared library extension
        if sys.platform.startswith("linux"):
            lib_filename = f"{lib_filename}.so"
        elif sys.platform == "darwin":
            lib_filename = f"{lib_filename}.dylib"
        elif sys.platform == "win32":
            lib_filename = f"{lib_filename}.dll"
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}")
        lib_path = os.path.join(PACKAGE_NAME, lib_filename)

        # Compile files to shared library
        fc   = os.environ.get("FC", "gfortran")
        args = ( 
            os.environ.get("FCARGS", "").split() 
            or 
            [ "-shared", "-fPIC", "-J", BUILD_DIR, "-fopenmp" ] 
        )
        print(f"using fortran compiler {fc} with args {args}"    , flush = True)
        print(f"compiling files {', '.join(files)} to {lib_path}", flush = True)
        subprocess.check_call([ fc, *args, *files, "-o", lib_path ])

        if os.path.exists(BUILD_DIR): shutil.rmtree(BUILD_DIR)
        return
    
    def run(self):

        # Shared lib for halo model
        path, lib_filename = os.path.join(SRC_DIR, "halo_model"), "libhm"
        files = topo_sort( build_dependency_graph(path) )
        self.compile_as_sharedlib(files, lib_filename)

        # Shared lib for estimators
        path, lib_filename = os.path.join(SRC_DIR, "estimators"), "libstm"
        files = topo_sort( build_dependency_graph(path) )
        self.compile_as_sharedlib(files, lib_filename)

        return

setup(
    name = PACKAGE_NAME,
    version = "0.1",
    packages = [ PACKAGE_NAME ],
    include_package_data = True,
    package_data = { PACKAGE_NAME: ["*.so", "*.dll", "*.dylib"] },
    cmdclass = {
        "build_ext": BuildShared,
    },
)