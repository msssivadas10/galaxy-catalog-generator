from setuptools import setup
from setuptools.command.build_ext import build_ext
import os, os.path, sys, subprocess

PACKAGE_NAME = "haloutils"
SRC_DIR      = "src"
BUILD_DIR    = "build_fortran"

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
        files = [ os.path.join(SRC_DIR, *p) for p in files ]
        subprocess.check_call([ fc, *args, *files, "-o", lib_path ])

        os.remove(BUILD_DIR)
        return
    
    def run(self):
        ...

# setup(
#     name="mymod",
#     version="0.1",
#     cmdclass={"build_ext": BuildShared},
# )

files = {
    ("utils",  "random"        ): [],            
    ("utils",  "constants"     ): [],               
    ("utils",  "integrate"     ): [],               
    ("utils",  "interpolate"   ): [],                 
    (          "zfunctions"    ): [],                
    (          "rfunctions"    ): [],                
    (          "halo_model"    ): [],                
    ("models", "power_spectrum"): [],                    
    ("models", "mass_function" ): [],                   
    ("models", "halo_bias"     ): [],               
    (          "galaxy_catalog"): [],                    
}