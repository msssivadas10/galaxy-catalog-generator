from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os, os.path, subprocess

class BuildShared(build_ext):
    def run(self):
        # compile Fortran to shared lib manually
        files = [
            os.path.join("src", "utils", "integrate.f90"),
            os.path.join("src", "power", "growthfactor.f90"),
            os.path.join("src", "power", "correlation.f90"),
            os.path.join("src", "power", "variance.f90"),
            os.path.join("src", "power", "models", "ingredients.f90"),
            os.path.join("src", "power", "models", "eisenstein98_zb.f90"),
            os.path.join("src", "power", "models", "eisenstein98_mnu.f90"),
            os.path.join("src", "power", "models", "eisenstein98_bao.f90"),
        ]
        args = [ "-shared", "-fPIC", "-J", ".include" ]
        subprocess.check_call(
            ["gfortran", *args, *files, "-o", "libpowerspectrum.so"]
        )

setup(
    name="mymod",
    version="0.1",
    cmdclass={"build_ext": BuildShared},
)
