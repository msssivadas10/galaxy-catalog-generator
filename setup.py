from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os, os.path, subprocess

class BuildShared(build_ext):
    def run(self):
        # compile Fortran to shared lib manually
        files = [
            os.path.join("src", "utils", "integrate.f90"),
            os.path.join("src", "growthfactor.f90"),
            os.path.join("src", "correlation.f90" ),
            os.path.join("src", "variance.f90"    ),
            os.path.join("src", "models", "power_spectrum.f90"),
            os.path.join("src", "models", "mass_function.f90" ),
            os.path.join("src", "models", "halo_bias.f90"     ),
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
