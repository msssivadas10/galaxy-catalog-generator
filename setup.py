from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess

class BuildShared(build_ext):
    def run(self):
        # compile Fortran to shared lib manually
        subprocess.check_call(
            ["gfortran", "-shared", "-fPIC", "-J", ".include", 
             "src/quadutils.f90", 
             "src/powerspectrum.f90", 
             "src/power_integrals.f90",
             "src/growthfactor.f90", 
             "-o", "libpowerspectrum.so"]
        )

setup(
    name="mymod",
    version="0.1",
    cmdclass={"build_ext": BuildShared},
)
