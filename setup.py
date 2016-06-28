from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Orbit compute',
  ext_modules = cythonize("c_orbits.pyx"),
)



