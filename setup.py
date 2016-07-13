from distutils.core import setup, Extension
#from Cython.Build import cythonize



setup(
  name = 'PNauty',
  ext_modules=[Extension('pnauty', ['pnauty.c'])],
)

#setup(
#  name = 'Orbit compute',
#  #ext_modules = cythonize("orbits_wrap.pyx"),
#  ext_modules=[Extension('c_orbits', ['c_orbits.c'])],
#  extra_compile_args = ["-O0"], 
#)



