from distutils.core import setup, Extension
import numpy

module =Extension('FunNorm', sources = ['FunNorm.c'],libraries=['gsl','gslcblas','m'],include_dirs=[numpy.get_include()+'/numpy'])
setup(name = 'Funciones creadas por Nestor: Extensiones de C/Python', version = '1.0', ext_modules = [module])
