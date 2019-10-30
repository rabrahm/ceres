from distutils.core import setup, Extension
import numpy
import os

#module = Extension('Marsh', sources = ['Marsh.c'],libraries=['gsl','gslcblas','m'], library_dirs=['/home/rabrahm/gsl/lib'], include_dirs=[numpy.get_include(),'/home/rabrahm/gsl/include'])

if os.access('../../gsl.temp',os.F_OK):
    path=open('../../gsl.temp','r').readlines()[0]
else:
    path = '/usr/local'

module =Extension('FunNorm', sources = ['FunNorm.c'],libraries=['gsl','gslcblas','m'],library_dirs=[path + '/lib'],include_dirs=[numpy.get_include()+'/numpy',path + '/include'])
setup(name = 'Funciones creadas por Nestor: Extensiones de C/Python', version = '1.0', ext_modules = [module])
