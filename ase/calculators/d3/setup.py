import os
from os.path import join

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils import fcompiler

    if fcompiler.find_executable('gfortran') is not None:
        extra_compile_args = '-fopenmp'
        extra_link_args = '-lgomp'
        fc = 'gfortran'
    elif fcompiler.find_executable('ifort') is not None:
        extra_compile_args = '-openmp'
        extra_link_args = '-lmkl_rt -lpthread -lm -liomp5'
        fc = 'ifort'
    else:
        return

    config = Configuration('d3', parent_package, top_path)

    config.add_extension('d3ef',
            sources=['d3ef.pyf','d3ef.f90'],
            extra_link_args=[extra_link_args],
            extra_f90_compile_args=[extra_compile_args],
#            f2py_options=['--fcompiler=' + fc],
#            extra_info={'fcompiler': fc},
            )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
