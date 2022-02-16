#!/usr/bin/python
#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy as np

cython_cflags = ["-Wno-unused"] # NUMPY EXTENSION MADNESS
cython_include_dirs = ['build/include']

setup_requires = [
    'cython'
]

extensions = cythonize([Extension("eikonal.solver", ['build/solver.pyx'],
                                  library_dirs=[],
                                  language="c++",
                                  include_dirs=cython_include_dirs,
                                  extra_compile_args=cython_cflags),
                        Extension("eikonal.raytrace", ['build/raytrace.pyx'],
                                  library_dirs=[],
                                  language="c++",
                                  include_dirs=cython_include_dirs,
                                  extra_compile_args=cython_cflags),
                        ])

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=extensions,
    packages=['eikonal', 'eikonal.vis'],
    name='Eikonal',
    version = '0.0.1',
    include_dirs=[np.get_include()],
    description="",
    author="J-Pascal Mercier",
    author_email="jean-pascal.mercier@agsis.com",
    license="All rights reserved."
)

