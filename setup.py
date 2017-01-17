#!/usr/bin/env python
import os
import sys
import pdb4amber

try:
    if '--no-setuptools' in sys.argv:
        sys.argv.remove('--no-setuptools')
        raise ImportError() # Don't import setuptools...
    from setuptools import setup
    kws = dict(entry_points={
        'console_scripts' : ['pdb4amber = pdb4amber.pdb4amber:main'],
    })
except ImportError:
    from distutils.core import setup
    kws = {'scripts' : [os.path.join('pdb4amber', 'pdb4amber'),]
    }

setup(name='pdb4amber',
      version=pdb4amber.__version__,
      description='PDB analyzer to prepare PDB files for Amber simulations.',
      author='Romain M. Wolf and AMBER developers',
      author_email='amber@ambermd.org',
      url='http://ambermd.org/',
      packages=['pdb4amber'],
      **kws
)
