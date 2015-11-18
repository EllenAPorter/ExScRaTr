#! /usr/bin/env python3

from distutils.core import setup

setup(name='lmprop',
      version='1.0',
      description='Light Mesh Propagation Visualizer',
      author='Bob Lewis',
      author_email='bobl@tricity.wsu.edu',
      #url='https://www.python.org/sigs/distutils-sig/',

      packages=[ 'lmprop' ],
      scripts=[ 'scripts/lmprop' ],
      data_files= [ ('share/lmprop', [ 'images/eye_rev.png' ] ) ],
)
