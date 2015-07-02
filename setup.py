#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Python installation file.
"""
from setuptools import setup

REQUIREMENTS = ['numpy', 'scipy']

CLASSIFIERS = ['Development Status :: 4 - Beta',
               'Intended Audience :: Science/Research',
               'Natural Language :: English',
               'License :: OSI Approved :: Apache Software License',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Programming Language :: Python :: 2.7',
               ]

setup(name='agilegeo',
      version=open('version.txt').read().rstrip(),
      author='Agile Geoscience',
      author_email='hello@agilegeoscience.com',
      packages=['agilegeo'],
      description='Useful geophysics functions',
      long_description=open('README.rst').read(),
      url='http://pypi.python.org/pypi/agilegeo/',
      install_requires=REQUIREMENTS,
      classifiers=CLASSIFIERS,
      license='Apache 2',
      )
