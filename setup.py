#!/usr/bin/env python

"""Setup script for the identifiability module distribution."""

from setuptools import setup

__version__ = '0.3.2'

# The directory containing this file
README_name = __file__.replace('setup.py', 'README.md')

# The text of the README file
with open(README_name) as README:
    identifiability_main_doc = README.read()

# Requirements
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

modules = ['identifiability']

setup(  # Distribution meta-data
    name='identifiability',
    version=__version__,
    description='Parameter identifiability analysis in Python',
    long_description_content_type="text/markdown",
    long_description=identifiability_main_doc,
    author='Johann M. Rohwer',
    author_email='j.m.rohwer@gmail.com',
    url='https://github.com/jmrohwer/identifiability/',
    download_url='https://github.com/jmrohwer/identifiability/archive/refs/tags/v0.1.zip',
    license='New BSD',
    py_modules=modules,
    python_requires='>=3.6',
    install_requires=requirements,
    platforms=['Windows', 'Linux', 'macOS'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
