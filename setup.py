#!/usr/bin/env python

from distutils.core import setup


def get_version():
    """
    Read the file to get current version
    """
    
    for line in open("BGClib/BGClib.py", "r"):
        if line.startswith('__version__'):
            return line.split("=")[-1].strip().strip('"')
        
    return ""


setup(name='BGClib',
      python_requires='>=3.5',
      version=get_version(),
      description='Biosynthetic Gene Cluster library utilities',
      author='Jorge Navarro',
      author_email='j.navarro@westerdijkinstitute.nl',
      packages=['BGClib'],
      package_data={'BGClib': ['data/*.tsv', 'data/Domain_models/*']},
      keywords=['secondary metabolites', 'natural products', 'genome mining'],
     )
