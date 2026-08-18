"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='scipion-chem-deeploc',  # Required
    version='2.1',  # Required
    description='Scipion plugin in order to use the DeepLoc software',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-chem/scipion-chem-deeploc',  # Optional
    author='Blanca Pueche',  # Optional
    author_email='blanca.pueche@cnb.csic.es',  # Optional
    keywords='scipion DeepLoc scipion-3.0 cheminformatics',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={'pyworkflow.plugin': 'deeploc = deeploc'},
    package_data={  # Optional
       'deeploc': [ 'protocols.conf'],
    }
)
