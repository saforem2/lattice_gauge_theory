# -*- coding: utf-8 -*-

# Learn more: https://github.com/saforem2/lattice_gauge_theory

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='lattice_gauge_theory',
    version='0.1.0',
    description='Monte Carlo simulation of Z(N) models in lattice gauge theory.',
    long_description=readme,
    author='Sam Foreman',
    author_email='samuel-foreman@uiowa.edu',
    url='https://github.com/saforem2/lattice_gauge_theory',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)


