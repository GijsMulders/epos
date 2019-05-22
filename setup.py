#!/usr/bin/env ipython
from setuptools import setup, find_packages
from EPOS import __version__

exclude= [
  'EPOS/scriptdir/examples/example_4_massradius.py'
  'EPOS/scriptdir/examples/example_6_q16_catalog.py'
  'EPOS/scriptdir/examples/example_7_spectral_type.py'
  'EPOS/scriptdir/examples/example_8_no_monte_carlo.py'
  ] # also in MANIFEST.in

setup(name='epospy',
      version=__version__,
      description='the Exoplanet Population Observation Simulator',
      url='https://github.com/GijsMulders/epos',
      author='Gijs Mulders',
      author_email='gdmulders@gmail.com',
      license='MIT',
      packages=find_packages(exclude=exclude),
      install_requires=[
          'pytest >= 2.8',
          'numpy >= 1.13',
          'scipy',
          'matplotlib >= 2.0',
          'astropy',
          'emcee >= 2.0',
          'corner >= 2.0',
          'h5py',
	       'shapely'
      ],
      classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
      ],
      include_package_data=True,
#      zip_safe=False
)
