#!/usr/bin/env ipython
from setuptools import setup, find_packages

#exclude=['EPOS/scriptdir/examples/example_{}_*.py'.format(k) for k in [3,4,5,6,7,8]]
#print exclude

exclude= [
  'EPOS/scriptdir/examples/example_3_population_synthesis.py'
  'EPOS/scriptdir/examples/example_4_massradius.py'
  'EPOS/scriptdir/examples/example_5_radial_velocity.py'
  'EPOS/scriptdir/examples/example_6_q16_catalog.py'
  'EPOS/scriptdir/examples/example_7_spectral_type.py'
  'EPOS/scriptdir/examples/example_8_no_monte_carlo.py'
  ] # also in MANIFEST.in

setup(name='epospy',
      version='1.0.3',
      description='the Exoplanet Population Observation Simulator',
      url='https://github.com/GijsMulders/epos',
      author='Gijs Mulders',
      author_email='gdmulders@gmail.com',
      license='MIT',
      #packages=['epos',find_packages()], # ??
      #packages=['epos','EPOS.plot'], # is this correct?
      packages=find_packages(exclude=exclude),
      install_requires=[
          'pytest >= 2.8',
          'numpy >= 1.13',
          'scipy',
          'matplotlib >= 2.0',
          'astropy',
          'emcee','corner'
      ],
      classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
      ],
      include_package_data=True,
#      zip_safe=False
)
