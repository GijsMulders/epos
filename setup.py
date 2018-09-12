#!/usr/bin/env ipython
from setuptools import setup, find_packages

setup(name='epospy',
      version='1.0.2',
      description='the Exoplanet Population Observation Simulator',
      url='https://github.com/GijsMulders/epos',
      author='Gijs Mulders',
      author_email='gdmulders@gmail.com',
      license='MIT',
      #packages=['epos',find_packages()], # ??
      #packages=['epos','EPOS.plot'], # is this correct?
      packages=find_packages(),
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
