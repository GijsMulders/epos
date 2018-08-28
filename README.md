# README #

Welcome to the Exoplanet Population Observation Simulator!
The documentation will be hosted [here](http://epos.readthedocs.io/en/latest)

## What is this repository for? ##

This repository hosts the EPOS source code

## How do I get set up? ##

### Summary of set up ###
You can install the latest version of epos with:
```
pip install epospy
```
And you should be ready to go!

Alternatively, download the source, preferably using git
```
git clone https://github.com//GijsMulders/epos
```
Then you can update to the latest version with
```
git pull
```
You can run EPOS from the main directory or add the EPOS/ directory to your path.

### Configuration ###
As of version 1.0.2, there is no configuration if you installed via pip and dependencies should be automatically resolved.

### Dependencies ###
* python 2.7
* numpy 1.13+
* scipy
* matplotlib 2.0+
* [emcee](http://dan.iel.fm/emcee) for the MCMC fitting
* [corner.py](http://corner.readthedocs.io/) for the MCMC plots
* astropy

### How do I learn how to use EPOS? ###
The best way is to run the test and example scripts as described in the [documentation](http://epos.readthedocs.io/en/latest/instructions.html#testing)

### Where can I find documentation? ###

EPOS is described in [Mulders et al. 2018](https://arxiv.org/abs/1805.08211)
The documentation is hosted at [readthedocs](http://epos.readthedocs.io/en/latest)

### Who do I talk to? ###

* Gijs Mulders: gdmulders@gmail.com
