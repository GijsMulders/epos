# README #

Welcome to the Exoplanet Population Observation Simulator!
The documentation will be hosted [here](http://epos.readthedocs.io/en/latest)

## What is this repository for? ##

This repository hosts the EPOS source code

## How do I get set up? ##

### Summary of set up ###
Download the source, preferably using git
```
git clone https://github.com//GijsMulders/epos
```
Then you can update to the latest version with
```
git pull
```

### Configuration ###
You can run EPOS from the main directory where the test and example scripts are located. 

If you want to run from a different directory, you can make a soft-link to the EPOS/ and files/ folder located within the main epos dir
```
cd run_epos_from_here/
ln -s path_to_epos/EPOS EPOS
ln -s path_to_epos/files files
```

### Dependencies ###
* python 2.7
* numpy 1.13+
* scipy
* matplotlib 2.0+
* [emcee](http://dan.iel.fm/emcee) for the MCMC fitting
* [corner.py](http://corner.readthedocs.io/) for the MCMC plots

### How to run tests ###
The test scripts test some basic functionality

```
./test_1_survey.py
./test_2_montecarlo.py
./test_3_mcmc.py
./test_4_multicore.py
./test_5_occurrence.py
```

if you don't have ipython installed, you can also run `python test_1_survey.py`

The example scripts describe some additional functionality:

```
./example_1_parametric_mode.py
./example_2_multiplanet_mode.py
```

### Where can I find documentation? ###

EPOS is described in [Mulders et al. 2018](https://arxiv.org/abs/1805.08211)
The documentation is hosted at [readthedocs](http://epos.readthedocs.io/en/latest)

### Who do I talk to? ###

* Gijs Mulders: gdmulders@gmail.com
