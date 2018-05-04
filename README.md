# README #

Welcome to the Exoplanet Population Observation Simulator!
The documentation will be hosted [here](http://epos.readthedocs.io/en/latest)

## What is this repository for? ##

This repository hosts the EPOS source code

## How do I get set up? ##

### Summary of set up ###
Download the source, preferably using git
```
git clone https://github.com//GijsMulders/epos.git
```
Then you can update to the latest version with
```
git pull
```

### Configuration ###
None. Optionally add the epos/EPOS dir to your python path

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

A draft paper describing the code can be [dowloaded here](https://www.dropbox.com/s/964mwknjdcueyj9/EPOS-draft.pdf?dl=1)

### Who do I talk to? ###

* Gijs Mulders: gdmulders@gmail.com
* Slack: eposua.slack.com
