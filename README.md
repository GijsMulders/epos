# README #

Welcome to the Exoplanet Population Observation Simulator!

## What is this repository for? ##

This repository hosts the EPOS source code

## How do I get set up? ##

### Summary of set up ###
Download the source, preferably using git
```
git clone https://username@bitbucket.org/GijsMulders/epos.git
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
Run the examples that test basic functionality

```
./example1.py
./example2.py
```

if you don't have ipython installed, you can also run `python example1.py`

### Where can I find documentation? ###

A draft paper describing the code can be [dowloaded here](https://www.dropbox.com/s/964mwknjdcueyj9/EPOS-draft.pdf?dl=1)

Different modules are described in _build/html/index.html

### Who do I talk to? ###

* Gijs Mulders: gdmulders@gmail.com
* Slack: eposua.slack.com