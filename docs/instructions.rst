Installation
============

Github
------

Download the source from `github <https://github.com/GijsMulders/epos>`_, preferably using git::

   git clone https://github.com//GijsMulders/epos

In the future you can update to the latest version with::

   git pull
   
Note that this may overwrite any changes you have made to the example scripts.

Configuration
-------------

You can run EPOS from the main directory (``epos/``) where the test and example scripts are located. 

If you want to run from a different directory, you can make a soft-link to the ``EPOS/`` and ``files/`` folder located within the main ``epos/`` dir::

   cd run_epos_from_here/
   ln -s path/to/epos/EPOS EPOS
   ln -s path/to/epos/files files

Dependencies
------------

EPOS has the following dependencies and will probably throw an error if you have not installed them:

* python 2.7
* numpy 1.13+
* scipy
* matplotlib 2.0+

The following software is required if you want to use the MCMC part of EPOS 

* `emcee <http://dan.iel.fm/emcee>`_ for the MCMC fitting
* `corner.py <(http://corner.readthedocs.io/>`_ for the corner plots

Testing
-------
To test the basic functionality of EPOS, you can run the test scripts in the main directory::

   ./test_1_survey.py
   ./test_2_montecarlo.py
   ./test_3_mcmc.py
   ./test_4_multicore.py
   ./test_5_occurrence.py

If you don't have ipython installed, you can also run ``python test_1_survey.py`` or copy-paste the script into a terminal.

Each command is described in the comments, and provides an introduction to the basic functionality. 

Examples
========

There are a number of scripts with examples in the root ``epos`` folder that demonstrate different functionality. The scripts contains comments describing each command.

Parametric mode
---------------
The first example shows how to fit a parametric function in planet radius and orbital period to Kepler data::

   ./example_1_parametric_mode.py




Multi-planet Mode
-----------------
::

   ./example_2_multiplanet_mode.py

Planet Ocurrence Rates
----------------------

Example 9 shows how to estimate occurrence rates using the inverse detection efficiency method. You can run the entire script with:
:: 

   ./example_9_occurrence_rate_inverse.py

Here, is a step-by-step description of each set of commands.
First, load EPOS and set the output directory name
::

   import EPOS
   epos= EPOS.epos(name='example_9')

Second, load the observations (Kepler DR25 planet candidate list) and survey dectetion efficiency
::

   obs, survey= EPOS.kepler.dr25(Huber=True, Vetting=True, score=0.9)
   epos.set_observation(**obs)
   epos.set_survey(**survey)

Next, define the occurrence rate bins for hot Jupiters and super-earths/mini-Neptunes:
::

   x_HJ= [1,10] # Orbital period range in days
   y_HJ= [7,20] # Planet size range in earth radii
   x_SEMN, y_SEMN= [2,150],[1.0,4.0] # super-Earths/mini-Neptunes
   epos.set_bins(xbins=[x_HJ, x_SEMN], ybins=[y_HJ, y_SEMN])

The rates are then calculated, plotted, and saved
::

   EPOS.occurrence.all(epos)
   EPOS.save.occurrence(epos)
   EPOS.plot.occurrence.all(epos)

The output appears in ``png/occurrence/bins.png`` and should look like this:

.. image:: fig_example_9.png

Alternatively, you can generate a 1D or 2D grid of bins, for example the SAG13 grid:
::

   import numpy as np
   epos.set_bins(xgrid=np.geomspace(10,640,7), 
   	ybins=np.geomspace(0.67,17,9))

.. image:: fig_example_9_SAG13.png

FAQ
===

Frequently asked questions
--------------------------

If you have any difficulties or questions running EPOS that are not addressed in the documentation or FAQ please contact gdmulders@gmail.com

I'm getting an AttributeError: 'module' object has no attribute 'geomspace'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please upgrade to numpy 1.13 or a more recent version

