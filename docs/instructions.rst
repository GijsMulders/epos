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

Examples
========

Occurrence rate mode
--------------------

Multi-planet Mode
-----------------

FAQ
===

Frequently asked questions
--------------------------

If you have any difficulties or questions running EPOS that are not addressed in the documentation or FAQ please contact gdmulders@gmail.com

I'm getting an AttributeError: 'module' object has no attribute 'geomspace'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Please upgrade to numpy 1.13 or a more recent version

