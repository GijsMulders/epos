.. epos documentation master file, created by
   sphinx-quickstart on Thu May  3 14:58:08 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EPOS
==============================================

.. image:: https://img.shields.io/badge/ascl-1909.013-blue.svg?colorB=262255
   :target: http://ascl.net/1909.013

.. image:: https://zenodo.org/badge/132065382.svg
   :target: https://zenodo.org/badge/latestdoi/132065382
   
the Exoplanet Population Observation Simulator
----------------------------------------------

EPOS is a software package to simulate observations of exoplanet populations. It provides an interface between planet formation simulations and exoplanet surveys such as Kepler. EPOS can also be used to estimate planet occurrence rates and the orbital architectures of planetary systems.

EPOS is written in Python and hosted on `github  <https://github.com/GijsMulders/epos>`_. 
Follow the instructions there to download and install.

**NEW** You can now install EPOS with pip:
::
   pip install epospy

(possibly preceded by ``sudo -H``)

How to use
----------

EPOS runs from the command line. The quickest way to familiarize yourself with the code is to copy the test and example scripts to a local directory:

>>> import EPOS
>>> EPOS.scripts.install()

These scripts demonstrate some of the basic functionality of EPOS.
Output wil appear in the terminal and plots in the ``png/`` folder.

.. toctree::
   :maxdepth: 2
   :caption: Detailed Instructions:

   instructions

API Documentation
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api

License / Attribution
---------------------
Copyright 2018 Gijs Mulders

Please cite `Mulders et al. 2018 <http://adsabs.harvard.edu/abs/2018arXiv180508211M>`_ if you use epos. You can also cite the github repository directly

.. image:: https://zenodo.org/badge/132065382.svg
   :target: https://zenodo.org/badge/latestdoi/132065382
   
List of Publications
--------------------
Mulders et al. 2018 
	`The Exoplanet Population Observation Simulator. I - The Inner Edges of Planetary Systems <http://adsabs.harvard.edu/abs/2018arXiv180508211M>`_
Pascucci et al. 2018 
	`A Universal Break in the Planet-to-star Mass-ratio Function of Kepler MKG Stars <http://adsabs.harvard.edu/abs/2018ApJ...856L..28P>`_
Kopparapu et al. 2018 
	`Exoplanet Classification and Yield Estimates for Direct Imaging Missions <http://adsabs.harvard.edu/abs/2018ApJ...856..122K>`_
Fernandes et al. 2019 
	`Hints for a Turnover at the Snow Line in the Giant Planet Occurrence Rate <http://adsabs.harvard.edu/abs/2018arXiv181205569F>`_

Version Notes:
--------------

:1.0.1: first public release
:1.0.2: pip installable version
:1.1.0: radial velocity without Monte Carlo |ImageLink|_

.. |ImageLink| image:: https://zenodo.org/badge/132065382.svg
.. _ImageLink: https://zenodo.org/badge/latestdoi/132065382

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
