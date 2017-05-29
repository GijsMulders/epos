def all():
	from __init__ import __all__ # better way to do this?
	import clearscreen
	print 'Description of EPOS module by Gijs Mulders'
	print '(Exoplanet Population Observation Simulator)'
	print
	print 'EPOS contains sub-modules:'
	print __all__
	print 
	print 'run EPOS.submodule.readme() for a description of each module'
	#print 'TODO: add automated way of listing all readmes?'
	print
	#print 'TODO: run EPOS.test.verify() to check code integrity'
	print 'run example_parametric.py to see a demonstration'
	print
	