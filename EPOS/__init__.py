__shortversion__= u'1.1'
__version__= u'1.1.1.dev1' # [N!]N(.N)*[{a|b|rc}N][.postN][.devN]

__all__ = ['epos','fitparameters','kepler','rv','run','population','plot','occurrence',
	'fitfunctions','pfmodel','massradius','regression','multi','analytics','save',
	'scripts', 'cgs']

import kepler, rv, run, plot, occurrence, population
import fitfunctions, pfmodel, regression, massradius, multi, analytics, save
import scripts, cgs

from classes import epos, fitparameters