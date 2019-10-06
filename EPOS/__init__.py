__shortversion__= u'3.0'
__version__= u'3.0.1.dev1' # [N!]N(.N)*[{a|b|rc}N][.postN][.devN]


__all__ = ['epos','fitparameters','kepler','rv','run','population','plot','occurrence',
	'fitfunctions','pfmodel','massradius','regression','multi','analytics','save',
	'scripts', 'cgs','classes']

from . import kepler, rv, run, plot, occurrence, population
from . import fitfunctions, pfmodel, regression, massradius, multi, analytics, save
from . import scripts, cgs

# is this correct python 3?
from . import classes
from EPOS.classes import epos, fitparameters
