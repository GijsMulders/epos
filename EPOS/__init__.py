__shortversion__= u'3.1'
__version__= u'3.1.0.dev0' # [N!]N(.N)*[{a|b|rc}N][.postN][.devN]

__all__ = ['epos','fitparameters','kepler','rv','run','population','plot','occurrence',
	'fitfunctions','massradius','regression','multi','analytics','save',
	'scripts', 'cgs','classes']

from . import kepler, rv, run, plot, occurrence, population
from . import fitfunctions, regression, massradius, multi, analytics, save
from . import scripts, cgs

# is this correct python 3?
from . import classes
from EPOS.classes import epos, fitparameters
