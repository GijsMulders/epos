import os.path, shutil
import glob

# is there a less clumsy way to do this?
import EPOS
fpath= os.path.dirname(EPOS.__file__)

def install(fdir='epos-scripts'):
	'''
	Installs the tests and examples to the directory epos-scripts/
	'''
	#if not os.path.isdir(fdir): os.makedirs(fdir)

	for subdir in ['tests','examples','papers']:
		scriptdir= '{}/{}'.format(fdir, subdir)
		if not os.path.isdir(scriptdir): os.makedirs(scriptdir)

		print '\nCopying scripts into {}'.format(scriptdir)

		flist= glob.glob('{}/scriptdir/{}/*.py'.format(fpath, subdir))
		for fname in flist:
			fname_stripped= fname.split('/')[-1]
			fdest= '{}/{}'.format(scriptdir, fname_stripped)
			shutil.copyfile(fname, fdest)
			# make them executable
			print '  copied {}'.format(fname_stripped)

	print '\nRun scripts with\n'
	print '  ./test_1_survey.py\n'
	print 'or\n'
	print '  ipython test_1_survey.py\n'