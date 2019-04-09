''' 
Save output to files
json/.../occurrence.??.json
'''
import numpy as np
import os, json

def survey(epos, Verbose=False):
	if hasattr(epos,'obs_zoom'):
		if not os.path.isdir(epos.jsondir): os.makedirs(epos.jsondir)

		print '\nStoring Observations'

		if Verbose:
			print '  keys in epos.obs_zoom:'
			print '  ',epos.obs_zoom.keys()

		keys=epos.obs_zoom.keys()
		save_to_json(epos,'obs_zoom',epos.obs_zoom, keys)

		if 'multi' in epos.obs_zoom:
			keys= epos.obs_zoom['multi'].keys()
			if Verbose: print '  keys in multi:'
			print '  ',keys
			save_to_json(epos,'obs_zoom.multi',epos.obs_zoom['multi'], keys)
	else:
		print 'No survey stored, did you run epos.set_survey() and set_observation() ?'

def occurrence(epos, Verbose=False):
	if hasattr(epos,'occurrence'):
		if not os.path.isdir(epos.jsondir): os.makedirs(epos.jsondir)
		
		focc= epos.occurrence
		
		if Verbose:
			print 'keys in epos.occurrence'
			for key in focc:
				print '\n{}'.format(key)
				for subkey in focc[key]:
					print '  {}'.format(subkey)

		if 'model' in epos.occurrence:
			save_to_json(epos,'occurrence.model',focc['model'], ['eta', 'completeness'])

		if 'planet' in epos.occurrence:
			save_to_json(epos,'occurrence.planet',focc['planet'],['xvar','yvar','completeness','occ'])

		if 'bin' in epos.occurrence:
			gridkeys= ['xc', 'yc', 'dlnx', 'dlny', 'y', 'x']
			
			# inverse detection efficiency
			if 'occ' in focc['bin']:
				keys=['occ', 'err', 'n']+gridkeys
				save_to_json(epos, 'occurrence.inverse', focc['bin'], keys)
			
			# initial guess
			if 'eta0' in focc['bin']:
				keys=['eta0','gamma0','area']+gridkeys
				save_to_json(epos, 'occurrence.initial', focc['bin'], keys)
			
			# posterior
			if 'eta' in focc['bin']:
				keys=['eta','eta+','eta-','gamma','gamma+','gamma-','area']+gridkeys
				save_to_json(epos, 'occurrence.posterior', focc['bin'], keys)

		# occurrence along x and y axes
		if 'xzoom' in epos.occurrence and 'yzoom' in epos.occurrence:
			pass
		
	else:
		print 'No occurrence rates found, did run use epos.occurrence.all() ?'

def synthetic_survey(epos, Verbose=False):
	if hasattr(epos,'synthetic_survey'):
		if not os.path.isdir(epos.jsondir): os.makedirs(epos.jsondir)

		ss= epos.synthetic_survey

		print '\nStoring Synthetic Survey'

		if Verbose:
			print '  keys in epos.synthetic_survey:'
			print '  ',ss.keys()

		keys=['P', 'Y', 'M', 'R', 'P zoom', 'Y zoom']
		save_to_json(epos,'synthetic_survey',ss, keys)

		if 'multi' in ss:
			keys= ss['multi'].keys()
			if Verbose: print '  keys in multi:'
			print '  ',keys
			save_to_json(epos,'synthetic_survey.multi',ss['multi'], keys)
	else:
		print 'No synthetic survey stored, did you run epos.run.once() ?'
	
def population(epos, Verbose=False):
	if hasattr(epos,'population'):
		if not os.path.isdir(epos.jsondir): os.makedirs(epos.jsondir)
		
		pop= epos.population

		print '\nStoring Planet Population'

		if Verbose:
			print 'keys in epos.population'
			for key in pop:
				print '\n{}'.format(key)
				for subkey in pop[key]:
					print '  {}'.format(subkey)

		if 'system' in epos.population:
			keys=['P', 'Y', 'ID', 'inc', 'detectable']
			save_to_json(epos,'population.systems',pop['system'], keys)
		
	else:
		print 'No population stored, did you run epos.run.once() ?'

def model(epos, Verbose=False):
	if hasattr(epos,'pfm'):
		if not os.path.isdir(epos.jsondir): os.makedirs(epos.jsondir)

		print '\nStoring Models'

		if Verbose:
			print '  keys in epos.pfm:'
			print '  ',epos.pfm.keys()

		keys=epos.pfm.keys()
		save_to_json(epos,'pfm',epos.pfm, keys)

	else:
		print 'No model stored, did you run epos.set_population() ?'
	
def save_to_json(epos, fname, bigdict, keys):
	tosave= {x: bigdict[x] for x in keys if x in bigdict}
	with open('{}{}.json'.format(epos.jsondir,fname), 'w') as f:
		json.dump(tosave, f, default=serialize_numpy_array)
	
def serialize_numpy_array(obj):
	if isinstance(obj, np.ndarray):
		return obj.tolist()
	else:
		raise TypeError('Not serializable')
