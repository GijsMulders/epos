''' 
Save output to files
json/.../occurrence.??.json
'''
import numpy as np
import os, json

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

		if 'planet' in epos.occurrence:
			save_to_json(epos,'occurrence.planet',focc,['x','y','completeness','obs'])

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
	
def population(epos, Verbose=False):
	if hasattr(epos,'population'):
		if not os.path.isdir(epos.jsondir): os.makedirs(epos.jsondir)
		
		pop= epos.population
		
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
	
def save_to_json(epos, fname, bigdict, keys):
	tosave= {x: bigdict[x] for x in keys if x in bigdict}
	with open('{}{}.json'.format(epos.jsondir,fname), 'w') as f:
		json.dump(tosave, f, default=serialize_numpy_array)
	
def serialize_numpy_array(obj):
	if isinstance(obj, np.ndarray):
		return obj.tolist()
		raise TypeError('Not serializable')
