import numpy as np

def save(dict,fname,tag='KID'):
	import numpy as np
	
	Mkeys=dict.keys()
	Munits=[]
	for key in Mkeys:
		if type(dict[key][0])==type(np.array([0.])[0]): Munits.append('float')
		elif type(dict[key][0])==type(np.array([0])[0]): Munits.append('int')
		elif type(dict[key][0])==type(np.array([''])[0]): Munits.append('string')
		elif type(dict[key][0])==type(np.array([True])[0]): Munits.append('bool')
		else: print 'no matching type for {} ({})'.format(key,type(dict[key][0]))
					
	with open(fname, 'w') as f:
		f.write(",".join(Mkeys)+'\n')
		f.write(",".join(Munits)+'\n')
		#for i in range(0,len(dict[tag])):
		for i in range(len(dict[Mkeys[0]])): 
			f.write(",".join([str(dict[ikey][i]) for ikey in Mkeys])+'\n')
			
def load(fname):
	import numpy as np
	import csv
	
	with open(fname, 'rb') as f:
		reader = csv.reader(f)
		Mkeys=reader.next()
		Munits=reader.next()
		dict={}
		for key in Mkeys: dict[key]=[]
		for row in reader:
			for j in range(0,len(Mkeys)):
				dict[Mkeys[j]].append(row[j])
		npdict={}
		for j,key in enumerate(Mkeys):
			if Munits[j]=='float': npdict[key]=np.array(dict[key],dtype=float)
			elif Munits[j]=='int': npdict[key]=np.array(dict[key],dtype=int)
			elif Munits[j]=='bool': npdict[key]=(np.array(dict[key])=='True')
			elif Munits[j]=='string': npdict[key]=np.array(dict[key])
			else: print 'no matching type for {} ({})'.format(key,type(dict[key][0]))
			
	return npdict
	
# 	updates fields in dictionary 1 with values from dictionary 2
#	The field being matched is the first entry in field1
# 	Fields updated can have different names, in this case supply field2 
# 	needs a check if field exists
def update(dict1,dict2,field1=[],field2=[],verbose=False,report_missing=False):
	import numpy as np
	
	if len(field1) < 2: print 'update.dictionary: need at least two fields!!'
	if len(field2) == 0: field2=field1
	if len(field1) != len(field2): print 'update.dictionary: fields need to be same length!!'


#	np.set_printoptions(threshold='nan')
		
	print "  matching field '{}' to '{}'".format(field1[0],field2[0])
	array1=dict1[field1[0]]
	array2=dict2[field2[0]]
	match1=[]  # get matching indices
	match2=[]
	for j in range(0,len(array2)):
		index=np.where(array2[j] == array1)
		if len(index[0]) == 0: 
			if verbose: print "field '{}' with value {} has no matches!".format(field2[0],array2[j])
			if report_missing: print "field '{}' with value {} has no matches!".format(field2[0],array2[j])
		elif len(index[0]) == 1: 
			match1.append(int(index[0]))
			match2.append(j)
			if verbose: print '	Couple {} to {}   (Matched {} == {})'.format(match1[-1],match2[-1],array1[match1[-1]],array2[match2[-1]])
		else: 
			match1.extend([int(i) for i in index[0]])
			match2.extend([j]*len(index[0]))
			if verbose:
				print 'update.dictionary: more than one match (???)'
				print '	Couple {} to {}'.format(match1[-1],match2[-1])
				for ii in range(0,len(index[0])):
					print '   (Matched {} == {})'.format(array1[match1[-1-ii]],array2[match2[-1-ii]])

	if len(match1) != len(match2): print 'update.dictionary: something went wrong!'
	print '  {} matching items out of {} and {}'.format(len(match1),len(array1),len(array2))
	
	for i in range(1,len(field1)):
		if field1[i] in dict1.keys():
			if verbose: print
			print "   copying field '{}' to '{}'".format(field2[i],field1[i])
			if verbose:
				for j in range(0,len(match1)): print '	Replace {} by {}'.format(dict1[field1[i]][match1[j]],dict2[field2[i]][match2[j]])
			dict1[field1[i]][match1]=dict2[field2[i]][match2]
		else:
			if (len(match1) == len(match2)):
				print "   copying field '{}' to new field '{}'".format(field2[i],field1[i])
				dict1[field1[i]]=np.zeros(len(match1),dtype=dict2[field2[i]].dtype)
				dict1[field1[i]][match1]=dict2[field2[i]][match2]
			else:
				print "  !! no field '{}' to merge into !!".format(field1[i])
		
	print
	
def remove_nonmatching_elements(dict1,dict2,field='KID',report_missing=False, ReturnNonMatch=False):
	# remove elements from (dict1) that are not found in (field) from (dict2)
	# ReturnNonMatch actually keeps those elements
	import numpy as np
	
	match=[]
	if ReturnNonMatch: nonmatch= []
	for i,KID in enumerate(dict1[field]):
		index=np.where(KID == dict2[field])
		if len(index[0]) == 0:
			if ReturnNonMatch: nonmatch.append(i)
			if report_missing: print "field '{}' with value {} has no matches!".format(field,KID)
		else:
			match.append(i)
	
	if ReturnNonMatch:
		print '  {} elements, {} matching, {} nonmatching'.format(len(dict1[field]), len(match), len(nonmatch))
	
	dict={}
	tomatch=nonmatch if ReturnNonMatch else match
	for key in dict1.keys():
		dict[key]=dict1[key][tomatch]
	
	print '  removed {} elements, {} left'.format(len(dict1[field])-len(tomatch), len(tomatch))
	
	return dict
	
def remove_elements(dict1,TFarray,verbose=True):
	# remove elements from (dict1) that are True in (TFarray)
	import numpy as np
	
	remove=np.where(TFarray == True)
	keep=np.where(TFarray == False)
	
	dict={}
	for key in dict1.keys():
		if type(dict1[key])==np.ndarray and len(dict1[key]) == len(TFarray):
			dict[key]=dict1[key][keep]
		else:
			print "could not remove elements from '{}'".format(key)
			dict[key]=dict1[key]
	
	if verbose: print '  removed {} elements, {} left'.format(len(remove[0]), len(keep[0]))
	
	return dict
	
def subset(dict1, bin, lenkey='KID', verbose=True, report_missing=True):
	# keep only elements from (dict1) based on np.where bin
	# returns a COPY
	import numpy as np
	
	length= len(dict1[lenkey])
	dict={}
	for key in dict1.keys():
		if (type(dict1[key])==np.ndarray) and len(dict1[key]) == length:
			dict[key]=dict1[key][bin]
		else:
			if report_missing: print "could not remove elements from '{}'".format(key)
			dict[key]=dict1[key]
	
	if verbose: print '  extracted {} elements out of {}'.format(len(bin[0]), length)
	
	return dict

def concatenate(dict1, dict2, verbose=False):
	ret={}
	for key in dict1:
		#if key in dict2 and type(dict1[key]) is np.ndarray:
		if key in dict2:
			ret[key]= np.concatenate( [dict1[key], dict2[key]] )
			# print key
		#else: print 'not found in dict2: ' + key
	
	if verbose: print 'Merged keys: {}'.format(ret.keys())
	return ret
		