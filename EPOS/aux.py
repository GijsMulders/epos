import os
def save(plt, name, dpi=150):

	# make sure path exists
	ipath= name.rfind('/')
	if ipath!= -1:
		#print name[:ipath]
		if not os.path.isdir(name[:ipath]): os.makedirs(name[:ipath])
	
	plt.savefig(name+'.png',bbox_inches='tight', dpi=dpi)
	plt.close()

def display_or_save(plt,name='',save='png', dpi=150):
	if save is not '': 
		plt.savefig(name,bbox_inches='tight', dpi=dpi)
		plt.close()
	else:
		plt.show()

def save_path(save, name, verbose=True, **kwargs):
	import os
	''' makes dir for saving files (if not present), returns path '''
	if save is not '':
		path='{}/{}/'.format(save,name)
#		if not os.path.exists(save): os.makedirs(save)
		if not os.path.exists(save+'/'+name): os.makedirs(save+'/'+name) # recursive
		if verbose: print 'Plots in '+path

	kwargs['name']=name
	kwargs['save']=save
	kwargs['path']=path

	return kwargs

def plot_defaults(Reset=False, fac=2., fontsize=16, legendfont='medium'):
	''' This function increases line width and text sizes for default plots.
	Reset=true reverts to default sizes / widths.  '''
	
	if Reset:
		from matplotlib import rcdefaults
		rcdefaults() 
	
	else:
		from matplotlib import rcParams
		
		#print rcParams
		
		rcParams.update({'font.size': fontsize})
		#print rcParams['legend.fontsize']
		#print rcParams['figure.dpi']
		rcParams.update({'legend.fontsize': legendfont})
		#rcParams.update({'figure.dpi': 150}) # currently overridden by display_or_save
	
		# patch.linewidth is for histograms
		for flag in ['axes.linewidth', \
			'lines.linewidth',\
			'patch.linewidth',\
			'xtick.major.size', 'xtick.minor.size', 'ytick.major.size', 'ytick.minor.size',\
			'xtick.major.width', 'xtick.minor.width', 'ytick.major.width', 'ytick.minor.width']:
			rcParams[flag] *= fac
		
		rcParams['lines.markeredgewidth']*= 2.*fac # error bar caps but not lines
