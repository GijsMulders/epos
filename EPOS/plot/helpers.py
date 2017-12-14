''' helpers '''
import os
from matplotlib.colors import Normalize
from matplotlib import cm

def set_axes(ax, epos, **args):
	set_axis_distance(ax, epos, **args)
	set_axis_size(ax, epos, **args)

def set_axis_distance(ax, epos, Trim=False, Eff=False, In=False, IsX=True):

	if IsX:
		ax.set_xlabel('Orbital Period [days]')
		if Trim:
			ax.set_xlim(epos.xtrim)
		elif Eff:
			ax.set_xlim(epos.eff_xlim)
		else:
			ax.set_xlim(epos.obs_xlim)
		
		ax.set_xscale('log')
	else:
		ax.set_ylabel('Orbital Period [days]')
		if Trim:
			ax.set_ylim(epos.xtrim)
		elif Eff:
			ax.set_ylim(epos.eff_xlim)
		else:
			ax.set_ylim(epos.obs_xlim)
		
		ax.set_yscale('log')

	#xticks=[3, 10, 30, 100, 300]
	#ax.set_xticks(xticks)
	#ax.set_xticklabels(xticks)
	
	pass

def set_axis_size(ax, epos, Trim=False, Eff=False, In=False, IsY=True):

	if In:			label= r'M [M$_\bigoplus$]'
	elif epos.RV:	label= r'M sin i [M$_\bigoplus$]'
	else:			label= r'Planet Radius [R$_\bigoplus$]'

	if IsY:
		ax.set_ylabel(label)

		if Trim:
			ax.set_ylim(epos.in_ytrim if In else epos.ytrim)
		elif Eff:
			ax.set_ylim(epos.eff_ylim)
		else:
			ax.set_ylim(epos.obs_ylim)
		
		ax.set_yscale('log')
	else:
		ax.set_xlabel(label)

		if Trim:
			ax.set_xlim(epos.in_ytrim if In else epos.ytrim)
		elif Eff:
			ax.set_xlim(epos.eff_ylim)
		else:
			ax.set_xlim(epos.obs_ylim)
		
		ax.set_xscale('log')

	ax.set_yticks(epos.yticks)
	ax.set_yticklabels(epos.yticks)
	
def set_pyplot_defaults():
	from matplotlib import rcParams
	#print rcParams
	rcParams.update({'font.size': 16})
	rcParams.update({'legend.fontsize': 'medium'})
	rcParams.update({'axes.linewidth': 2.0})
	rcParams.update({'lines.linewidth': 2.0})
	rcParams.update({'patch.linewidth': 2.0})
	rcParams.update({'xtick.major.size': 8.0})
	rcParams.update({'xtick.minor.size': 4.0})
	rcParams.update({'ytick.major.size': 8.0})
	rcParams.update({'ytick.minor.size': 4.0})
	rcParams.update({'xtick.major.width': 1.0})
	rcParams.update({'xtick.minor.width': 1.0})
	rcParams.update({'ytick.major.width': 1.0})
	rcParams.update({'ytick.minor.width': 1.0})
	rcParams.update({'lines.markeredgewidth': 2.0})

def save(plt, name, dpi=150):

	# make sure path exists
	ipath= name.rfind('/')
	if ipath!= -1:
		#print name[:ipath]
		if not os.path.isdir(name[:ipath]): os.makedirs(name[:ipath])
	
	plt.savefig(name+'.png',bbox_inches='tight', dpi=dpi)
	plt.close()

def color_array(vals, vmin=None, vmax=None, cmap='jet'):
	# creates colors usable by matplotlib colorbar
	norm = Normalize(vmin=vmin, vmax=vmax)
	# plt.cm.get_cmap
	#Can put any colormap you like here.
	colours = cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(vals)
	return colours, norm