''' helpers '''
import os
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib import gridspec

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
		if hasattr(epos,'xticks'):
			ax.set_xticks(epos.xticks)
			ax.set_xticklabels(epos.xticklabels)
	
	else:
		ax.set_ylabel('Orbital Period [days]')
		if Trim:
			ax.set_ylim(epos.xtrim)
		elif Eff:
			ax.set_ylim(epos.eff_xlim)
		else:
			ax.set_ylim(epos.obs_xlim)
		
		ax.set_yscale('log')
		if hasattr(epos,'xticks'):
			ax.set_yticks(epos.xticks)
			ax.set_yticklabels(epos.xticklabels)

def set_axis_size(ax, epos, Trim=False, Eff=False, In=False, IsY=True):

	ticks=None
	if In:
		label= r'M [M$_\bigoplus$]'
		if hasattr(epos,'y_inticks'):
			ticks= epos.y_inticks
			ticklabels= epos.y_inticklabels
	elif epos.RV:
		label= r'M sin i [M$_\bigoplus$]'
		if hasattr(epos,'yticks'):
			ticks= epos.yticks
			ticklabels= epos.yticklabels
	else:
		label= r'Planet Radius [R$_\bigoplus$]'
		if hasattr(epos,'yticks'):
			ticks= epos.yticks
			ticklabels= epos.yticklabels

	if IsY:
		ax.set_ylabel(label)

		if Trim:
			ax.set_ylim(epos.in_ytrim if In else epos.ytrim)
		elif Eff:
			ax.set_ylim(epos.eff_ylim)
		else:
			ax.set_ylim(epos.obs_ylim)
		
		ax.set_yscale('log')
		if hasattr(epos,'yticks'):
		#if ticks is not None:
			ax.set_yticks(ticks)
			ax.set_yticklabels(ticklabels)
		
	else:
		ax.set_xlabel(label)

		if Trim:
			ax.set_xlim(epos.in_ytrim if In else epos.ytrim)
		elif Eff:
			ax.set_xlim(epos.eff_ylim)
		else:
			ax.set_xlim(epos.obs_ylim)
		
		ax.set_xscale('log')
		if hasattr(epos,'yticks'):
			ax.set_xticks(ticks)
			ax.set_xticklabels(ticklabels)

def make_panels(plt, Fancy=False):
	if Fancy:
		gs = gridspec.GridSpec(2, 2,
	                       width_ratios=[20,4],
	                       height_ratios=[3, 10]
	                       )
		f= plt.figure()
		f.subplots_adjust(wspace=0, hspace=0)

		axR = plt.subplot(gs[1, 1])
		axP = plt.subplot(gs[0, 0])

		ax = plt.subplot(gs[1, 0], sharex=axP, sharey=axR)	
		#axb = plt.subplot(gs[0, 2])	



		#axh.tick_params(direction='in', which='both', left=False, right=True, labelleft=False)
		#axh.yaxis.set_label_position('right')
		axR.axis('off')
		axP.axis('off')
		
		ax.tick_params(direction='out', which='both', top=False, right=False, 
			bottom=True, left=True)
		axP.tick_params(direction='out', which='both', top=False, right=False, 
			bottom=True, left=False)
		axR.tick_params(direction='out', which='both', top=False, right=False, 
			bottom=False, left=True)
	else:		
		gs = gridspec.GridSpec(2, 2,
	                       width_ratios=[4, 20],
	                       height_ratios=[10, 3]
	                       )
		f= plt.figure()
		f.subplots_adjust(wspace=0, hspace=0)
		
		ax = plt.subplot(gs[0, 1])	
		#axb = plt.subplot(gs[0, 2])	

		axR = plt.subplot(gs[0, 0], sharey=ax)
		axP = plt.subplot(gs[1, 1], sharex=ax)
		
		ax.tick_params(direction='in', which='both', top=True, right=True, 
			bottom=False, left=False)

	return f, (ax, axR, axP)

def make_panels_right(plt):
	gs = gridspec.GridSpec(2, 2,
                       width_ratios=[20, 4],
                       height_ratios=[10, 3]
                       )
	f= plt.figure()
	f.subplots_adjust(wspace=0, hspace=0)
	
	ax = plt.subplot(gs[0, 0])	
	#axb = plt.subplot(gs[0, 2])	

	axR = plt.subplot(gs[0, 1], sharey=ax)
	axP = plt.subplot(gs[1, 0], sharex=ax)

	ax.tick_params(direction='in', which='both', top=True)
	ax.tick_params(direction='out', which='both', left=True, top=False)
	
	axR.yaxis.tick_right()
	#axR.tick_params(direction='out', which='both', right=True, left=False)
	#axR.yaxis.set_label_position('right')
			
	return f, (ax, axR, axP)

def make_panels_clrbar(plt):
	gs = gridspec.GridSpec(2, 3,
                       width_ratios=[4, 20, 1],
                       height_ratios=[10, 3]
                       )
	f= plt.figure()
	f.subplots_adjust(wspace=0, hspace=0)
	
	ax = plt.subplot(gs[0, 1])	
	axb = plt.subplot(gs[0, 2])	

	axR = plt.subplot(gs[0, 0], sharey=ax)
	axP = plt.subplot(gs[1, 1], sharex=ax)

	return f, (ax, axb, axR, axP)

def set_pyplot_defaults():
	from matplotlib import rcParams
	#print rcParams
	rcParams.update({'font.size': 14}) # 16
	rcParams.update({'legend.fontsize': 12}) # medium
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
	rcParams.update({'axes.formatter.min_exponent': 3})
	#rcParams.update({'hatch.linewidth':  1.0}) # error in pyplot <2.0?

def default_pyplot2_colors(colors):
	colors.ColorConverter.colors['C0']='#1f77b4'
	colors.ColorConverter.colors['C1']='#ff7f0e'
	colors.ColorConverter.colors['C2']='#2ca02c'
	colors.ColorConverter.colors['C3']='#d62728'
	colors.ColorConverter.colors['C4']='#9467bd'
	colors.ColorConverter.colors['C5']='#8c564b'
	colors.ColorConverter.colors['C6']='#e377c2'
	colors.ColorConverter.colors['C7']='#7f7f7f'
	colors.ColorConverter.colors['C8']='#bcbd22'
	colors.ColorConverter.colors['C9']='#17becf'

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