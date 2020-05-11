import fnmatch
import datetime
import subprocess
from pylab import *

import os

def get_colors():
    return ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','8c564b']

def bfigure():
    F = cfigure()
    ax=  F.add_axes([0.25,0.25,0.7,0.7])
    return F,ax

def cfigure(w=15,h=10,police=33,d=0):
	#fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
	#inches_per_pt = 1.0/72.27               # Convert pt to inch
	#golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	#fig_width = fig_width_pt*inches_per_pt  # width in inches
	#fig_height = fig_width*golden_mean      # height in inches
	wi = 0.39*w
	hi = 0.39*h
	fig_size =  [wi,hi]
	params = {#'backend': 'PS',
		  #'ps.usedistiller' : 'xpdf',
		  'axes.labelsize': police,
		  'axes.formatter.limits' : [-4,4], 
		  'lines.linewidth': 2,
		  'lines.markersize': 11,
		  'font.size': police,
		  'legend.fontsize': 27,
		  'xtick.labelsize': police-7,
		  'ytick.labelsize': police-7,
		  'ytick.major.pad' : 10,
		  'xtick.major.pad' : 10,
		  'text.usetex': True,
	#	  'font' : {'family': 'serif', 'serif': ['Computer Modern']},
		  #3'font.family':'sans-serif',
		  #'font.sans-serif':'Helvetica',
		  #'mathtext.fontset' : 'stix',
		  'font.family' : 'STIXGeneral',
		  #'font.serif' : 'C',
		  #'font.family' : 'serif',
		  'font.weight' : 'light',
		  'figure.figsize': fig_size
}
	rcParams.update(params)
	F = figure(facecolor='w',edgecolor='k')
	if d==1:
		suptitle(datetime.date.today(),fontsize=police-10,horizontalalignment='right',x=0.92)	
	return F



def save_figure(fname):
	savefig(fname)
	if fnmatch.fnmatch(fname, '*.eps'):
			subprocess.Popen([r"ps2eps",fname]).wait()
			subprocess.Popen([r"mv",fname+".eps",fname])
def crop_dir(dname):
	I = os.listdir(dname)
	for f in I:
		fname = dname+f
		print(fname)
		if fnmatch.fnmatch(fname, '*.eps'):
				subprocess.Popen([r"ps2eps",fname]).wait()
				subprocess.Popen([r"mv",fname+".eps",fname])


#def subtext(cx,cy,t):
#    xmin,xmax = xlim()
#    ymin,ymax = ylim()
#    police = rcParams['font.size']	
#    if gca().get_xscale()=='linear':
#    x = xmin+cx*(xmax-xmin)
#    else:
#    x = xmin*(xmax/xmin)**cx
#    if gca().get_yscale()=='linear':
#    y = ymin+cy*(ymax-ymin)
#    else:
#    y = ymin*(ymax/ymin)**cy
#    return text(x,y,t,fontsize=police+2)

def create_figure_eps(w,h,police,d):
	#fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
	#inches_per_pt = 1.0/72.27               # Convert pt to inch
	#golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	#fig_width = fig_width_pt*inches_per_pt  # width in inches
	#fig_height = fig_width*golden_mean      # height in inches
	wi = 0.39*w
	hi = 0.39*h
	fig_size =  [wi,hi]
	params = {'backend': 'PS',
		  'ps.usedistiller' : 'xpdf',
		  'axes.labelsize': police,
		  'axes.formatter.limits' : [-4,4], 
		  'lines.linewidth': 2,
		  'lines.markersize': 11,
		  'font.size': police,
		  'legend.fontsize': police-7,
		  'xtick.labelsize': police-2,
		  'ytick.labelsize': police-2,
		  'ytick.major.pad' : 10,
		  'xtick.major.pad' : 10,
		  'text.usetex': True,
		  #'font' : {'family': 'serif', 'serif': ['Computer Modern']},
		  'font.family' : 'serif',
		  'font.serif' : 'Computer Modern',
		  'figure.figsize': fig_size}
	rcParams.update(params)
	F = figure(facecolor='w',edgecolor='k')
	if d==1:
		suptitle(datetime.date.today(),fontsize=police-10,horizontalalignment='right',x=0.92)	
	return F
