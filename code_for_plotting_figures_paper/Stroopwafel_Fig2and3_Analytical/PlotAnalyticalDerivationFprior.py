#  Last updated 31-12-2018
#  Plots Figure 2 and 3 in STROOPWAFEL paper: 
#  plots of the analytical gain and fexpl with colored bars inside plot to give example of rates / rareness. 
#  




from __future__ import division 
# import seaborn as sns
#Import functions
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import matplotlib
from scipy.stats import multivariate_normal
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D 

from scipy.integrate import quad
from scipy import stats
from itertools import repeat
import matplotlib.patches as mpatches
import matplotlib.patches as patches
# from matplotlib.patches import Ellipse
import random
from scipy import stats

from matplotlib import rc                                                                                                                                                                                                                    
from matplotlib import rcParams

import pylab as pl
# la = pl.matplotlib.font_manager.FontManager()
# lu = pl.matplotlib.font_manager.FontProperties(family = 'palatino')
# la.findfont(lu)


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


# matplotlib.rcParams['text.usetex'] = True 
# matplotlib.rcParams['text.latex.preamble'] = [
#        r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
#        r'\usepackage{helvet}',    # set the normal font here
#        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#        r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
# ]




rc('font', family='serif', weight = 'bold')
# # rc('weight', 'bold')

rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
# matplotlib.rcParams['text.latex.preamble'] = ['bold']

matplotlib.rcParams['text.latex.unicode'] = True

rc('axes', linewidth=2)

matplotlib.rcParams['xtick.major.size'] = 12
matplotlib.rcParams['ytick.major.size'] = 12
matplotlib.rcParams['xtick.minor.size'] = 8
matplotlib.rcParams['ytick.minor.size'] = 8
# colors  =   ['#189ef8' , '#ff7f00', '#4daf4a', '#f781bf', 'gray', 'gray'] ##a65628', '#984ea3', '#dede00', '#999999', '#e41a1c'] #'#377eb8'
colors  =   ['#189ef8' , '#ff7f00', '#4daf4a', '#f781bf', '#a65628',  '#984ea3'] # '#dede00' ]
# matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['font.weight']= 'bold'
# matplotlib.rcParams['text.weight']= 'bold'
matplotlib.rcParams.update({'font.weight': 'bold'})

# matplotlib.rcParams.update({'font.size': 22})

# the rates for the Z = 0.001 simulations are obtained from the NhitsVSNbinariesPlot_Z0_001 notebook. 
ratess = np.asarray([6711, 5163, 655, 893, 544,  332     ]) / 1E6 # from bdMC simulation   #double check 600
ratessZ0_002 = [0.007339, 0.005712,0.001593, 0.000759 , 0.00055, 0.00027] 
FracUnc = [0.00315, 0.00221, 0.0100, 0.00820, 0.00410, 0.0179 ]

fexpls = [0.234980765019, 0.270022919908, 0.660361018917, 0.593282440305, 0.689711861729, 0.772744090606 ]

# true efficiencies in the refinement phase:
Gains = [ 35, 53, 39, 45, 203, 24]

labelss = [r"\textbf{ALL }  ", r"\textbf{BH--BH} ", r"\textbf{BH--NS} ", r"\textbf{NS--NS} ", r" $\rm{\textbf{BH--BH}} \, \geq 50 \, \rm{M}_{\odot} $ ", r" $\rm{\textbf{NS--NS}} \, {\leq 50 \, \rm{Myr}} $ "]
# ax.tick_params(axis='both', which='major', length=12,width=2, pad=10,labelsize=30)
# ax.tick_params(axis='both', which='minor', length=12,width=2, pad=10)

fs = 34
fs_text = fs







def CalculateFractionPrior(rate, Ntot):
	d2 = 1./Ntot
	d1 = rate 

	numerator = d1 * (np.sqrt(1.-d1) - np.sqrt(d2))
	denominator = (np.sqrt(1-d1)*(np.sqrt(d2*(1.-d1))+d1) )
	fref = 	numerator / denominator
	return 1- fref


def CalculateFractionPriorIterated(rate, Ntot):
    d2 = 1./Ntot
    d1 = rate

    numerator = d1 * (np.sqrt(1.-d1) - np.sqrt(d2))
    denominator = (np.sqrt(1-d1)*(np.sqrt(d2*(1.-d1))+d1) )
    fref =  numerator / denominator

    fprior = fref
    # print fprior, 'fprior'
    for i in range(100):

        d2 = 1./(fprior*Ntot)
        d1 = rate

        numerator = d1 * (np.sqrt(1.-d1) - np.sqrt(d2))
        denominator = (np.sqrt(1-d1)*(np.sqrt(d2*(1.-d1))+d1) )
        fexpl =  numerator / denominator
        fprior = 1-fexpl


    return fprior





# From https://stackoverflow.com/questions/38208700/matplotlib-plot-lines-with-colors-through-colormap
def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale 
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc









############ BELOW IS WITH ITERATIVE FEXPL, FIG 1 #############


f, (axarr) = plt.subplots(1,1, figsize=(16,10))




rates = np.logspace(-5, -1, 200)


axarr.set_xlim(np.min(rates), np.max(rates))
axarr.set_ylim(0, 1)

## PLOT 
n_lines = 4
ListNtot = [10**4, 10**5, 10**6, 10**7]
LABELS = [r'$N= 10^4$', r'$N = 10^5$', r'$N = 10^6$', r'$N = 10^7$']
ALPHA = [1,1,1,1]#[0.2, 0.4, 0.6, 1]
LINESTYLES = [ ':','-.', '--','-',  ]


for  ind, Ntot in enumerate(ListNtot):
    # fexpl = CalculateFractionPrior(rates, Ntot)
    fexplIterated = CalculateFractionPriorIterated(rates, Ntot)
    axarr.plot(rates, fexplIterated, c='k', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100, linestyle = LINESTYLES[ind])
    # axarr.plot(rates, fexplIterated, c='b', lw = 6, alpha = ALPHA[ind])


axarr.grid(True, alpha = 0.5)



axarr.tick_params(labelsize=fs)

axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs, fontweight = 'bold') #$ \mathcal{R}_{\mathrm{x}}$


axarr.set_ylabel(r'\textbf{Fraction exploratory phase} $f_{\rm{expl}}$', fontsize=fs)

#cbar.ax
axarr.tick_params(labelsize=fs, pad = 10.1) 

axarr.text(0.000012, 0.2, r"\textbf{more rare}" +" \n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left', weight = 'bold' )
axarr.text(0.051, 0.2,    r"\textbf{more common}" + "\n "  +r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
axarr.legend(loc = 'upper right',  fontsize = fs, bbox_transform=axarr.transData) #, loc = 'upper center') #frameon=False, 

# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
plt.xscale('log')

for i in range(6):
    # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
    axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [1, 1], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
    axarr.text((ratess[i]), 0.015,  "  " + labelss[i] + " " , fontsize = fs-10, rotation=90, va = 'bottom', ha = 'center', zorder = 2)	
    axarr.scatter(ratess[i], fexpls[i],  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 500,  linewidth=4.5 )


plt.tight_layout()

plt.savefig('AnalyticalFpriorIterated.pdf')
plt.savefig('AnalyticalFpriorIterated.png', dpi = 300)




######### PLOT 2 ##########


f, (axarr) = plt.subplots(1,1, figsize=(16,10))

for  ind, Ntot in enumerate(ListNtot):
    fexpl = CalculateFractionPriorIterated(rates, Ntot)


    NhitsMC = (rates*Ntot)
    NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


    ratioNhits = NhitsAIS / NhitsMC

    axarr.plot(rates, ratioNhits, c='k', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100, linestyle = LINESTYLES[ind])
    





axarr.grid(True, alpha = 0.5)

axarr.tick_params(labelsize=fs)

axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs) #$ \mathcal{R}_{\mathrm{x}}$
axarr.set_ylabel(r'\textbf{Max. possible gain}', fontsize=fs)

axarr.tick_params(labelsize=fs, pad = 10.1) 

axarr.text(0.000012,20.3, r"\textbf{more rare}" + "\n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left' , weight ='bold')
axarr.text(0.051, 20.3,   r"\textbf{more common}" +" \n " + r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
axarr.legend( loc = 'upper right',  fontsize = fs, bbox_transform=axarr.transData) #, loc = 'upper center') #frameon=False, bbox_to_anchor=(0.1,3000)

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
plt.xscale('log')

plt.yscale('log')

for i in range(6):
    # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
    axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [10**4, 10**4], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
    axarr.text((ratess[i]), 0.9* 10**4 ,  labelss[i] + "   " , fontsize = fs-10, rotation=90, va = 'bottom', ha = 'center', zorder = 2)


    NhitsMCtemp = (ratess[i]*1E6)
    NhitsAIStemp = (ratess[i]*1E6*fexpls[i] + (1-fexpls[i])*1E6) 


    ratioNhitstemp = NhitsAIStemp / NhitsMCtemp

    axarr.scatter(ratess[i], ratioNhitstemp,  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 1000,  linewidth=4.5 ) 





plt.ylim(10**1, 10**4)
plt.xlim(10**-5, 10**-1)
plt.tight_layout()
plt.savefig('AnalyticalFpriorRatioIterated.pdf')
plt.savefig('AnalyticalFpriorRatioIterated.png', dpi = 300)








# 
######### PLOT 3c: UNCERTAINTY with N = 1E4 AND N = 1E6 ##########
#  this is the bottom plot of Figure 2 in STROOPWAFEL paper 

f, (axarr) = plt.subplots(1,1, figsize=(16,10))

for  ind, Ntot in enumerate([1E4]):   #(ListNtot):
    fexpl = CalculateFractionPriorIterated(rates, Ntot)
    ind = 0

    NhitsMC = (rates*Ntot)
    NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


    uncertaintyMC = 1 / np.sqrt(NhitsMC)
    uncertaintyAIS = (1 / np.sqrt(NhitsAIS)) + ((1 / (fexpl*Ntot))/rates) 

    uncertaintyRatio = uncertaintyAIS / uncertaintyMC

    axarr.plot(rates, uncertaintyAIS   , c='gray', lw = 6,  zorder = 100, linestyle = LINESTYLES[3])
    axarr.plot(rates, uncertaintyMC    , c='gray', lw = 6, zorder = 100, linestyle = LINESTYLES[0])
    axarr.fill_between(rates,  uncertaintyAIS, uncertaintyMC, color='gray',  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)

for  ind, Ntot in enumerate([1E6]):   #(ListNtot):
    fexpl = CalculateFractionPriorIterated(rates, Ntot)
    ind = 2

    NhitsMC = (rates*Ntot)
    NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


    uncertaintyMC = 1 / np.sqrt(NhitsMC)
    uncertaintyAIS = (1 / np.sqrt(NhitsAIS)) #+ ((1 / (fexpl*Ntot))/rates) 

    uncertaintyRatio = uncertaintyAIS / uncertaintyMC

    axarr.plot(rates, uncertaintyAIS   , c='gray', lw = 6, label = 'this study',  zorder = 100, linestyle = LINESTYLES[3])
    axarr.plot(rates, uncertaintyMC    , c='gray', lw = 6, label = 'traditional', zorder = 100, linestyle = LINESTYLES[0])
    axarr.fill_between(rates,  uncertaintyAIS, uncertaintyMC, color='gray',  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)


axarr.text(0.000018,0.065, r"\textbf{traditional uncertainty}", fontsize = fs-6,rotation=338, va = 'bottom', ha = 'left' )
axarr.text(0.00002,0.0083, r"\textbf{STROOPWAFEL}" + "\n" + r"\textbf{uncertainty}", fontsize = fs-6,rotation=340, va = 'bottom', ha = 'left' )


axarr.grid(True, alpha = 0.5)

# axarr.pines()

axarr.tick_params(labelsize=fs)

axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs) #$ \mathcal{R}_{\mathrm{x}}$


axarr.set_ylabel(r'\textbf{Fractional statistical uncertainty}', fontsize=fs)

#cbar.ax
axarr.tick_params(labelsize=fs, pad = 10.1) 

axarr.text(0.000012,0.0017203, r"\textbf{more rare}" +" \n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left' )
axarr.text(0.051, 0.0017203,   r"\textbf{more common}" +" \n " + r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
# axarr.legend( loc = 'upper left',  fontsize = fs, bbox_transform=axarr.transData, ncol = 1) #, loc = 'upper center') #frameon=False, bbox_to_anchor=(0.1,3000)

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
plt.xscale('log')

plt.yscale('log')

for i in range(6):

    # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
    axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [10**4, 10**4], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
    axarr.text((ratess[i]), 0.9*10**0,  labelss[i] + "   " , fontsize = fs-10, rotation=90, va = 'top', ha = 'center', zorder = 2)
    axarr.scatter(ratess[i], FracUnc[i],  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 500,  linewidth=4.5 )


axarr.text(0.009,0.033, r"$N = 10^4$", color = 'k',  fontsize = fs, va = 'bottom', ha = 'left' )    
axarr.text(0.0107,0.002, r"$N = 10^6$", color = 'k',  fontsize = fs, va = 'bottom', ha = 'left' )


plt.xlim(10**-5, 10**-1)
plt.ylim(10**-3, 10**0)
plt.tight_layout()
plt.savefig('AnalyticalFpriorIteratedUncertaintyCOMBINED.pdf')
plt.savefig('AnalyticalFpriorIteratedUncertaintyCOMBINED.png', dpi = 300)




######### PLOT 4 ##########
##  This is the Maximum gain VS Rate of target population plot (Fig 3 in STROOPWAFEL paper)

LABELS2 = [r'$1$', r'$1/3$', r'$0.1$', r'$1/20$']

f, (axarr) = plt.subplots(1,1, figsize=(16,10))

for  ind, Ntot in enumerate(ListNtot):
	if ind ==1:

	    fexpl = CalculateFractionPriorIterated(rates, Ntot)
	    # print 'fexpl = ', fexpl
	    # NhitsMC = (rates*Ntot)
	    # NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot)
	    # Nexpl = np.zeros_like(rates)
	    # maskNexplTooLow = (rates*Ntot*fexpl < 4)

	    # Nexpl[maskNexplTooLow] = 4 / rates[[maskNexplTooLow]] 

	    # makszero = (Ntot - Nexpl) <= 0
	    # NhitsAIS[makszero] = np.zeros_like(NhitsAIS)[makszero]

	    # maksLeft = ((Ntot - Nexpl) > 0) & (Nexpl >0)
	    # NhitsAIS[maksLeft] = (Ntot - Nexpl )[maksLeft]



	    NhitsMC = (rates*Ntot)
	    NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


	    ratioNhits = NhitsAIS / NhitsMC

	    axarr.plot(rates,      ratioNhits, c='r', lw = 6, label = LABELS2[0], alpha = ALPHA[ind], zorder = 100)
	    # axarr.plot(rates, 0.33* ratioNhits, c='r', lw = 6, label = LABELS2[1], alpha = ALPHA[ind], zorder = 100, linestyle = ':')
	    axarr.plot(rates, 0.1* ratioNhits, c='r', lw = 6, label = LABELS2[2], alpha = ALPHA[ind], zorder = 100, linestyle = '--')
	    # axarr.plot(rates, 0.05*ratioNhits, c='r', lw = 6, label = LABELS2[3], alpha = ALPHA[ind], zorder = 100, linestyle = '-.')
	    # axarr.plot(rates, NhitsAIS, c='r', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100)
	    # axarr.plot(rates, NhitsMC,  c='k', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100)




axarr.grid(True, alpha = 0.5)

# axarr.pines()

axarr.tick_params(labelsize=fs)

axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs) #$ \mathcal{R}_{\mathrm{x}}$


axarr.set_ylabel(r'\textbf{Gain} ', fontsize=fs)

#cbar.ax
axarr.tick_params(labelsize=fs, pad = 10.1) 

axarr.text(0.000012,20.3, r"\textbf{more rare}" +" \n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left' )
axarr.text(0.051, 20.3,   r"\textbf{more common}" +" \n " + r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
lg = axarr.legend( loc = 'upper left',  bbox_transform=axarr.transData,  title="title", fontsize = fs-5, ncol = 2 ) #, loc = 'upper center') #frameon=False, bbox_to_anchor=(0.1,3000)



axarr.text(0.07, 6.3*1E3, r"$N = 10^6$", fontsize = fs+15, rotation=0, va = 'top', ha = 'right' )

lg.set_title(r'\textbf{Efficiency of}' + '\n' + r'\textbf{refinement phase}')
lg.get_title().set_fontsize(fs-5)

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
plt.xscale('log')

plt.yscale('log')

for i in range(6):
    # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
    axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [10**4, 10**4], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
    axarr.text((ratess[i]), 0.9* 10**4,  labelss[i] + "   " , fontsize = fs-10, rotation=90, va = 'top', ha = 'center', zorder = 2)
    axarr.scatter(ratess[i], Gains[i],  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 1000,  linewidth=4.5 ) 

plt.ylim(10**1, 10**4)
plt.xlim(10**-5, 10**-1)
plt.tight_layout()
plt.savefig('AnalyticalFpriorRatioIteratedVaryEfficiency.pdf')
plt.savefig('AnalyticalFpriorRatioIteratedVaryEfficiency.png', dpi = 300)










def CalculateFractionPriorIterated_imperfect(rate, Ntot, v_factor):
    d2 = 1./Ntot
    d1 = rate


    numerator = v_factor*d1 - v_factor**2 *d1 - d1**2 +2*v_factor*(d1**2) + (v_factor**2)*d2 
    denominator = (-1+v_factor)*(-v_factor*d1 +v_factor**2*d1 + d1**2 +v_factor**2 * d2)

    


    # numerator = d1 * (np.sqrt(1.-d1) - np.sqrt(d2))
    # denominator = (np.sqrt(1-d1)*(np.sqrt(d2*(1.-d1))+d1) )
    fref =  1 - (numerator / denominator)

    fprior = (numerator / denominator)
    # print fprior, 'fprior'
    for i in range(100):

        d2 = 1./(fprior*Ntot)
        d1 = rate

        numerator = v_factor*d1 - v_factor**2 *d1 - d1**2 +2*v_factor*(d1**2) + (v_factor**2)*d2 
        denominator = (-1+v_factor)*(-v_factor*d1 +v_factor**2*d1 + d1**2 +v_factor**2 * d2)


        fprior =  numerator / denominator
        # fprior = 1-fexpl


    return fprior








#######################################################################
#######     IMPERFECT SAMPLING
#######################################################################

############ BELOW IS WITH ITERATIVE FEXPL, FIG 1 #############

# Vimperfect  =1

# f, (axarr) = plt.subplots(1,1, figsize=(16,10))




# rates = np.logspace(-5, -1, 200)


# axarr.set_xlim(np.min(rates), np.max(rates))
# axarr.set_ylim(0, 1)

# ## PLOT 
# n_lines = 4
# ListNtot = [10**4, 10**5, 10**6, 10**7]
# LABELS = [r'$N= 10^4$', r'$N = 10^5$', r'$N = 10^6$', r'$N = 10^7$']
# ALPHA = [1,1,1,1]#[0.2, 0.4, 0.6, 1]
# LINESTYLES = [ ':','-.', '--','-',  ]


# for  ind, Ntot in enumerate(ListNtot):
#     fexplIterated = CalculateFractionPriorIterated_imperfect(rates, Ntot,  Vimperfect)
#     axarr.plot(rates, fexplIterated, c='k', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100, linestyle = LINESTYLES[ind])
#     # axarr.plot(rates, fexplIterated, c='b', lw = 6, alpha = ALPHA[ind])


# axarr.grid(True, alpha = 0.5)
# axarr.tick_params(labelsize=fs)
# axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs, fontweight = 'bold') #$ \mathcal{R}_{\mathrm{x}}$
# axarr.set_ylabel(r'\textbf{Fraction exploratory phase} $f_{\rm{expl}}$', fontsize=fs)
# axarr.tick_params(labelsize=fs, pad = 10.1) 

# axarr.text(0.000012, 0.2, r"\textbf{more rare}" +" \n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left', weight = 'bold' )
# axarr.text(0.051, 0.2,    r"\textbf{more common}" + "\n "  +r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
# axarr.legend(loc = 'upper right',  fontsize = fs, bbox_transform=axarr.transData) #, loc = 'upper center') #frameon=False, 

# # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
# plt.xscale('log')

# for i in range(6):
#     # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
#     axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [1, 1], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
#     axarr.text((ratess[i]), 0.015,  "  " + labelss[i] + " " , fontsize = fs-10, rotation=90, va = 'bottom', ha = 'center', zorder = 2)  
#     axarr.scatter(ratess[i], fexpls[i],  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 500,  linewidth=4.5 )


# plt.tight_layout()

# plt.savefig('AnalyticalFpriorIterated_imperfect.pdf')
# plt.savefig('AnalyticalFpriorIterated_imperfect.png', dpi = 300)




######### PLOT 2 ##########


# f, (axarr) = plt.subplots(1,1, figsize=(16,10))

# for  ind, Ntot in enumerate(ListNtot):
#     fexpl = CalculateFractionPriorIterated_imperfect(rates, Ntot,  Vimperfect)


#     NhitsMC = (rates*Ntot)
#     NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


#     ratioNhits = NhitsAIS / NhitsMC

#     axarr.plot(rates, ratioNhits, c='k', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100, linestyle = LINESTYLES[ind])
    





# axarr.grid(True, alpha = 0.5)

# axarr.tick_params(labelsize=fs)

# axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs) #$ \mathcal{R}_{\mathrm{x}}$
# axarr.set_ylabel(r'\textbf{Max. possible gain}', fontsize=fs)

# axarr.tick_params(labelsize=fs, pad = 10.1) 

# axarr.text(0.000012,20.3, r"\textbf{more rare}" + "\n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left' , weight ='bold')
# axarr.text(0.051, 20.3,   r"\textbf{more common}" +" \n " + r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
# axarr.legend( loc = 'upper right',  fontsize = fs, bbox_transform=axarr.transData) #, loc = 'upper center') #frameon=False, bbox_to_anchor=(0.1,3000)

# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
# plt.xscale('log')

# plt.yscale('log')

# for i in range(6):
#     # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
#     axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [10**4, 10**4], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
#     axarr.text((ratess[i]), 0.9* 10**4 ,  labelss[i] + "   " , fontsize = fs-10, rotation=90, va = 'bottom', ha = 'center', zorder = 2)


#     NhitsMCtemp = (ratess[i]*1E6)
#     NhitsAIStemp = (ratess[i]*1E6*fexpls[i] + (1-fexpls[i])*1E6) 


#     ratioNhitstemp = NhitsAIStemp / NhitsMCtemp

#     axarr.scatter(ratess[i], ratioNhitstemp,  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 1000,  linewidth=4.5 ) 





# plt.ylim(10**1, 10**4)
# plt.xlim(10**-5, 10**-1)
# plt.tight_layout()
# plt.savefig('AnalyticalFpriorRatioIterated_imperfect.pdf')
# plt.savefig('AnalyticalFpriorRatioIterated_imperfect.png', dpi = 300)








# # 
# ######### PLOT 3c: UNCERTAINTY with N = 1E4 AND N = 1E6 ##########
# #  this is the bottom plot of Figure 2 in STROOPWAFEL paper 

# f, (axarr) = plt.subplots(1,1, figsize=(16,10))

# for  ind, Ntot in enumerate([1E4]):   #(ListNtot):
#     fexpl = CalculateFractionPriorIterated_imperfect(rates, Ntot,  Vimperfect)
#     ind = 0

#     NhitsMC = (rates*Ntot)
#     NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


#     uncertaintyMC = 1 / np.sqrt(NhitsMC)
#     uncertaintyAIS = (1 / np.sqrt(NhitsAIS)) + ((1 / (fexpl*Ntot))/rates) 

#     uncertaintyRatio = uncertaintyAIS / uncertaintyMC

#     axarr.plot(rates, uncertaintyAIS   , c='gray', lw = 6,  zorder = 100, linestyle = LINESTYLES[3])
#     axarr.plot(rates, uncertaintyMC    , c='gray', lw = 6, zorder = 100, linestyle = LINESTYLES[0])
#     axarr.fill_between(rates,  uncertaintyAIS, uncertaintyMC, color='gray',  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)

# for  ind, Ntot in enumerate([1E6]):   #(ListNtot):
#     fexpl = CalculateFractionPriorIterated_imperfect(rates, Ntot,  Vimperfect)
#     ind = 2

#     NhitsMC = (rates*Ntot)
#     NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


#     uncertaintyMC = 1 / np.sqrt(NhitsMC)
#     uncertaintyAIS = (1 / np.sqrt(NhitsAIS)) #+ ((1 / (fexpl*Ntot))/rates) 

#     uncertaintyRatio = uncertaintyAIS / uncertaintyMC

#     axarr.plot(rates, uncertaintyAIS   , c='gray', lw = 6, label = 'this study',  zorder = 100, linestyle = LINESTYLES[3])
#     axarr.plot(rates, uncertaintyMC    , c='gray', lw = 6, label = 'traditional', zorder = 100, linestyle = LINESTYLES[0])
#     axarr.fill_between(rates,  uncertaintyAIS, uncertaintyMC, color='gray',  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)


# axarr.text(0.000018,0.065, r"\textbf{traditional uncertainty}", fontsize = fs-6,rotation=338, va = 'bottom', ha = 'left' )
# axarr.text(0.00002,0.0083, r"\textbf{STROOPWAFEL}" + "\n" + r"\textbf{uncertainty}", fontsize = fs-6,rotation=340, va = 'bottom', ha = 'left' )


# axarr.grid(True, alpha = 0.5)

# # axarr.pines()

# axarr.tick_params(labelsize=fs)

# axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs) #$ \mathcal{R}_{\mathrm{x}}$


# axarr.set_ylabel(r'\textbf{Fractional statistical uncertainty}', fontsize=fs)

# #cbar.ax
# axarr.tick_params(labelsize=fs, pad = 10.1) 

# axarr.text(0.000012,0.0017203, r"\textbf{more rare}" +" \n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left' )
# axarr.text(0.051, 0.0017203,   r"\textbf{more common}" +" \n " + r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
# # axarr.legend( loc = 'upper left',  fontsize = fs, bbox_transform=axarr.transData, ncol = 1) #, loc = 'upper center') #frameon=False, bbox_to_anchor=(0.1,3000)

# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
# plt.xscale('log')

# plt.yscale('log')

# for i in range(6):

#     # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
#     axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [10**4, 10**4], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
#     axarr.text((ratess[i]), 0.9*10**0,  labelss[i] + "   " , fontsize = fs-10, rotation=90, va = 'top', ha = 'center', zorder = 2)
#     axarr.scatter(ratess[i], FracUnc[i],  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 500,  linewidth=4.5 )


# axarr.text(0.009,0.033, r"$N = 10^4$", color = 'k',  fontsize = fs, va = 'bottom', ha = 'left' )    
# axarr.text(0.0107,0.002, r"$N = 10^6$", color = 'k',  fontsize = fs, va = 'bottom', ha = 'left' )


# plt.xlim(10**-5, 10**-1)
# plt.ylim(10**-3, 10**0)
# plt.tight_layout()
# plt.savefig('AnalyticalFpriorIteratedUncertaintyCOMBINED_imperfect.pdf')
# plt.savefig('AnalyticalFpriorIteratedUncertaintyCOMBINED_imperfect.png', dpi = 300)




######### PLOT 4 ##########
##  This is the Maximum gain VS Rate of target population plot (Fig 3 in STROOPWAFEL paper)

# LABELS2 = [r'$1$', r'$1/3$', r'$0.1$', r'$1/20$']

# f, (axarr) = plt.subplots(1,1, figsize=(16,10))

# for  ind, Ntot in enumerate(ListNtot):
#     if ind ==1:

#         fexpl = CalculateFractionPriorIterated_imperfect(rates, Ntot,  Vimperfect)
#         # print 'fexpl = ', fexpl
#         # NhitsMC = (rates*Ntot)
#         # NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot)
#         # Nexpl = np.zeros_like(rates)
#         # maskNexplTooLow = (rates*Ntot*fexpl < 4)

#         # Nexpl[maskNexplTooLow] = 4 / rates[[maskNexplTooLow]] 

#         # makszero = (Ntot - Nexpl) <= 0
#         # NhitsAIS[makszero] = np.zeros_like(NhitsAIS)[makszero]

#         # maksLeft = ((Ntot - Nexpl) > 0) & (Nexpl >0)
#         # NhitsAIS[maksLeft] = (Ntot - Nexpl )[maksLeft]



#         NhitsMC = (rates*Ntot)
#         NhitsAIS = (rates*Ntot*fexpl + (1-fexpl)*Ntot) 


#         ratioNhits = NhitsAIS / NhitsMC

#         axarr.plot(rates,      ratioNhits, c='r', lw = 6, label = LABELS2[0], alpha = ALPHA[ind], zorder = 100)
#         # axarr.plot(rates, 0.33* ratioNhits, c='r', lw = 6, label = LABELS2[1], alpha = ALPHA[ind], zorder = 100, linestyle = ':')
#         axarr.plot(rates, 0.1* ratioNhits, c='r', lw = 6, label = LABELS2[2], alpha = ALPHA[ind], zorder = 100, linestyle = '--')
#         # axarr.plot(rates, 0.05*ratioNhits, c='r', lw = 6, label = LABELS2[3], alpha = ALPHA[ind], zorder = 100, linestyle = '-.')
#         # axarr.plot(rates, NhitsAIS, c='r', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100)
#         # axarr.plot(rates, NhitsMC,  c='k', lw = 6, label = LABELS[ind], alpha = ALPHA[ind], zorder = 100)




# axarr.grid(True, alpha = 0.5)

# # axarr.pines()

# axarr.tick_params(labelsize=fs)

# axarr.set_xlabel(r'\textbf{Rate of target population}', fontsize=fs) #$ \mathcal{R}_{\mathrm{x}}$


# axarr.set_ylabel(r'\textbf{Gain} ', fontsize=fs)

# #cbar.ax
# axarr.tick_params(labelsize=fs, pad = 10.1) 

# axarr.text(0.000012,20.3, r"\textbf{more rare}" +" \n " + r"\textbf{events}", fontsize = fs-6,rotation=90, va = 'bottom', ha = 'left' )
# axarr.text(0.051, 20.3,   r"\textbf{more common}" +" \n " + r"\textbf{events}", fontsize = fs-6, rotation=90, va = 'bottom', ha = 'left')
# lg = axarr.legend( loc = 'upper left',  bbox_transform=axarr.transData,  title="title", fontsize = fs-5, ncol = 2 ) #, loc = 'upper center') #frameon=False, bbox_to_anchor=(0.1,3000)



# axarr.text(0.07, 6.3*1E3, r"$N = 10^6$", fontsize = fs+15, rotation=0, va = 'top', ha = 'right' )

# lg.set_title(r'\textbf{Efficiency of}' + '\n' + r'\textbf{refinement phase}')
# lg.get_title().set_fontsize(fs-5)

# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), fs = fs)
# plt.xscale('log')

# plt.yscale('log')

# for i in range(6):
#     # mpl.patches.Rectangle((ratess[i], 10**2), width, height, angle=0.0, **kwargs)
#     axarr.fill_between([(ratess[i])-0.100*(ratess[i]), (ratess[i])+0.100*(ratess[i])],  [10**4, 10**4], color=colors[i],  alpha = 0.55, edgecolor="k", linewidth=0.2, zorder=1)
#     axarr.text((ratess[i]), 0.9* 10**4,  labelss[i] + "   " , fontsize = fs-10, rotation=90, va = 'top', ha = 'center', zorder = 2)
#     axarr.scatter(ratess[i], Gains[i],  s = 730, marker = '*', edgecolors = 'k', color = colors[i], zorder = 1000,  linewidth=4.5 ) 

# plt.ylim(10**1, 10**4)
# plt.xlim(10**-5, 10**-1)
# plt.tight_layout()
# plt.savefig('AnalyticalFpriorRatioIteratedVaryEfficiency_imperfect.pdf')
# plt.savefig('AnalyticalFpriorRatioIteratedVaryEfficiency_imperfect.png', dpi = 300)















