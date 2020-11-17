# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:15:52 2020

@author: dv516
"""

import numpy as np
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl

from scipy.optimize import curve_fit

import warnings

from odes_and_curve_fitting_functions import datafitting_transcription_experimental

from data_functions import plot_experimental_data, get_experimental_data

X_new, C_RNA_new, C_RNA_stdev, MAE_exp, MSE_exp = get_experimental_data()
        
# Plot experimental data set
fig = plt.figure(figsize = (7, 4))
ax = fig.add_subplot(111)      
ax = plot_experimental_data(ax, C_RNA_new, C_RNA_stdev)
ax.set_ylabel('Effective RNA yield CQA [g/L]')
ax.set_xlabel('Index of dataset')
# Print experimental errors for reference
print('Experimental mean absolute error: %.4f' % np.average(MAE_exp))
print('Experimental mean squared error: %.4f' % np.average(MSE_exp))

## Ignore Runtime warnings about fsolve() convergence being slow, 
## as initguess() ignores 'bad' solutions anyway
warnings.filterwarnings("ignore", category= RuntimeWarning)


## Parameter estimation/ Curve fitting: 1st entry is funcction, 2nd is input, 3rd is y value to fit to
## p0 is initial guess, and bounds are parameter bounds (bounds could also be specified for each parameter individually)
pexp, pcov__ = curve_fit(datafitting_transcription_experimental, X_new, C_RNA_new, \
                          p0 = [13000*1e-7, 2e1, 100, 1e6, 1e6, 2, 0, 0, 0], \
                          bounds = (0, 1e8))

print('Optimal parameters: ', pexp)

StDev = np.sqrt(np.diag(pcov__))
print('Standard deviation: ', StDev)

outer_stdev = np.outer(StDev, StDev)
corr = pcov__ / outer_stdev
corr[pcov__ == 0] = 0
print('Correlation matrix: ', corr)

## Entries in arrays correspond to kapp, K1, K2, k_ac, k_ba, k_Mg, K3, K4, K5
## Correlations matrix works similarly, for example, 2nd row, first element gives covariance between kapp and K1


## Test optimal parameters

fig = plt.figure(figsize = (14, 11))
ax = fig.add_subplot(211)

ax = plot_experimental_data(ax, C_RNA_new, C_RNA_stdev)

#rcParams['axes.linewidth'] = 5 # set the value globally
# Thicken the axes lines and labels
pl.rc('axes', linewidth=2) # setting the width of the figure border lines thicker (2)


# Set list of colour of each sample point
color_list = [plots.get_color()  for plots in ax.get_lines()] # Get colors of 8 different experimental regions
col = [None] * len(C_RNA_new) # Allocate size
for i in range(11):
    col[i] = color_list[0] ; col[11+i] = color_list[1] ; col[11*2+i] = color_list[2]
for i in range(4):
    col[33+i] = color_list[3] ; col[33+4 + i] = color_list[4] ; col[33+4*2 + i] = color_list[5] 
for i in range(3):
    col[-6+i] = color_list[6] ; col[-3+i] = color_list[7]

# Run simulation and record duration
t0 = time.time()
y_sim = datafitting_transcription_experimental(X_new, *pexp)
t1 = time.time()
print(t1 - t0, 'seconds for ', len(X_new), 'simulations')

# Compare fit using pexp with experimental data set
ax.scatter(np.arange(len(y_sim)), y_sim, c = col, marker = 'x')
ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0., fontsize = 16, markerscale = 5.1)
ax.set_xlabel('Index of experimental dataset')
ax.set_ylabel('Effective RNA yield CQA [g/L]')

print('Simulated mean absolute error: ', np.average(abs(y_sim - C_RNA_new)))
print('Simulated mean squared error: ', np.average((y_sim - C_RNA_new)**2))

# Exploration space
ax1 = fig.add_subplot(224, projection = '3d')
scatter_size = np.zeros(len(C_RNA_new)) + 100
scatter_size[-6:] = 30 # Size of points is smaller if only sampled after 2 hrs rather than 2, 4 and 6 hrs
ax1.scatter(X_new[:,1], X_new[:,2], X_new[:,4], c=col, s = scatter_size)    
ax1.set_xlim3d(0.01, 0.14)
ax1.set_ylim3d(0.01, 0.16)
ax1.set_zlim3d(0.125, 1)

for axis in ['bottom','left']:
  ax1.spines[axis].set_linewidth(2) # Setting the width of the top and right border lines to lower thickness (1)

##### Code for adjusting the figure appearance
tick_x = np.arange(0.01, 0.14)
formatter_x = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax1.xaxis.set_ticks(tick_x)
ax1.xaxis.set_major_formatter(formatter_x)
ax1.xaxis.set_tick_params(labelsize=13, width=4) 

tick_y = np.arange(0.01, 0.16)
formatter_y = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax1.yaxis.set_ticks(tick_y)
ax1.yaxis.set_major_formatter(formatter_y)
ax1.yaxis.set_tick_params(labelsize=13, labelrotation=-20, width=4)
#for label in ax.yaxis.get_ticklabels():
#    label.set_horizontalalignment('center')

tick_z = np.arange(0.125, 1)
formatter_z = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax1.zaxis.set_ticks(tick_z)
ax1.zaxis.set_major_formatter(formatter_z)
ax1.zaxis.set_tick_params(labelsize=13, pad=8, width=4) 
#ax1.view_init(20, 145)
#####

ax1.set_xlabel('Total Mg conc. [M]', fontsize = 13, labelpad = 14)
ax1.set_xticks(np.linspace(0.01, 0.14, 4))
ax1.set_ylabel('Initial NTP conc. [M]', fontsize = 13, labelpad = 16)
ax1.set_yticks(np.linspace(0.03, 0.15, 3))
ax1.set_zlabel('T7RNAP conc. x $10^{-8}$ [M]', fontsize = 13, labelpad = 16)
ax1.set_zticks(np.linspace(0.2, 1, 4))


# Prediction error plot
ax2 = fig.add_subplot(223)
ax2.scatter(C_RNA_new, y_sim, c = col, s = scatter_size)
ax2.plot([0, max(max(C_RNA_new), max(y_sim))], [0, max(max(C_RNA_new), max(y_sim))], '-k', label = 'Identity line')

for axis in ['top','right']:
  ax2.spines[axis].set_linewidth(1)

###
#ax2.set_ylim(0, 0.10)
#ax2.set_xlim(0, 0.14)
##ax2.set_ylim3d(0, 0.1)
##ax2.set_xlim3d(0, 0.16)

##### Code for adjusting the figure appearance
tick_x_pDS = np.arange(0, 3.1, 1)
formatter_x_pDS = mpl.ticker.StrMethodFormatter('{x:1.0f}')
ax2.xaxis.set_ticks(tick_x_pDS)
ax2.xaxis.set_major_formatter(formatter_x_pDS)
#ax2.xaxis.set_tick_params(labelsize=14, width=2) 

tick_y_pDS = np.arange(0, 3.1, 1)
formatter_y_pDS = mpl.ticker.StrMethodFormatter('{x:1.0f}')
ax2.yaxis.set_ticks(tick_y_pDS)
ax2.yaxis.set_major_formatter(formatter_y_pDS)
#ax2.yaxis.set_tick_params(labelsize=14, labelrotation=0, width=2)

ax2.xaxis.set_tick_params(labelsize=18, width=2) 
ax2.yaxis.set_tick_params(labelsize=18, width=2) 
#####

#ax2.scatter(label = 'color_list')
ax2.legend(loc = 'upper left', fontsize = 16)
ax2.set_ylabel('Predicted RNA yield CQA [g/L]', fontsize = 18, labelpad = 16)
ax2.set_xlabel('Actual RNA yield CQA [g/L]', fontsize = 18, labelpad = 16);