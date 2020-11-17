# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:43:12 2020

@author: dv516
"""

import numpy as np
import time

import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl

import warnings

from odes_and_curve_fitting_functions import datafitting_transcription_experimental

# opt_param = [10.9, 1.4e6, 4.89e5, 1.20e6, 0, 0, 0, 0, 0] # If run without optimal model parameters
opt_param = [4.34, 5.55e+05, 1.94e+05, 1.20e+06, 0, 0, 0, 0, 0]

new = opt_param

## 3d deterministic design space
# Set 3d grid
N_Mg = 41
N_NTP = 41
N_T7RNAP = 11
grid_Mg = np.linspace(0.01, 0.15, N_Mg)
grid_NTP = np.linspace(0.01, 0.15, N_NTP)
grid_T7RNAP = np.linspace(0.5, 1.5, N_T7RNAP)
X_grid = np.zeros((N_Mg*N_NTP*N_T7RNAP, 5))
X_grid[:,0] = 6 # Time fixed to 6 hours
X_grid[:,3] = 0.0002 # Constant spermidine concentration
N_total = 0
for i in range(0, N_Mg):
    for j in range(0, N_NTP):
        for k in range(0,N_T7RNAP):
            X_grid[N_total, 1] = grid_Mg[i]
            X_grid[N_total, 2] = grid_NTP[j]
            X_grid[N_total, 4] = grid_T7RNAP[k]
            N_total += 1

## Ignore Runtime warnings about fsolve() convergence being slow, as initguess() etc. ignores 
## 'bad' solutions anyway
warnings.filterwarnings("ignore", category= RuntimeWarning)


### EITHER: Simulate and record simulation time
# t0 = time.time()
# y_designSpace = datafitting_transcription_experimental(X_grid, *new)
# t1 = time.time()
# print("Computation time for " + str(len(y_designSpace)) + " simulations: ", int((t1 - t0)/6)/10, " min")
# np.savetxt('3d_DesSpace_new.csv', y_designSpace, delimiter=",") # Save simulation to avoid running time
###

### OR: Load previous simulation results
y_designSpace = np.loadtxt('3d_DesSpace_new.csv', delimiter=",") # Load simulation to avoid running time
###

fig = plt.figure(figsize = (8, 6))
ax = fig.add_subplot(111, projection='3d')
pl.rc('axes', linewidth=2) # setting the width of the figure border lines thicker (2)

for z in range(len(y_designSpace)): # colour code with the higher the RNA yield, the greener
    # filter out solutions that don't meet design criteria or are unstable
    if (y_designSpace[z]>1.5) and (y_designSpace[z]<1e3):
        if y_designSpace[z]<5: 
            col = (1-(y_designSpace[z]-1.5)/(5-1.5),0.7, 0)
        else:
            col = (0,0.7,0)
        ax.scatter(X_grid[z,1], X_grid[z,2], X_grid[z,4], color = col, s = 50)    

##### Code for adjusting the figure appearance
tick_x = np.arange(0.00, 0.165, 0.04)
formatter_x = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax.xaxis.set_ticks(tick_x)
ax.xaxis.set_major_formatter(formatter_x)
ax.xaxis.set_tick_params(labelsize=15, width=4) 

tick_y = np.arange(0.00, 0.165, 0.04)
formatter_y = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax.yaxis.set_ticks(tick_y)
ax.yaxis.set_major_formatter(formatter_y)
ax.yaxis.set_tick_params(labelsize=15, width=4)
#for label in ax.yaxis.get_ticklabels():
#    label.set_horizontalalignment('center')

tick_z = np.arange(2.5, 7.62, 1.5)
formatter_z = mpl.ticker.StrMethodFormatter('{x:1.1f}')
ax.zaxis.set_ticks(tick_z)
ax.zaxis.set_major_formatter(formatter_z)
ax.zaxis.set_tick_params(labelsize=15, pad=8, width=4) 
#####


ax.set_xlim3d(0, 0.165)
ax.set_ylim3d(0, 0.165)
ax.set_zlim3d(0.4, 1.50)

ax.set_xlabel('Total Mg conc. [M]', fontsize = 15, labelpad = 15)
ax.set_xticks(np.linspace(0.02, 0.16, 4))
ax.set_ylabel('Initial NTP conc. [M]', fontsize = 15, labelpad = 16)
ax.set_yticks(np.linspace(0.02, 0.16, 4))
ax.set_zlabel('T7RNAP conc. x $10^{-8}$ [M]', fontsize = 15, labelpad = 14)
ax.set_zticks(np.linspace(0.5, 1.5, 4))
#ax.set_title('3d Deterministic Design space - Initial -  Batch \n CQA of > 1.5 g/L RNA', \
#             fontsize = 14, fontweight = "bold")

#ax.view_init(30, 120)
ax.view_init(15, 142)

## Create colourmap and colourbar

fig, ax = plt.subplots(figsize=(6.3, 1))
fig.subplots_adjust(top=0.5)

cmap = mpl.colors.ListedColormap([((60 - i)/(60 - 15), 0.7, 0) for i in range(15,61)]) # Create colourmap
norm = mpl.colors.Normalize(vmin=1.5, vmax=6)

cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Effective RNA yield CQA [g/L]', fontsize = 15, fontweight = "bold", labelpad = 10)
cb1.set_ticks(np.linspace(1.5, 6, 4))
cb1.ax.tick_params(labelsize=15, width=2) 
fig.show()
#fig.savefig('Colorbar_DesSp.svg', bbox_inches = 'tight')