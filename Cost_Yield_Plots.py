# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:49:09 2020

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

## Cost per yield plot

fig = plt.figure(figsize = (8, 6))
ax = fig.add_subplot(111, projection='3d')
cm = plt.get_cmap("copper") # Colour map

cost_component = {'T7RNAP': 13474567795.68, 'Mg': 754.20, 'NTP': 250000.00} # cost per mol of components


# Construct the NTP and T7RNAP grid for constant Mg (Magnesium is inexpensive)
N_NTP = 51
N_T7RNAP = 51
grid_NTP = np.linspace(0.01, 0.08, N_NTP)
grid_T7RNAP = np.linspace(0.5, 1.5, N_T7RNAP)
X_cost = np.zeros((N_T7RNAP*N_NTP, 5))
X_cost[:,0] = 6 # Time fixed to 6 hours
X_cost[:,3] = 0.0002 # Constant spermidine concentration
X_cost[:,1] = 0.085
N_total = 0
for j in range(0, N_NTP):
    for k in range(0,N_T7RNAP):
        X_cost[N_total, 2] = grid_NTP[j]
        X_cost[N_total, 4] = grid_T7RNAP[k]
        N_total += 1

warnings.filterwarnings("ignore", category= RuntimeWarning)

## Simulate and time
t0 = time.time()
y_yield = datafitting_transcription_experimental(X_cost, *new)
t1 = time.time()
print("Computation time for " + str(len(y_yield)) + " simulations: ", int((t1 - t0)), " sec")        
        
cost = lambda NTP, T7RNAP: NTP*cost_component['NTP']+ T7RNAP*cost_component['T7RNAP']*1e-8 # cost function
costPerYield = cost(X_cost[:,2], X_cost[:,4]) / y_yield

threshold = 25000 # upper limit to be shown (as yield approaches 0, costPerYield approaches infinity)

for z in range(0, len(y_yield)): # colour code and plot only yield is below threshold
    y = costPerYield[z]
    if y < threshold:
        col = cm(y_yield[z]/max(y_yield))
        ax.scatter(X_cost[z,2], X_cost[z,4]*1e-8, y, color = col, s = 30)  
   
## Find and plot minimum costperYield
ind = np.where(costPerYield == min(costPerYield))
ax.scatter(X_cost[ind, 2], X_cost[ind, 4]*1e-8, min(costPerYield), color = 'k', s = 140)

##### Code for adjusting the figure appearance
tick_x = np.arange(0.01, 0.07, 0.04)
formatter_x = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax.xaxis.set_ticks(tick_x)
ax.xaxis.set_major_formatter(formatter_x)
ax.xaxis.set_tick_params(labelsize=15, width=4) 

#tick_y = np.arange(0.5e-8, 1.5e-8, 0.5e-8)
#formatter_y = mpl.ticker.StrMethodFormatter('{x:1.2f}')
#ax.yaxis.set_ticks(tick_y)
#ax.yaxis.set_major_formatter(formatter_y)
ax.yaxis.set_tick_params(labelsize=15, width=4, labelrotation=0)
#for label in ax.yaxis.get_ticklabels():
#    label.set_horizontalalignment('center')

#tick_z = np.arange(1000, threshold)
formatter_z = mpl.ticker.StrMethodFormatter('{x:1.0f}')
#ax.zaxis.set_ticks(tick_z)
ax.zaxis.set_major_formatter(formatter_z)
ax.zaxis.set_tick_params(labelsize=15, pad=8, width=4) 
#####

ax.set_zlim3d(1000, threshold)
ax.set_xlabel('Initial NTP conc. [M]', fontsize = 15, labelpad = 15)
ax.set_xticks(np.linspace(0.01, 0.07, 4))
ax.set_ylabel('T7RNAP conc. x $10^{-8}$ [M]', fontsize = 15, labelpad = 19)
ax.set_yticks(np.linspace(0.5e-8, 1.6e-8, 4))
ax.set_zlabel('Cost of T7RNAP & NTPs \n per g of RNA [US$/g]', fontsize = 15, labelpad = 22)


ax.view_init(8, 142)

#fig.savefig("costPerYield.svg", bbox_inches = 'tight')


fig, ax = plt.subplots(figsize=(6.3, 1))
fig.subplots_adjust(top=0.5)

pl.rc('axes', linewidth=2) # setting the width of the figure border lines thicker (2)

cmap =  plt.get_cmap("copper") # Create colourmap which is the reverse of the copper color map used in previous cell
norm = mpl.colors.Normalize(vmin= 0, vmax= 4.5)

cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Effective RNA yield CQA [g/L]', fontsize = 15, fontweight = "bold", labelpad = 10)
cb1.set_ticks(np.linspace(0, 4.5, 4))
cb1.ax.tick_params(labelsize=15, width=2) 
fig.show()
#fig.savefig('Colorbar_CostPerYield.svg', bbox_inches = 'tight')


print('Optimal NTP: ', X_cost[ind, 2], 'mol/L')
print('Optimal T7RNAP: ', X_cost[ind, 4], 'U/microL')
print('Yield at optimal cost point: %.2f' % y_yield[ind], " g/L")
print('Maximum yield: %.2f' % max(y_yield), " g/L")
print('Optimal cost per yield: ', costPerYield[ind], "$ / g RNA")
NTP_cont = int(X_cost[ind, 2]*cost_component['NTP']                / (X_cost[ind, 2]*cost_component['NTP']+ X_cost[ind, 4]*cost_component['T7RNAP']*1e-8) * 100)
print('Contibution of NTP: ', int(X_cost[ind, 2]*cost_component['NTP'] / y_yield[ind]), " $ / g RNA or ", NTP_cont, " %")
