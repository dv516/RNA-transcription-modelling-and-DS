# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:52:00 2020

@author: dv516
"""

import numpy as np
import time

import matplotlib as mpl
import matplotlib.pyplot as plt

import pylab as pl

import warnings

from odes_and_curve_fitting_functions import datafitting_transcription_prob_stDev

## These are our optimal parameters. The reader should also ave gotten similar
## values to these in curve_fitting_exploration_predicton. If the optimal 
## parameters are different due to different dependencies (which shouldn't 
## happen), the reader may wish to change these optimal parameters
## These parameters can also be changed for exploration
opt_param = [4.34, 5.55e+05, 1.94e+05, 1.20e+06, 0, 0, 0, 0, 0]

new = opt_param

perr = list(np.array(new) * 0.2) # Assume 20% standard deviation

## Probabilistic design space
## Load simulations, otherwise: Will take 3 hrs to run. Uncomment to run simulations.


# Create 2d grid at constant T7RNAP
N_Mg = 51
N_NTP = 41
grid_Mg3d = np.linspace(0.01, 0.14, N_Mg)
grid_NTP3d = np.linspace(0.01, 0.10, N_NTP)
X_prob = np.zeros((N_Mg*N_NTP, 5))
X_prob[:,0] = 6 # Time fixed to 6 hours
X_prob[:,3] = 0.0002 # Constant spermidine concentration
X_prob[:, 4] = 1
N_total = 0
for i in range(0, N_Mg):
    for j in range(0, N_NTP):
        X_prob[N_total, 1] = grid_Mg3d[i]
        X_prob[N_total, 2] = grid_NTP3d[j]
        N_total += 1

N = 50 # 50 Monte Carlo simulations over grid        
        
proba = np.zeros(len(X_prob))
y_sim_prob = np.zeros((N, len(proba)))
        
threshold = 1.5 # The RNA yield CQA threshold

warnings.filterwarnings("ignore", category= RuntimeWarning)

### EITHER: Run simulation
# t0 = time.time()
# for k in range(N): # for each simulation, randomly pick parameters from distribution                 
#     y_prob1 = datafitting_transcription_prob_stDev(X_prob, perr, *new) # Run and time simulation
#     t1 = time.time()   
#     y_sim_prob[k] = y_prob1 # Store yield in matrix
#     for j in range(len(proba)):
#         if y_prob1[j] >= threshold: # for each grid point, if it meets the threshold, add 1 to proba
#             proba[j] += 1           
#     print('Simulation ', k+1, ' complete: ',int((k+1)/N*100), '%')
#     print('Time for ' + str((k+1)*len(y_prob1)) + ' simulations: ', int((t1-t0)/6)/10, ' min')
#     print('Estimated waiting time: ', int((t1-t0)/(k+1)*N/60 - (t1-t0)/60), ' min')   
# np.savetxt('y_sim_prob_Design.csv', y_sim_prob, delimiter=",")
# np.savetxt('prob_Design.csv', proba, delimiter=",")
###

### OR: Load previous simulation results
proba = np.loadtxt('prob_Design.csv', delimiter=",")
###

fig = plt.figure(figsize = (6, 10))
ax5 = fig.add_subplot(2, 1, 1)

pl.rc('axes', linewidth=2) # setting the width of the figure border lines thicker (2)

cm = plt.get_cmap("RdYlBu") # Colour map

for axis in ['top','right']:
  ax5.spines[axis].set_linewidth(1)  # Setting the width of the top and right border lines to lower thickness (1)

#col1 = []

proba_tot = proba # + proba2 + proba3 + proba4 + proba5

for i in range(len(proba_tot)):
    if (proba_tot[i] == 0):
        col = (1, 1, 1)
    else:
        col = cm(proba_tot[i]/N)
    ax5.scatter(X_prob[i,1], X_prob[i,2], color = col, cmap = cm, s=100, alpha = 0.85, marker='o')

ax5.scatter(0.085, 0.04, color = 'k', marker = 'X', s = 150)

ax5.set_xlim(0, 0.14)
ax5.set_ylim(0, 0.10)
#ax5.set_ylim3d(0, 0.1)
#ax5.set_xlim3d(0, 0.16)


##### Code for adjusting the figure appearance
tick_x_pDS = np.arange(0, 0.14, 5)
formatter_x_pDS = mpl.ticker.StrMethodFormatter('{x:1.2f}')
ax5.xaxis.set_ticks(tick_x_pDS)
ax5.xaxis.set_major_formatter(formatter_x_pDS)
#ax2.xaxis.set_tick_params(labelsize=14, width=2) 

tick_y_pDS = np.arange(0, 0.1, 6)
#formatter_y_pDS = mpl.ticker.StrMethodFormatter('{x:1.2f}')
#ax5.yaxis.set_ticks(tick_y_pDS)
#ax5.yaxis.set_major_formatter(formatter_y_pDS)
#ax5.yaxis.set_tick_params(labelsize=14, labelrotation=0, width=2)

ax5.xaxis.set_tick_params(labelsize=16, width=2) 
ax5.yaxis.set_tick_params(labelsize=16, width=2) 
#####



ax5.set_xlabel('Total Mg conc. [M]', fontsize = 16, labelpad = 16)
ax5.set_xticks(np.linspace(0, 0.12, 5))
ax5.set_yticks(np.linspace(0, 0.1, 6))
ax5.set_ylabel('Initial NTP conc. [M]', fontsize = 16, labelpad = 16)
#ax5.set_title('2d - Probabilistic design space - Batch - Initial', fontsize = 14, fontweight = "bold")

plt.show()

#fig.savefig("2dProbDesSpace.svg")


fig, ax = plt.subplots(figsize=(8, 1))
fig.subplots_adjust(top=0.5)

cmap =  plt.get_cmap("RdYlBu") # Create colourmap
norm = mpl.colors.Normalize(vmin=0, vmax= 1)

cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Probability of achieving 1.5 g/L effective RNA yield CQA', fontsize = 14, fontweight = "bold", labelpad = 10)
cb1.set_ticks(np.linspace(0, 1, 6))
cb1.ax.tick_params(labelsize=16, width=2) 
fig.show()
#fig.savefig('Colorbar_probDesSpace.svg', bbox_inches = 'tight')
