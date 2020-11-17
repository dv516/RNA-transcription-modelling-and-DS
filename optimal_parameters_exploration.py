# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:38:42 2020

@author: dv516
"""

import numpy as np
import time

import matplotlib.pyplot as plt
import pylab as pl

import warnings

from odes_and_curve_fitting_functions import datafitting_transcription_experimental
from data_functions import get_experimental_data

X_new, C_RNA_new, C_RNA_stdev, MAE_exp, MSE_exp = get_experimental_data()

# opt_param = [10.9, 1.4e6, 4.89e5, 1.20e6, 0, 0, 0, 0, 0] # If run without optimal model parameters
opt_param = [4.34, 5.55e+05, 1.94e+05, 1.20e+06, 0, 0, 0, 0, 0]

new = opt_param

## Plot dependence of RNA yield on Mg

plt.figure(figsize=(8, 6))
pl.rc('axes', linewidth=2) # setting the width of the figure border lines thicker (2)
warnings.filterwarnings("ignore", category= RuntimeWarning)
## Ignore Runtime warnings about fsolve() convergence being slow, as initguess() etc. ignores 
## 'bad' solutions anyway

# Set input matrix X: 30 samples with varying Mg at constant spermidine, T7RNAP and NTP after 2 hrs
N = 30
X_fitting = np.zeros((N, 5))
Mg_linSp = np.linspace( 0.025, 0.125, N)
X_fitting[:,3] = 0.0002 # Constant spermidine concentration
X_fitting[:, 4] = 1
X_fitting[:, 2] = 0.04
X_fitting[:,0] = 2
X_fitting[:, 1] = Mg_linSp
plt.plot(X_fitting[:, 1], datafitting_transcription_experimental(X_fitting, *new),'--b',          label = 'Simulation: 2 hrs') # Plot 2 hr dependence

X_fitting[:, 0] = 4 # Same after 4 hrs
plt.plot(X_fitting[:, 1], datafitting_transcription_experimental(X_fitting, *new), '--g',          label = 'Simulation: 4 hrs')

X_fitting[:, 0] = 6 # Same after 6 hrs
plt.plot(X_fitting[:,1], datafitting_transcription_experimental(X_fitting, *new), '--m',          label = 'Simulation: 6 hrs') # Same after 6 hrs

# Plot experimental data plus 1 standard deviation
plt.errorbar(X_new[:11, 1], C_RNA_new[:11], C_RNA_stdev[:11], \
             label = 'Experimental after 2 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'b')
plt.errorbar(X_new[11:11*2, 1], C_RNA_new[11:11*2], C_RNA_stdev[11:11*2], \
             label = 'Experimental after 4 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'g')
plt.errorbar(X_new[11*2:11*3, 1], C_RNA_new[11*2:11*3], C_RNA_stdev[11*2:11*3], \
             label = 'Experimental after 6 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'm')

plt.ylim(0, 3.5)
plt.xlim(0.01, 0.13)

#plt.title("Data fitting of RNA yield dependence on Mg \n concentration at different reaction times", \
#          fontsize = 14, fontweight = "bold")
plt.xlabel('$[Mg]_{total}$ [M]', fontsize = 12)
plt.ylabel('$[RNA]_{effective}$ [g/L]', fontsize = 12)

plt.legend(loc = "upper left");


## Plot dependence of RNA yield on T7RNAP

plt.figure(figsize=(8, 6))

# Set input matrix X: 15 samples with varying T7RNAP at constant spermidine, Mg and NTP after 2 hrs
N = 15
X_fitting = np.zeros((N, 5))
Tr_linSp = np.linspace( 0.125, 1, N)
X_fitting[:,3] = 0.0002 # Constant spermidine concentration
X_fitting[:, 4] = Tr_linSp
X_fitting[:, 2] = 0.04
X_fitting[:,0] =  2
X_fitting[:, 1] = 0.075
plt.plot(X_fitting[:, 4]*5e-7, datafitting_transcription_experimental(X_fitting, *new),'--b',          label = 'Simulation: 2 hrs') # Plot

X_fitting[:, 0] = 4 # Same after 4 hrs
plt.plot(X_fitting[:, 4]*5e-7, datafitting_transcription_experimental(X_fitting, *new), '--g',          label = 'Simulation: 4 hrs')

X_fitting[:, 0] = 6 # Same after 6 hrs
plt.plot(X_fitting[:,4]*5e-7, datafitting_transcription_experimental(X_fitting, *new), '--m',          label = 'Simulation: 6 hrs')

plt.errorbar(X_new[33:33+4, 4]*5e-7, C_RNA_new[33:33+4], C_RNA_stdev[33:33+4], \
             label = 'T7RNAP dependence after 2 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'b')
plt.errorbar(X_new[33+4:33+4*2, 4]*5e-7, C_RNA_new[33+4:33+4*2], C_RNA_stdev[33+4:33+4*2], \
             label = 'T7RNAP dependence after 4 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'g')
plt.errorbar(X_new[33+4*2:33+4*3, 4]*5e-7, C_RNA_new[33+4*2:33+4*3], C_RNA_stdev[33+4*2:33+4*3], \
             label = 'T7RNAP dependence after 6 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'm')

plt.ylim(0, 3.5)
plt.xlim(0.01*5e-7, 1.25*5e-7)

#plt.title("Data fitting of RNA yield dependence on \n [T7RNAP] at different reaction times", fontsize = 14, fontweight = "bold")
plt.xlabel('$[T7RNAP]$ [M]', fontsize = 12)
plt.ylabel('$[RNA]_{effective}$ [g/L]', fontsize = 12)

plt.legend(loc = "upper left");



## Plot dependence of RNA yield on NTP at medium and high Mg

plt.figure(figsize=(8, 6))

# Set input matrix X: 20 samples with varying NTP at constant spermidine, Mg and T7RNAP after 2 hrs
N = 20
X_fitting = np.zeros((N, 5))
NTP_linSp = np.linspace( 0.02, 0.08, N)
X_fitting[:,3] = 0.0002 # Constant spermidine concentration
X_fitting[:, 4] = 1
X_fitting[:, 2] = NTP_linSp
X_fitting[:,0] =  2
X_fitting[:, 1] = 0.075
plt.plot(X_fitting[:, 2], datafitting_transcription_experimental(X_fitting, *new),'--b',          label = 'Simulation: 75 mM Mg2+')

X_fitting[:, 1] = 0.14 # Same with higher Mg
plt.plot(X_fitting[:, 2], datafitting_transcription_experimental(X_fitting, *new), '--g',          label = 'Simulation: 140 mM Mg2+')

plt.errorbar(X_new[-6:-3, 2], C_RNA_new[-6:-3], C_RNA_stdev[-6:-3], \
             label = 'NTP dependence at low Mg after 2 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'b')
plt.errorbar(X_new[-3:, 2], C_RNA_new[-3:], C_RNA_stdev[-3:], \
             label = 'NTP dependence at low Mg after 2 hrs', linestyle = '', \
             marker = 'o', markersize = 3, color = 'g')

plt.ylim(0, 2)
plt.xlim(0.01, 0.09)

#plt.title("Data validation of RNA yield dependence on \n [NTP] at different [Mg] after 2 hrs", \
#          fontsize = 14, fontweight = "bold")
plt.xlabel('$[NTP_{total}]$ [M]', fontsize = 12)
plt.ylabel('$[RNA]_{effective}$ [g/L]', fontsize = 12)

plt.legend(loc = "upper left");


## 3d plot of RNA yield dependence on Mg and NTP

N_Mg = 31
N_NTP = 31

## Set 2d grid of Mg and NTP dependence
grid_Mg = np.linspace(0.01, 0.15, N_Mg)
grid_NTP = np.linspace(0.01, 0.15, N_NTP)
X_2d = np.zeros((N_Mg*N_NTP, 5))
X_2d[:,0] = 6 # Time fixed to 6 hours
X_2d[:,3] = 0.0002 # Constant spermidine concentration
X_2d[:,4] = 1

N_total = 0
for i in range(0, N_Mg):
    for j in range(0, N_NTP):
        X_2d[N_total, 1] = grid_Mg[i]
        X_2d[N_total, 2] = grid_NTP[j]
        N_total += 1
        
# Run and record simulation time      
t0 = time.time()
y_2d = datafitting_transcription_experimental(X_2d, *new)
t1 = time.time()
print("Computation time for " + str(len(X_2d)) + " simulations: ", int(t1 - t0), " seconds")

fig1 = plt.figure(figsize = (8, 6))
ax1 = fig1.add_subplot(111, projection = '3d')

# Colour code such that high RNA yield is green, and low RNA yield is yellow
for z in range(len(y_2d)):
    col = (1-max(y_2d[z], 0)/6,0.7, 0)
    ax1.scatter(X_2d[z,1], X_2d[z,2], y_2d[z], color = col, s = 100)

ax1.set_xlabel('Initial total Mg concentration [mol/L]')
ax1.set_ylabel('Initial total NTP concentration [mol/L]')
ax1.set_zlabel('Effective RNA yield CQA [g/L]');

#plt.xlim(0, 0.14)
#plt.ylim(0, 0.1)