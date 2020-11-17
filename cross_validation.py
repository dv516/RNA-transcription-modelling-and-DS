# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:27:46 2020

@author: dv516
"""

# K-fold cross validation

import numpy as np

import matplotlib.pyplot as plt
import pylab as pl

from scipy.optimize import curve_fit

import warnings

from odes_and_curve_fitting_functions import datafitting_transcription_experimental
from data_functions import plot_experimental_data, get_experimental_data

# pexp = [10.9, 1.4e6, 4.89e5, 1.20e6, 0, 0, 0, 0, 0] # If run without optimal model parameters
pexp = [4.34, 5.55e+05, 1.94e+05, 1.20e+06, 0, 0, 0, 0, 0]

X_new, C_RNA_new, C_RNA_stdev, MAE_exp, MSE_exp = get_experimental_data()

# Randomize and shuffle the dataset
np.random.seed(123)
X_copy = X_new.copy()
C_RNA_copy = C_RNA_new.copy()
idx = np.random.permutation(len(X_copy))
X_shuffled, Y_shuffled = X_copy[idx], C_RNA_copy[idx]


prop_train = 0.9 # 10 cross-fold
idx_train = int(prop_train*len(X_shuffled))
K = int(1/(1-prop_train))

param_matrix = np.zeros((K, len(pexp)))

## Ignore Runtime warnings about fsolve() convergence being slow, as initguess() etc. ignores 
## 'bad' solutions anyway
warnings.filterwarnings("ignore", category= RuntimeWarning)

for i in range(K):
    # for K iterations, set the first idx_train samples as training, and the rest as testing set.
    X_train = X_shuffled.copy()[:idx_train]
    X_test = X_shuffled.copy()[idx_train:]
    Y_train = Y_shuffled.copy()[:idx_train]
    Y_test = Y_shuffled.copy()[idx_train:]
    
    
    p_train, pcov_train = curve_fit(datafitting_transcription_experimental, X_train, Y_train, \
                        p0 = pexp, bounds = (0, 1e8)) # Do the parameter estimation on training set
    
    param_matrix[i] = p_train
    y_sim_train = datafitting_transcription_experimental(X_train, *p_train)
    y_sim_test = datafitting_transcription_experimental(X_test, *p_train)
    
    
    print('Fold ', i+1, ': ')
    print('Parameters: ', p_train)
    print('Training MAE: ', np.average(np.abs(y_sim_train - Y_train)))
    print('Testing MAE: ', np.average(np.abs(y_sim_test - Y_test)))
    
    col = [] # Each original sample is allocated blue if it was in the training, and red if it was in testing set
    for i in range(len(X_copy)):
        if list(X_copy[i]) in X_train.tolist():
            col.append('b')
        else:
            col.append('r')
    
    fig = plt.figure(figsize = (14, 6))
    pl.rc('axes', linewidth=2) # setting the width of the figure border lines thicker (2)
    ax = fig.add_subplot(121)
    
    ax = plot_experimental_data(ax, C_RNA_new, C_RNA_stdev) # Set experimental data as background
    
    ## Create prediction plots as in previous cell
    y_sim = datafitting_transcription_experimental(X_copy, *p_train)
    ax.scatter(np.arange(len(y_sim)), y_sim, c = col, label = 'Simulation')
    ax.set_xlabel('Index of experimental dataset')
    ax.set_ylabel('Effective RNA yield CQA [g/L]')
    
    ax2 = fig.add_subplot(122)
    ax2.plot([0, max(max(C_RNA_new), max(y_sim))], [0, max(max(C_RNA_new), max(y_sim))], '-k', label = 'Identity line')
    ax2.scatter(C_RNA_new, y_sim, c = col)
    ax2.set_ylabel('Predicted RNA yield CQA [g/L]')
    ax2.set_xlabel('Actual RNA yield CQA [g/L]')
    
    plt.show(); # Show plot
    plt.clf(); # Clear plot
    
    X_shuffled = np.concatenate((X_test, X_train), axis = 0) # Queue testing set to the front of the array and repeat
    