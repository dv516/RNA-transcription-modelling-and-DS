# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 18:04:47 2020

@author: dv516
"""

import numpy as np
import pandas as pd

def plot_experimental_data(ax, data, stdev):
## function that takes in an axis to plot the experimental RNA yield plus standard deviation and returns the axis
    
    r = np.arange(len(data))

    ax.errorbar(r[:11], data[:11], stdev[:11],                 label = 'Mg dependence after 2 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[11:11*2], data[11:11*2], stdev[11:11*2],                 label = 'Mg dependence after 4 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[11*2:11*3], data[11*2:11*3], stdev[11*2:11*3],                 label = 'Mg dependence after 6 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[33:33+4], data[33:33+4], stdev[33:33+4],                 label = 'T7RNAP dependence after 2 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[33+4:33+4*2], data[33+4:33+4*2], stdev[33+4:33+4*2],                 label = 'T7RNAP dependence after 4 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[33+4*2:33+4*3], data[33+4*2:33+4*3], stdev[33+4*2:33+4*3],                 label = 'T7RNAP dependence after 6 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[-6:-3], data[-6:-3], stdev[-6:-3],                 label = 'NTP dependence at low Mg after 2 hrs', linestyle = '', marker = 'o', markersize = 3)
    ax.errorbar(r[-3:], data[-3:], stdev[-3:],                  label = 'NTP dependence at low Mg after 2 hrs', linestyle = '', marker = 'o', markersize = 3)
    
    return ax

def get_experimental_data():
    ## Extract the data from specific columns in a certain sheet in a given Excel file as dataframe
    data_DoE = pd.read_excel('Experimental_Data.xlsx', 'ModelData')
    dataFrame_DoE = pd.DataFrame(data_DoE, columns = ['t [hr]', 'Mg2+ [M]', 'NTP [M]', 'Spermidine [M]',                                                   'T7RNAP ratio', 'template ratio', 'RNA [g/L]'])
    # Last column of dataframe is RNA yield, while all the other columns are 'input'
    indVar_DoE = dataFrame_DoE.values
    C_RNA_exp2 = indVar_DoE[:,-1]
    X2 = indVar_DoE[:,:-1]

    ## Each data sample appears three times, so take the average, st dev, 
    ## mean absolute error and mean square error of data set
    C_RNA_new = np.zeros(int(len(C_RNA_exp2)/3))
    X_new = np.zeros((int(len(X2)/3), len(X2[0])))
    C_RNA_stdev = np.zeros(int(len(C_RNA_exp2)/3))
    MAE_exp = []
    MSE_exp = []
    for i in range(len(C_RNA_exp2)):
        if i%3 == 0:
            C_RNA_new[int(i/3)] = np.average(C_RNA_exp2[i:i+3])
            X_new[int(i/3)] = X2[i]
            C_RNA_stdev[int(i/3)] = np.std(C_RNA_exp2[i:i+3])
        
            MAE_exp.append(np.sum(abs(C_RNA_exp2[i:i+3]-np.average(C_RNA_exp2[i:i+3]))/3))
            MSE_exp.append(np.sum((C_RNA_exp2[i:i+3]-np.average(C_RNA_exp2[i:i+3]))**2/3))
    return X_new, C_RNA_new, C_RNA_stdev, MAE_exp, MSE_exp