# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 21:06:03 2020

@author: dv516
"""

import numpy as np

from scipy.optimize import fsolve

from mass_balances import solution_proton, solution

import warnings


def initguess(C_Mg_0, C_NTP_0, C_PPi_0, C_HEPES_0):
## Return an initial guess that gives a physical solution to solution()  
    
    # Initial guess, hopefully in the right order of magnitude
    guess = [C_Mg_0, C_NTP_0, C_HEPES_0, C_PPi_0, C_HEPES_0]

    guess_array = np.logspace(0, 2, 151) # Array of constants to multiply guess with

    tot_conc0 = (C_Mg_0, C_NTP_0, 10**(-7.5), C_PPi_0, C_HEPES_0) # Parameters needed for solution() 
    
    # Check if dividing or multiplying guess by a constant gives a feasible solution and if so, return it
    for g in guess_array:
        gu1 = guess/g
        all_vals = fsolve(solution, gu1, args = tot_conc0)
        if (all_vals>0).all():
            return gu1
        gu2 = np.array(guess)*g
        all_vals = fsolve(solution, gu2, args = tot_conc0)
        if (all_vals>0).all():
            return gu2
    
    # Warning if no solution was found
    warnings.warn("Warning: No valid initial guess found at initguess() for: ",
                  C_Mg_0, C_NTP_0, C_PPi_0, C_HEPES_0)
    return guess

def guess_tr(tot_conc, guess):
## Look for a multiple of guess, the free solution concentration of the previous time step, 
## that gives a physical solution

    # Array of multiples
    guess_array = np.logspace(0, 4, 151)
    
    # Check if dividing or multiplying guess by a constant gives a feasible solution, and if so, return guess
    for g in guess_array:
        gu1 = np.array(guess)*g
        all_vals = fsolve(solution_proton, gu1, args = tot_conc)
        if (all_vals>0).all():
            return gu1
        gu2 = np.array(guess)/g
        all_vals = fsolve(solution_proton, gu2, args = tot_conc)
        if (all_vals>0).all():
            return gu2
    
    # If none was found, repeat for guesses at different orders of magnitude
    guess_array = np.logspace(2, 8, 151)
    for g in guess_array:
        gu1 = np.array(guess)*g
        all_vals = fsolve(solution_proton, gu1, args = tot_conc)
        if (all_vals>0).all():
            return gu1
        gu2 = np.array(guess)/g
        all_vals = fsolve(solution_proton, gu2, args = tot_conc)
        if (all_vals>0).all():
            return gu2
    
    # Warning if no solution was found (It shouldn't come to this)
    warnings.warn("Warning: No physical solution to free solution concentration computation found at: "                   + str(tot_conc))
    
    return guess

def all_sol_conc(C_Mg, C_NTP, C_H, C_PPi):
## Give concentrations of complex components from free solution components given by association equilibria 

  C_HNTP = (C_H*C_NTP) * 10**(6.95)
  C_MgNTP = (C_Mg*C_NTP) * 10**(4.42)
  C_Mg2NTP = (C_Mg*C_MgNTP) * 10**(1.69)
  C_MgHNTP = (C_Mg*C_HNTP) * 10**(1.49)
  C_MgPPi = (C_Mg*C_PPi) * 10**(5.42)
  C_Mg2PPi = (C_Mg*C_MgPPi) * 10**(2.33)
  C_HPPi = (C_H*C_PPi) * 10**(8.94)
  C_H2PPi = (C_H*C_HPPi) * 10**(6.13)
  C_MgHPPi = (C_Mg*C_HPPi) * 10**(3.05)

  return C_HNTP, C_MgNTP, C_Mg2NTP, C_MgHNTP, C_MgPPi, C_Mg2PPi, C_HPPi, C_H2PPi, C_MgHPPi


def tr(t, y, trans_constants, guess):
## tr() takes in the timestep t (not really used), the vector of concentrations y at time step t, constants
## to be used in differential expressions, and the initial guess to the algebraic solver guess
## and returns the differential expressions, the changes in the 6 component concentrations at that time step

    # Retrieve concentrations at current timestep
    C_RNA = y[0]
    C_PPi_0 = y[1]
    C_NTP_0 = y[2]
    C_H_0 = y[3]
    C_T7RNAP = y[4]
    C_Mg_0 = y[5]
    C_HEPES_0 = 0.04
    
    # Retrieve constants to be used in kinetic expressions
    k_app, K1, K2, alpha, Nall, k_d, k_ac, n_ac, k_ba, n_ba, k_Mg, n_Mg, n_RNA, K3, K4, K5, k_precip = trans_constants

    # total concentrations to be used to calculate free solution concentrations
    tot_conc = (C_Mg_0, C_NTP_0, C_H_0, C_PPi_0, C_HEPES_0)
    
    # get free solution concentrations associated with first guess
    free_sol_conc = fsolve(solution_proton, guess, args = tot_conc)
    
    # Only update guess if the solution is non-negative, otherwise call function to find valid guess
    if (free_sol_conc >= 0).all():
        guess = free_sol_conc # Only update guess if the solution makes physical sense
    else:
        guess = guess_tr(tot_conc, guess)
        free_sol_conc = fsolve(solution_proton, guess, args = tot_conc) 
    
    C_Mg, C_NTP, C_H, C_PPi, C_HEPES = free_sol_conc

    complex_conc = all_sol_conc(C_Mg, C_NTP, C_H, C_PPi) # get complex component concentrations
    C_HNTP, C_MgNTP, C_Mg2NTP, C_MgHNTP, C_MgPPi, C_Mg2PPi, C_HPPi, C_H2PPi, C_MgHPPi = complex_conc
    
    V_tr = k_app*C_T7RNAP*alpha*C_Mg*C_MgNTP            / (K5*C_MgNTP**2*C_Mg +  K4*C_Mg**2 + K3*C_MgNTP**2                + K2*C_MgNTP + K1*alpha*C_Mg + 1) # Rate of transcription
    # K3 to K5 are 'investigative' terms, whose influence was tested during previous parameter estimations
    # They are however not used in simulations

    V_deg = (k_ac*C_H**n_ac + k_ba*(10**(-14)/(C_H))**n_ba + k_Mg*C_Mg**n_Mg)             * C_RNA**n_RNA # Rate of degradation
    
    # Rate of precipitation. Only for investigative purposes. In simulations, k_precip = 0
    C_Mg2PPi_equ = 0.014e-3
    V_precip = k_precip*max(0 , C_Mg2PPi - C_Mg2PPi_equ) # = 0 for the following simulations
    
    # Differential expressions
    df0dt = V_tr - V_deg  # Change in total RNA
    df1dt = (Nall-1)*V_tr - 1 * V_precip # Change in total PPi
    df2dt = - Nall*V_tr # Change in total NTP
    df3dt = (Nall-1)*V_tr # Change in total H
    df4dt = - k_d*C_T7RNAP # Change in T7RNAP ## k_d = 0
    df5dt = - 2 * V_precip # Change in total Mg
    
    # Return both differential expressions and updated guess
    return np.array([df0dt, df1dt, df2dt, df3dt, df4dt, df5dt]), guess


def datafitting_transcription_experimental(X, k_app, K1, K2, k_ac, k_ba, k_Mg, K3, K4, K5):
## Takes the N x D vector X, where each of the N rows denotes one set of D-dimensional inputs
## and various parameters to return the array of RNA concentrations associated with each of the N samples
# Columns of X: 1st column: end time of simulation, 2nd: C_Mg_0 initial, 3rd: C_NTP_0 initial
# 4th: C_spermidine (not used after all), 5th: C_T7RNAP

# Note that this function is not efficient for calculating different RNA yields at different times during the same
# trajectory (at same Mg, NTP and spermidine). Even for the same initial concentrations every time step along the 
# same trajectory needs to be called separately
    
    t_array = X[:,0] # size N: 1st column of X vector, the time of each of the N samples
    N = len(t_array)
    
    # The following parameters were set to avoid overfitting or through a priori knowledge
    k_d = 0
    alpha = 1
    Nall = 10000
    n_ac = 1
    n_ba = 1
    n_Mg = 1
    n_RNA = 1
    
    ## Mainly used for parameter estimation. 
    ## Uncomment one of these to fix a certain parameter no matter what its input
    #K1 = 0
    #K2 = 0
    #k_ac = 0
    #k_ba = 0
    K3 = 0
    K4 = 0
    K5 = 0
    k_precip = 0
    
    C_array = np.zeros(N)
    
    for j in range(N):

        C_PPi_0 = 1e-18 # [M], Always assme non-zero, very small initial C_PPi_0 for numerical stability
        C_H = 10**(-7.5) # [M], Assume initial pH is 7.5, at pKa of HEPES buffer
        C_Mg_0 = X[j,1] # [M]
        C_NTP_0 = X[j,2] # [M]
        C_T7RNAP_0 = X[j,4] # [U/L]
        C_HEPES_0 = 0.04 # Always assume same concentration for initial HEPES
        
        ## At first, find C_H_0 given C_H
        ## Find valid initial guess, especially important at first time step
        tot_conc0 = (C_Mg_0, C_NTP_0, C_H, C_PPi_0, C_HEPES_0) 
        guess = initguess(C_Mg_0, C_NTP_0, C_PPi_0, C_HEPES_0)
        all_vals = fsolve(solution, guess, args = tot_conc0)
        if (all_vals<0).any():
            warnings.warn("Warning: No valid guess found for solution() at initial conditions of: " + str(tot_conc0))
        
        # After knowing C_H_0, all total concentrations are known to find guess free solution concentrations
        C_H_0 = all_vals[2]
        guess[2] = C_H
        tot_conc0 = (C_Mg_0, C_NTP_0, C_H_0, C_PPi_0, C_HEPES_0)
        guess = fsolve(solution_proton, guess, args = tot_conc0) # Set guess for future timesteps
        
        # Retrieve parameters used for kinetics
        trans_cons = (k_app, K1, K2, alpha, Nall, k_d, k_ac, n_ac, k_ba, n_ba,                       k_Mg, n_Mg, n_RNA, K3, K4, K5, k_precip)
        
        N_step = 101 # 101 time steps for numerical integration was found to give OK stability
        U_0 = [0, C_PPi_0, C_NTP_0, C_H_0, C_T7RNAP_0, C_Mg_0] # Set array of initial concentrations

        t1 = np.linspace(0, t_array[j], N_step) # Set time grid for current simulation
        delta_t = t1[1]-t1[0] # Time step
        U1 = np.zeros((N_step, len(U_0))) # Matrix storing concentrations at all intermediate time steps for given D-dimensional input
        U1[0] = U_0 # Initial conditions
        
        ## Runge-Kutta 4 solver to iterate over each time step
        for n in range(N_step-1):
            y_1 = U1[n]
            dy_1, guess = tr(t1[n], y_1, trans_cons, guess)
            y_2 = U1[n] + 0.5 * delta_t * dy_1
            dy_2, guess = tr(t1[n] + 0.5 * delta_t, y_2, trans_cons, guess)
            y_3 = U1[n] + 0.5 * delta_t * dy_2
            dy_3, guess = tr(t1[n] + 0.5 * delta_t, y_3, trans_cons, guess)
            y_4 = U1[n] + delta_t * dy_3
            dy_4, guess = tr(t1[n] + delta_t, y_4, trans_cons, guess)
            U1[n+1] = U1[n] + delta_t / 6.0 * (dy_1+ 2.0 * dy_2 + 2.0 * dy_3 + dy_4) 
        
        if (U1[-1]<0).any():
            warnings.warn("Warning: Negative final concentration for the following initial concentrations: " + str(tot_conc0))
            C_array[j] = -0.001
        elif np.isfinite(U1[-1]).all():
            C_array[j] = U1[-1][0]*5e6 # RNA molecular weight
        else:
            warnings.warn("Warning: Unstable solution for the following initial concentrations: " + str(tot_conc0))
            C_array[j] = 1e4
    # Return N-dimensional RNA yield corresponding to X    
    return C_array


def datafitting_transcription_prob_stDev(X, stDev_per, k_app, K1, K2, k_ac, k_ba, k_Mg, K3, K4, K5):
## Same function as datafitting_transcription_experimental() but with the additional input of a standartd deviation
## and consequent random sampling of parameters
    
    # Normal distribution around k_app with stDev_per standard deviation
    k_app = np.random.normal(k_app, stDev_per[0])
    K1 = np.random.normal(K1, stDev_per[1])
    K2 = np.random.normal(K1, stDev_per[2])
    k_ac = np.random.normal(k_ac, stDev_per[3])
    k_ba = np.random.normal(k_ba, stDev_per[4])
    k_Mg = np.random.normal(k_Mg, stDev_per[5])
    K3 = np.random.normal(K3, stDev_per[6])
    K4 = np.random.normal(K4, stDev_per[7])
    K5 = np.random.normal(K5, stDev_per[8])
    
    k_d = 0
    alpha = 1
    Nall = 10000
    n_ac = 1
    n_ba = 1
    n_Mg = 1
    n_RNA = 1
    k_precip = 0
    
    t_array = X[:,0]
    N = len(t_array)
    C_array = np.zeros(N)
    
    for j in range(N):

        C_PPi_0 = 1e-18 # [M], Always assme non-zero, very small initial C_PPi_0 for numerical stability
        C_H = 10**(-7.5) # [M], Assume initial pH is 7.5, at pKa of HEPES buffer
        C_Mg_0 = X[j,1] # [M]
        C_NTP_0 = X[j,2] # [M]
        C_T7RNAP_0 = X[j,4] # [U/L]
        C_HEPES_0 = 0.04 # Always assume same concentration for initial HEPES
        
        ## At first, find C_H_0 given C_H
        ## Find valid initial guess, especially important at first time step
        tot_conc0 = (C_Mg_0, C_NTP_0, C_H, C_PPi_0, C_HEPES_0) 
        guess = initguess(C_Mg_0, C_NTP_0, C_PPi_0, C_HEPES_0)
        all_vals = fsolve(solution, guess, args = tot_conc0)
        if (all_vals<0).any():
            warnings.warn("Warning: No valid guess found for solution() at initial conditions of: "                   + str(tot_conc0))
        
        # After knowing C_H_0, all total concentrations are known to find guess free solution concentrations
        C_H_0 = all_vals[2]
        guess[2] = C_H
        tot_conc0 = (C_Mg_0, C_NTP_0, C_H_0, C_PPi_0, C_HEPES_0)
        guess = fsolve(solution_proton, guess, args = tot_conc0) # Set guess for future timesteps
        
        # Retrieve parameters used for kinetics
        trans_cons = (k_app, K1, K2, alpha, Nall, k_d, k_ac, n_ac, k_ba, n_ba,                       k_Mg, n_Mg, n_RNA, K3, K4, K5, k_precip)
        
        N_step = 101 # 101 time steps for numerical integration was found to give OK stability
        U_0 = [0, C_PPi_0, C_NTP_0, C_H_0, C_T7RNAP_0, C_Mg_0] # Set array of initial concentrations

        t1 = np.linspace(0, t_array[j], N_step) # Set time grid for current simulation
        delta_t = t1[1]-t1[0] # Time step
        U1 = np.zeros((N_step, len(U_0))) # Matrix storing concentrations at all intermediate time steps for given D-dimensional input
        U1[0] = U_0 # Initial conditions
        
        ## Runge-Kutta 4 solver to iterate over each time step
        for n in range(N_step-1):
            y_1 = U1[n]
            dy_1, guess = tr(t1[n], y_1, trans_cons, guess)
            y_2 = U1[n] + 0.5 * delta_t * dy_1
            dy_2, guess = tr(t1[n] + 0.5 * delta_t, y_2, trans_cons, guess)
            y_3 = U1[n] + 0.5 * delta_t * dy_2
            dy_3, guess = tr(t1[n] + 0.5 * delta_t, y_3, trans_cons, guess)
            y_4 = U1[n] + delta_t * dy_3
            dy_4, guess = tr(t1[n] + delta_t, y_4, trans_cons, guess)
            U1[n+1] = U1[n] + delta_t / 6.0 * (dy_1+ 2.0 * dy_2 + 2.0 * dy_3 + dy_4) 
        
        if (U1[-1]<0).any():
            warnings.warn("Warning: Negative final concentration for the following initial concentrations: " + str(tot_conc0))
            C_array[j] = -0.001
        elif np.isfinite(U1[-1]).all():
            C_array[j] = U1[-1][0]*5e6 
        else:
            warnings.warn("Warning: Unstable solution for the following initial concentrations: " + str(tot_conc0))
            C_array[j] = 1e4
    # Return N-dimensional RNA yield corresponding to X    
    return C_array