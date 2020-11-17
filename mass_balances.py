# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 20:54:59 2020

@author: dv516
"""


def solution(variables, *constants):
## Mass balances for free solution concentrations
## C_H is given, and need to find C_H_0


    C_Mg_0, C_NTP_0, C_H, C_PPi_0, C_HEPES_0 = constants # Retrieve total concentrations
    
    C_Mg, C_NTP, C_H_0, C_PPi, C_HEPES = variables # Free solution concentrations to be solved
    
    ## E1-E10 substituted into M1-M5
    
    # Mg mass balance
    mass1 = C_Mg + (C_Mg*C_PPi)*10**(5.42) + \
                   (C_Mg*(C_H*C_PPi)*10**(8.94))*10**(3.05) + \
                   (C_Mg*C_NTP)*10**(4.42) + \
                 2*(C_Mg*(C_Mg*C_NTP)*10**(4.42))*10**(1.69) + \
                 2*(C_Mg*(C_Mg*C_PPi)*10**(5.42))*10**(2.33) + \
                   (C_Mg*(C_H*C_NTP) *10**(6.95))*10**(1.49) - C_Mg_0
    
    # NTP mass balance
    mass2 = C_NTP + (C_H*C_NTP) *10**(6.95) + \
                    (C_Mg*C_NTP)*10**(4.42) + \
                    (C_Mg*(C_Mg*C_NTP)*10**(4.42))*10**(1.69) + \
                    (C_Mg*(C_H*C_NTP) *10**(6.95))*10**(1.49) - C_NTP_0
    
    # H mass balance
    mass3 = C_H + (C_Mg*(C_H*C_NTP) *10**(6.95))*10**(1.49) + \
                  (C_H*C_NTP) *10**(6.95) + \
                  (C_Mg*(C_H*C_PPi)*10**(8.94))*10**(3.05) + \
                  (C_H*C_PPi)*10**(8.94) + \
                2*(C_H*(C_H*C_PPi)*10**(8.94))*10**(6.13) + \
                   (C_HEPES*C_H)*10**(7.5) - C_H_0
    
    # PPi mass balance
    mass4 = C_PPi + (C_Mg*C_PPi)*10**(5.42) + \
                    (C_Mg*(C_Mg*C_PPi)*10**(5.42))*10**(2.33) + \
                    (C_H*C_PPi)*10**(8.94) + \
                    (C_Mg*(C_H*C_PPi)*10**(8.94))*10**(3.05) + \
                    (C_H*(C_H*C_PPi)*10**(8.94))*10**(6.13) - C_PPi_0
    
    # HEPES mass balance
    mass5 = C_HEPES + C_HEPES*C_H*10**(7.5) - C_HEPES_0
    
    return (mass1, mass2, mass3, mass4, mass5)

def solution_proton(variables, *constants):
## Mass balances for free solution concentrations
## C_H_0 is given, and need to find C_H

    # Same as solution() but with total C_H_0 as input and C_H to be solved

    C_Mg_0, C_NTP_0, C_H_0, C_PPi_0, C_HEPES_0 = constants
    
    C_Mg, C_NTP, C_H, C_PPi, C_HEPES = variables
    
    # Mg mass balance
    mass1 = C_Mg + (C_Mg*C_PPi)*10**(5.42) + \
                   (C_Mg*(C_H*C_PPi)*10**(8.94))*10**(3.05) + \
                   (C_Mg*C_NTP)*10**(4.42) + \
                 2*(C_Mg*(C_Mg*C_NTP)*10**(4.42))*10**(1.69) + \
                 2*(C_Mg*(C_Mg*C_PPi)*10**(5.42))*10**(2.33) + \
                   (C_Mg*(C_H*C_NTP) *10**(6.95))*10**(1.49) - C_Mg_0
    
    # NTP mass balance
    mass2 = C_NTP + (C_H*C_NTP) *10**(6.95) + \
                    (C_Mg*C_NTP)*10**(4.42) + \
                    (C_Mg*(C_Mg*C_NTP)*10**(4.42))*10**(1.69) + \
                    (C_Mg*(C_H*C_NTP) *10**(6.95))*10**(1.49) - C_NTP_0
    
    # H mass balance
    mass3 = C_H + (C_Mg*(C_H*C_NTP) *10**(6.95))*10**(1.49) + \
                  (C_H*C_NTP) *10**(6.95) + \
                  (C_Mg*(C_H*C_PPi)*10**(8.94))*10**(3.05) + \
                  (C_H*C_PPi)*10**(8.94) + \
                2*(C_H*(C_H*C_PPi)*10**(8.94))*10**(6.13) + \
                   (C_HEPES*C_H)*10**(7.5) - C_H_0
    
    # PPi mass balance
    mass4 = C_PPi + (C_Mg*C_PPi)*10**(5.42) + \
                    (C_Mg*(C_Mg*C_PPi)*10**(5.42))*10**(2.33) + \
                    (C_H*C_PPi)*10**(8.94) + \
                    (C_Mg*(C_H*C_PPi)*10**(8.94))*10**(3.05) + \
                    (C_H*(C_H*C_PPi)*10**(8.94))*10**(6.13) - C_PPi_0
    
    # HEPES mass balance
    mass5 = C_HEPES + C_HEPES*C_H*10**(7.5) - C_HEPES_0
    
    return (mass1, mass2, mass3, mass4, mass5)
