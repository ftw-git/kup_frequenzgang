# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:32:35 2019

@author: FWass
"""

import statistics

def data_validity_test(crit_data):
    
    crit_var = statistics.stdev(crit_data)
    
    if crit_var < 10**(-2):
        
        test = 'false'
    
    else: test = 'true'   
    
    return test