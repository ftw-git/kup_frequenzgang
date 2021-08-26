# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 18:57:41 2019

@author: FWass
"""

def common_delta_volt(data, blade):
    
    keys_in_common = []
    
    for key in data[blade[0]]:
        
        if key in data[blade[1]].keys():

            keys_in_common.append(key)
    
    return keys_in_common