# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 12:43:00 2019

@author: FWass
"""
import numpy as np

def write_data_lex(blade_name, data):
    
    #data_lex_file = open(blade_name, 'a+')
    
    data_keys = []
    
    for key in data.keys():
        
        data_keys.append(key)
    
    for key in data.keys():
    
        data_formated = np.empty(np.shape(data[key])).T
        
        for column in range(0, len(data[key]-1)):
            
            print(column)
            data_formated[:, column] = data[key][column]
    

    with open(blade_name, 'a+') as data_lex_file:

        np.savetxt(data_lex_file, data_formated, fmt = '%s', delimiter = ' ', header = 'Frequency Voltage Deflection')
    
    data_lex_file.close()
    

    return data_formated