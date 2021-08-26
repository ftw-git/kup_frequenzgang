# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 12:43:00 2019

@author: FWass
"""
import numpy as np
import os

from functions.fg_norm import fg_norm

def write_data_lex(blade_name, data):
    
    data_keys = []
    check = False
    counter = 0
    eigen = {}
    
    for key in data.keys():
        
        data_keys.append(key)
        data_formated = np.empty(np.shape(data[key])).T
        
        for column in range(0, len(data[key])):

            data_formated[:, column] = data[key][column]
            eigen[key] = fg_norm(data[key][0], data[key][2])[0]
    
    filepath = os.getcwd()
    
    # Loop over cwd to check if already existing
    for (dirpath, dirs, files) in os.walk(filepath + '/data_lex'): #Schleife über Ordner und eventuelle Unterverzeichnisse
        
        for file in files: #Schleife über Messdateien

            if blade_name in file:
                check = True
                counter += 1 # Counts preexisting files of same name
                                
    if not check:
        with open('data_lex/' + blade_name, 'a+') as data_lex_file:    
            np.savetxt(data_lex_file, data_formated, fmt = '%s', delimiter = ' ', header = 'Frequency Voltage Deflection')
            
    else:        
         with open('data_lex/' + blade_name + '_{}'.format(counter), 'a+') as data_lex_file:    
            np.savetxt(data_lex_file, data_formated, fmt = '%s', delimiter = ' ', header = 'Frequency Voltage Deflection')
      
    # Close file to be sure
    data_lex_file.close()
    
    print(eigen)

    return data_formated