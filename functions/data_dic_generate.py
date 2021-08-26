# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 15:00:16 2019

@author: FWass
"""

import numpy as np

def data_dic_generate(channels, slicing_lim_tdms):
    
    key = {}
    data_dic = {} #Initialisieren des dic und der keys
    
    for data_dic_key in range(int(np.min(channels[1].data)/50), int(np.max(channels[1].data)/50)+1): #Schleife über alle Spannungszustände
        
        key[data_dic_key]="%dV" %((data_dic_key)*50) #Erstellt key
        
        #Erstellt dictionnary mit allen Spannungszuständen und den dazugehörigen relevanten channels
        data_dic[key[data_dic_key]] = [channels[0].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
                  channels[1].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
                  channels[3].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1]]
            
    return data_dic