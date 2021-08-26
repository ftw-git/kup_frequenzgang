# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 13:24:54 2019

@author: FWass
"""

import numpy as np

def data_slicing(frequency, tension):
    
    #Initialisierung der Anregungsspannung-Slice-Grenzen
    slicing_lim_tdms = np.zeros((int(np.max(tension)/50)+1,1))
    
    #Durchläuft Anregungsspannung und liest Zeilen mit Spannungssprüngen ein
    for spg in range(int(np.min(tension)), int(np.max(tension))+50, 50):
        
            for zeile in range(0, np.size(frequency, 0)):
                
                if tension[zeile] == spg:
    
                    slicing_lim_tdms[int(spg/50)] = zeile
    
    return slicing_lim_tdms