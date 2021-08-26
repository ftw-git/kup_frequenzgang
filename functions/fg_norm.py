# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:08:32 2019

@author: FWass
"""

import numpy as np

def fg_norm(frequency, deflection):
    
    EF = np.array(np.flip([np.argmax(deflection), frequency[np.argmax(deflection)]]))
    
    frequency_norm = frequency/EF[0]
    deflection_norm = deflection/EF[1]
    
    return (EF, frequency_norm, deflection_norm)
        
        