# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:17:52 2019

@author: FWass
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

from functions.blade_evaluation import blade_evaluation
from functions.fg_norm import fg_norm
from functions.comparison_ask import comparison_ask
from functions.common_delta_volt import common_delta_volt
from functions.write_data_lex import write_data_lex

plt.close('all') #Zum schließen eventuell offener Plots bei erneutem kompilieren

#%% Auswertung - Datenvergleich

(counter_lim, cond_comp) = comparison_ask() #Abfrage - Vergleich oder einfache Auswertung

data_to_compare = {} #Initialisierung der zu vergleichenden Daten in einem Dictionnary
data_key = {}
filename = {}

for counter in range(0, counter_lim): #Schleife über die zu vergleichenden Schaufeln
    
    data_key[counter] = "%d. Schaufel" %(counter+1) #Key = Nummer der Schaufel
    
    (data_to_compare[data_key[counter]], filename[data_key[counter]]) = blade_evaluation(cond_comp) #Holt Frequenzgang aus Auswerteskript

if cond_comp == 'y':

    keys_in_common = common_delta_volt(data_to_compare, data_key)

#%% Relevante Plots

for zaehler in range(0, len(data_key)-1):
    
    for voltage in keys_in_common:
        
        plt.figure()
        plt.plot(data_to_compare[data_key[zaehler]][voltage][0], data_to_compare[data_key[zaehler]][voltage][2], 'g', label = r'1.Schaufel')
        plt.plot(data_to_compare[data_key[zaehler+1]][voltage][0], data_to_compare[data_key[zaehler+1]][voltage][2], 'r', label = r'2.Schaufel')
        plt.title(voltage + '_Frequenzgang')
        plt.xlabel('Frequenz [Hz]')
        plt.ylabel('Auslenkung  [mm]')
        plt.legend(loc='best')
        plt.grid()
        plt.tight_layout()
        plt.show()

#Normieren der Frequenzgänge

if cond_comp == 'y':   
     
    data_to_compare_norm = {}
        
    data_to_compare_norm = {str(key): {str(voltage): fg_norm(data_to_compare[key][voltage][0], data_to_compare[key][voltage][2]) for voltage in keys_in_common} for key in data_to_compare.keys()}
    #Normiert in Schleife über alle Spannungszustände und alle Schaufeln den Frequenzgang auf die Eigenfrequenz
    
    #Plot der normierten Frequenzgänge beider Schaufeln im Vergleich für jeden Spannungszustand
    for key in range(0, len(data_key)-1):
        
        for voltage in keys_in_common:
    
            plt.figure()
            plt.plot(data_to_compare_norm[data_key[key]][voltage][1], data_to_compare_norm[data_key[key]][voltage][2], 'g', label =  r'%s - normiert' %(data_key[key]))
            plt.plot(data_to_compare_norm[data_key[key+1]][voltage][1], data_to_compare_norm[data_key[key+1]][voltage][2], 'r', label =  r'%s - normiert' %(data_key[key+1]))
            plt.title(str(voltage))
            plt.xlabel('Frequenz [Hz]')
            plt.ylabel('Auslenkung  [mm]')
            plt.legend(loc='best')
            plt.grid()
            plt.tight_layout()
            plt.show()


#%% Relevante Daten zu einzelnen Schaufeln in eigene Datei schreiben

# =============================================================================
# for key in filename.keys():
# 
#     data_formated = write_data_lex(os.path.splitext(filename[key])[0], data_to_compare[key])
# =============================================================================
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
