# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:13:10 2019

@author: FWass
"""

def comparison_ask():
    
    test = 'false' #Initialisierung Test-Variable

    while (test != 'true'): #Abfrage ob Vergleich oder einfache Auswertung
          
        cond_comp = input('Compare? [y/n]\n') #Abfrage
        
        if (cond_comp != 'y') & (cond_comp != 'n'): #Fehleingabe abfangen
            
            print('Warning: y or n expected!')
        
        else:
            
            test = 'true' 
    
    if cond_comp == 'y': #Vergleich
        
        counter_lim = 2
        
    else: counter_lim = 1 #Einfache Auswertung
    
    return (counter_lim, test)