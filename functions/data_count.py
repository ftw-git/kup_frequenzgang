# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 12:56:24 2019

@author: FWass
"""
import os

def data_count(filepath):

    #%% Zählt die zu der Messkampagne dazugehörigen Dateien 
    num_of_files = 0
    num_of_measure = [] #Initialisieren Anzahl der Dateien und Nummerierung der Messung
    
    for (dirpath, dirs, files) in os.walk(filepath): #Schleife über Ordner und eventuelle Unterverzeichnisse
        
        print(files)
        
        for file in files: #Schleife über Messdateien
            
            print(file)
            print(os.path.splitext(file)[0][:-1])
            
            if (os.path.splitext(file)[0][:-1] == 'Msg') & (os.path.splitext(file)[1] == '.dat'): #Wenn Root_Name und Dateityp übereinstimmen
            
                num_of_files = num_of_files + 1 #Zählt Anzahl der zur Messkampagne dazugehörigen Dateien
                num_of_measure.append(int(os.path.splitext(file)[0][-1])) #Generiert die Nummerierung der selbigen Datei aus Dateinamen
                    
    return (num_of_files, num_of_measure)