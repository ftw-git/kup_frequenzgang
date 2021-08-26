# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 14:16:24 2019

@author: FWass
"""
import os
import tkinter as tk #GUI, ermöglicht Fenster Dialog

from tkinter import filedialog
from nptdms import TdmsFile

def data_read_tdms(DATA_FILE_tdms = 'nan'):
    
    ## Importieren der TDMS-Datei        
    
    test = 'false' #Initialisierung Testvariable
    
    if DATA_FILE_tdms == 'nan':
    
        while (test!='true'): #Fehler abfangen
        
            root = tk.Tk()
            root.withdraw()
            DATA_FILE_tdms = filedialog.askopenfilename(title='Choose .tdms file to be read.') #Öffnet Dialogfenster zum auswählen der Datei
            
            #DATA_FILE_tdms = "C:/Users/FWass/Desktop/TUB_HiWi/Schaufelmessung/FG_Rot_Mit_Deckel/Frequenzgang_NEU_Schaufel_senkrecht_02_08_2019_Messkampagne_7.tdms" #Alternativ: absoluten Pfad angeben
        
            if os.path.splitext(DATA_FILE_tdms)[1]=='.tdms': #Falls .tdms-Datei
                
                test = 'true'
                
            else:
                
                print('Warning: .tdms-Data expected!') #Falls keine .tdms-Datei
        
        
    tdms_file = TdmsFile(DATA_FILE_tdms) #tdms-file einlesen
    
    channels = tdms_file.group_channels(tdms_file.groups()[0]) #Kanäle der tdms werden der Var channels zugewiesen
        
    
    return (DATA_FILE_tdms, channels)