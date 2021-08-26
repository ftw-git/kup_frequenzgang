# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 15:38:41 2019

@author: FWass
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import pandas as pd #Einlesen von EXCEL-Dateien
import xlrd #Einlesen von EXCEL-Dateien
import tkinter as tk #GUI, ermöglicht Fenster Dialog

from tkinter import filedialog
from nptdms import TdmsFile

def compare_FG():
    
    plt.close('all') #Zum schließen eventuell offener Plots bei erneutem kompilieren
    
    ## Importieren der TDMS-Datei        
    
    test = 'false' #Initialisierung Testvariable
    
    while (test!='true'): #Fehler abfangen

        root = tk.Tk()
        root.withdraw()
        DATA_FILE_tdms_1 = filedialog.askopenfilename() #Öffnet Dialogfenster zum auswählen der Datei
        #WARUM? Hat vorher funktioniert?
        #DATA_FILE_tdms = "C:/Users/FWass/Desktop/TUB_HiWi/Schaufelmessung/FG_Rot_Mit_Deckel/Frequenzgang_NEU_Schaufel_senkrecht_02_08_2019_Messkampagne_7.tdms" #Alternativ: absoluten Pfad angeben

    if os.path.splitext(DATA_FILE_tdms_1)[1]=='.tdms': #Falls .tdms-Datei
        
        test = 'true'
        
    else:
        
        print('Warning: .tdms-Data expected!') #Falls keine .tdms-Datei

    
    tdms_file_1 = TdmsFile(DATA_FILE_tdms_1) #tdms-file einlesen

    channels_1 = tdms_file_1.group_channels(tdms_file_1.groups()[0]) #Kanäle der tdms werden der Var channels zugewiesen
    
    #Speichern der relevanten channels in dic
    data_dic_tdms = {"Anregungsfrequenz": channels_1[0].data, 
                     "Anregungsspannung": channels_1[1].data, 
                     "Laser-Mittlere Amplitude": channels_1[3].data
                     }
    
    #Initialisierung der Anregungsspannung-Slice-Grenzen
    slicing_lim_tdms_1 = np.zeros((int(np.max(channels_1[1].data)/50)+1,1))
    
    #Durchläuft Anregungsspannung und liest Zeilen mit Spannungssprüngen ein
    for spg in range(int(np.min(channels_1[1].data)), int(np.max(channels_1[1].data))+50, 50):
        
            for zeile in range(0, np.size(channels_1[0].data, 0)):
                
                if channels_1[1].data[zeile] == spg:
    
                    slicing_lim_tdms_1[int(spg/50)]=zeile
    
    key = {} #Initialisierung Spannung als dic-key
    
    #Schleife über alle Anregungsspannungen
    
    data= []
    
    for data_dic_key in range(int(np.min(channels_1[1].data)/50), int(np.max(channels_1[1].data)/50)+1):
        
        key[data_dic_key]="%dV" %((data_dic_key)*50) #Erstellt key
        
        #Erstellt dictionnary mit allen Spannungszuständen und den dazugehörigen relevanten channels
        
        data.append({key[data_dic_key]: (channels_1[0].data[int(slicing_lim_tdms_1[data_dic_key-1,0]):int(slicing_lim_tdms_1[data_dic_key,0])+1], 
                  channels_1[1].data[int(slicing_lim_tdms_1[data_dic_key-1,0]):int(slicing_lim_tdms_1[data_dic_key,0])+1], 
                  channels_1[3].data[int(slicing_lim_tdms_1[data_dic_key-1,0]):int(slicing_lim_tdms_1[data_dic_key,0])+1])})
    
    return 'true'