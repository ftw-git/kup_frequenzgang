# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:23:48 2019

@author: FWass
"""

#------------------------------------------------------------------------------
# import packages
#------------------------------------------------------------------------------

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import pandas as pd #Einlesen von EXCEL-Dateien
import xlrd #Einlesen von EXCEL-Dateien
import tkinter as tk #GUI, ermöglicht Fenster Dialog

from tkinter import filedialog
from nptdms import TdmsFile
from Functions.data_count import data_count
from Functions.data_validity_test import data_validity_test
from Functions.data_slicing import data_slicing
from Functions.data_read_tdms import data_read_tdms
from Functions.data_dic_generate import data_dic_generate

def blade_evaluation(cond_comp): 

    ## Importieren der TDMS-Datei        
        
    (DATA_FILE_tdms, channels) = data_read_tdms() #Liest die tdms-Datei ein und generiert channels als array
    
    (num_of_files, num_of_measure) = data_count(DATA_FILE_tdms) #Zählt die zu untersuchenden Messdateien
        
    #Speichern der relevanten channels in dic
    data_dic_tdms = {"Anregungsfrequenz": channels[0].data, 
                     "Anregungsspannung": channels[1].data, 
                     "Laser-Mittlere Amplitude": channels[3].data
                     }
    
    slicing_lim_tdms = data_slicing(channels[0].data, channels[1].data) #Ordnet die Spannungszuständen den Frequenz- und Auslenkungsverläufe zu
    
    key = {} #Initialisierung Spannung als dic-key
    
    data_dic = {}
    
    for data_dic_key in range(int(np.min(channels[1].data)/50), int(np.max(channels[1].data)/50)+1): #Schleife über alle Spannungszustände
        
        key[data_dic_key]="%dV" %((data_dic_key)*50) #Erstellt key
        
        #Erstellt dictionnary mit allen Spannungszuständen und den dazugehörigen relevanten channels
        
        data_dic[key[data_dic_key]] = [channels[0].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
                  channels[1].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
                  channels[3].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1]]
        
        if data_validity_test(data_dic[key[data_dic_key]][2]) == 'false': #Wenn fehlerhaftes Signal enthalten
            
            for zaehler in range(0, num_of_files): #Schleife über alle Messdateien mit selbem Root_Name
                
                file_new = os.path.splitext(DATA_FILE_tdms)[0][:-1] + '%d' %num_of_measure[zaehler] + '.tdms'
                
                (file_new, channels_new) = data_read_tdms(file_new) #Einlesen
                
                slicing_lim_tdms_new = data_slicing(channels_new[0].data, channels_new[1].data) #Spannungswechselindizes generieren
                
                data_dic_new = data_dic_generate(channels_new, slicing_lim_tdms_new) #Entsprechendes dic anlegen
                
                if key[data_dic_key] in data_dic_new: #Wenn Ersatzsignal zu entsprechendem fehlerhaften Signal vorhanden
                    
                    if data_validity_test(data_dic_new[key[data_dic_key]][2]) == 'true': #Wenn außerdem nicht fehlerhaft
                        
                        data_dic[key[data_dic_key]][0] = data_dic_new[key[data_dic_key]][0]
                        data_dic[key[data_dic_key]][2] = data_dic_new[key[data_dic_key]][2] #Ersetzen durch neues Signal
     
#%% Vervollständigen mit Spannungszuständen aus anderen Messdateien
    
    for zaehler in range(0, num_of_files):
        
        file_to_check = os.path.splitext(DATA_FILE_tdms)[0][:-1] + '%d' %num_of_measure[zaehler] + '.tdms'
        
        (file_to_check, channels_to_check) = data_read_tdms(file_to_check) #Einlesen
        
        slicing_lim_tdms_to_check = data_slicing(channels_to_check[0].data, channels_to_check[1].data) #Spannungswechselindizes generieren
        
        data_dic_to_check = data_dic_generate(channels_to_check, slicing_lim_tdms_to_check) #Entsprechendes dic anlegen
        
        for key in data_dic_to_check.keys():
            
            if (data_validity_test(data_dic_to_check[key][2]) == 'true') & (not key in data_dic.keys()):
                
                data_dic[key] = data_dic_to_check[key]
     
#%% Plot der Frequenzgänge
    
    for key in data_dic.keys():
    
        #Plot: Auslenkung über Frequenz
        plt.figure()
        plt.plot(data_dic[key][0], data_dic[key][2], 'g')
        plt.title(key + '_Frequenzgang')
        plt.xlabel('Frequenz [Hz]')
        plt.ylabel('Auslenkung  [mm]')
        plt.grid()
        plt.tight_layout()
        plt.show()
    
    return (data_dic, os.path.split(DATA_FILE_tdms)[-1])
