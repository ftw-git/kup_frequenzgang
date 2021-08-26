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

from Functions.compare_FG import compare_FG

##

plt.close('all') #Zum schließen eventuell offener Plots bei erneutem kompilieren

## Importieren der TDMS-Datei        

test = 'false' #Initialisierung Testvariable


while (test!='true'): #Fehler abfangen

    root = tk.Tk()
    root.withdraw()
    DATA_FILE_tdms = filedialog.askopenfilename() #Öffnet Dialogfenster zum auswählen der Datei
    #WARUM? Hat vorher funktioniert?
    #DATA_FILE_tdms = "C:/Users/FWass/Desktop/TUB_HiWi/Schaufelmessung/FG_Rot_Mit_Deckel/Frequenzgang_NEU_Schaufel_senkrecht_02_08_2019_Messkampagne_7.tdms" #Alternativ: absoluten Pfad angeben

    if os.path.splitext(DATA_FILE_tdms)[1]=='.tdms': #Falls .tdms-Datei
        
        test = 'true'
        
    else:
        
        print('Warning: .tdms-Data expected!') #Falls keine .tdms-Datei

tdms_file = TdmsFile(DATA_FILE_tdms) #tdms-file einlesen

channels = tdms_file.group_channels(tdms_file.groups()[0]) #Kanäle der tdms werden der Var channels zugewiesen

#Speichern der relevanten channels in dic
data_dic_tdms = {"Anregungsfrequenz": channels[0].data, 
                 "Anregungsspannung": channels[1].data, 
                 "Laser-Mittlere Amplitude": channels[3].data
                 }

#Initialisierung der Anregungsspannung-Slice-Grenzen
slicing_lim_tdms = np.zeros((int(np.max(channels[1].data)/50)+1,1))

#Durchläuft Anregungsspannung und liest Zeilen mit Spannungssprüngen ein
for spg in range(int(np.min(channels[1].data)), int(np.max(channels[1].data))+50, 50):
    
        for zeile in range(0, np.size(channels[0].data, 0)):
            
            if channels[1].data[zeile] == spg:

                slicing_lim_tdms[int(spg/50)]=zeile

key = {} #Initialisierung Spannung als dic-key

#Schleife über alle Anregungsspannungen

data= []

for data_dic_key in range(int(np.min(channels[1].data)/50), int(np.max(channels[1].data)/50)+1):
    
    key[data_dic_key]="%dV" %((data_dic_key)*50) #Erstellt key
    
    #Erstellt dictionnary mit allen Spannungszuständen und den dazugehörigen relevanten channels
    data_dic = {key[data_dic_key]: (channels[0].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
              channels[1].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
              channels[3].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1])}
    
    data.append({key[data_dic_key]: (channels[0].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
              channels[1].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1], 
              channels[3].data[int(slicing_lim_tdms[data_dic_key-1,0]):int(slicing_lim_tdms[data_dic_key,0])+1])})
        
    #Plot: Auslenkung über Frequenz
    plt.figure()
    #plt.plot(data_dic[key[data_dic_key]][0], data_dic[key[data_dic_key]][2], 'g')
    plt.plot(data[data_dic_key-1][key[data_dic_key]][0], data[data_dic_key-1][key[data_dic_key]][2], 'g')
    plt.title(key[data_dic_key] + '_Frequenzgang')
    plt.xlabel('Frequenz [Hz]')
    plt.ylabel('Auslenkung  [mm]')
    plt.grid()
    plt.tight_layout()
    plt.show()


## EXCEL-DATEIEN
    
    
# =============================================================================
# # Read Excel Data Function
# def read_excel_file(filename, sheet_index):
#     
#     book = xlrd.open_workbook(filename, encoding_override = "utf-8") #Öffnet Excel-Datei
#     
#     sheet_name=book.sheet_names()
#         
#     sheet = book.sheet_by_index(sheet_index) #Öffnet sheet
#     
#     data=np.zeros((sheet.nrows, sheet.ncols))
#     
#     for zeile in range(1,sheet.nrows): #Schleife über alle Zeilen des XL-Sheets
#         
#         for spalte in range(0,sheet.ncols): #Schleife über alle Spalten des XL-Sheets
#     
#             data[zeile, spalte] = sheet.cell(zeile, spalte).value #Speichern des Sheets in Matrix-Form
#             
#     return (data, sheet_name)
# 
# ## Importieren der Excel-Datei
# 
# DATA_FILE = "FG_Rot_mit_Deckel_final.xlsx" #Definition des filename
# 
# #os.makedirs("Img/" + "%s" % os.path.splitext(DATA_FILE)[0])
# 
# for sheet in range(0,3):
# 
#     data, sheet_name = read_excel_file(DATA_FILE, sheet) #Einlesen
#     
#     slicing_lim=np.zeros((int(np.max(data[:,1])/50)+1,1))
#         
#     for zaehler in range(50, 300, 50):
#     
#         for zaehler2 in range(0, np.size(data, 0)):
#             
#             if data[zaehler2, 1] == zaehler:
#     
#                 slicing_lim[int(zaehler/50)]=zaehler2
#                 
#     data_dic={"50V": data[int(slicing_lim[0])+1:int(slicing_lim[1])+1], "100V": data[int(slicing_lim[1])+1:int(slicing_lim[2])+1], "150V": data[int(slicing_lim[2])+1:int(slicing_lim[3])+1], "200V": data[int(slicing_lim[3])+1:int(slicing_lim[4])+1], "250V": data[int(slicing_lim[4])+1:int(slicing_lim[-1])+1]}
#         
#     Voltage={}
#     
#     
#     for zaehler3 in range(1, np.size(data,1)):
# 
#         Voltage[zaehler3]="%dV" %(zaehler3*50)
#         
#         plt.figure()
#         plt.plot(data_dic[Voltage[zaehler3]][:,0], data_dic[Voltage[zaehler3]][:,3], 'g')
#         plt.title('Frequenzgang')
#         plt.xlabel(r'$1/s$')
#         plt.ylabel(r'$mm$')
#         #plt.legend(loc='best')
#         plt.grid()
#         plt.tight_layout()
#         plt.show()
#         #plt.savefig("Img/" + os.path.splitext(DATA_FILE)[0] + "/%s_%s.png" % (Voltage[zaehler3], sheet_name[sheet]))
#         plt.savefig("Img/" + "/%s_%s.png" % (Voltage[zaehler3], sheet_name[sheet]))
# =============================================================================
# =============================================================================
# data_dic={"50V": (channels[0].data[int(slicing_lim_tdms[0,0]):int(slicing_lim_tdms[1,0])+1], channels[1].data[int(slicing_lim_tdms[0,0]):int(slicing_lim_tdms[1,0])+1], channels[3].data[int(slicing_lim_tdms[0,0]):int(slicing_lim_tdms[1,0])+1]), 
#           "100V": (channels[0].data[int(slicing_lim_tdms[1,0])+1:int(slicing_lim_tdms[2,0])+1], channels[1].data[int(slicing_lim_tdms[1,0])+1:int(slicing_lim_tdms[2,0])+1], channels[3].data[int(slicing_lim_tdms[1,0])+1:int(slicing_lim_tdms[2,0])+1]), 
#           "150V": (channels[0].data[int(slicing_lim_tdms[2,0])+1:int(slicing_lim_tdms[3,0])+1], channels[1].data[int(slicing_lim_tdms[2,0])+1:int(slicing_lim_tdms[3,0])+1], channels[3].data[int(slicing_lim_tdms[2,0])+1:int(slicing_lim_tdms[3,0])+1]), 
#           "200V": (channels[0].data[int(slicing_lim_tdms[3,0])+1:int(slicing_lim_tdms[4,0])+1], channels[1].data[int(slicing_lim_tdms[3,0])+1:int(slicing_lim_tdms[4,0])+1], channels[3].data[int(slicing_lim_tdms[3,0])+1:int(slicing_lim_tdms[4,0])+1]), 
#           "250V": (channels[0].data[int(slicing_lim_tdms[4,0])+1:int(slicing_lim_tdms[5,0])+1], channels[1].data[int(slicing_lim_tdms[4,0])+1:int(slicing_lim_tdms[5,0])+1], channels[3].data[int(slicing_lim_tdms[4,0])+1:int(slicing_lim_tdms[5,0])+1])
#           }
# =============================================================================

# =============================================================================
# voltage={}
#      
#      
# for zaehler3 in range(1, (int(np.max(channels[1].data)/50)+1)):
#  
#      voltage[zaehler3]="%dV" %(zaehler3*50)
#      
#      plt.figure()
#      plt.plot(data_dic[voltage[zaehler3]][0], data_dic[voltage[zaehler3]][2], 'g')
#      plt.title('Frequenzgang')
#      plt.xlabel(r'$1/s$')
#      plt.ylabel(r'$mm$')
#      #plt.legend(loc='best')
#      plt.grid()
#      plt.tight_layout()
#      plt.show()
#      #plt.savefig("Img/" + os.path.splitext(DATA_FILE)[0] + "/%s_%s.png" % (Voltage[zaehler3], sheet_name[sheet]))
#      #plt.savefig("Img/" + "/%s_%s.png" % (Voltage[zaehler3], sheet_name[sheet]))
# =============================================================================
