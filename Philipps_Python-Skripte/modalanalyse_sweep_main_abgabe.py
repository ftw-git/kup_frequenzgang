# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 16:59:19 2020

@author: phil SSD
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
from nptdms import TdmsFile as tdms
import tkinter as tk
from tkinter import filedialog
import os

'''Import Punktmessung Files im untergeordneten Ordner von main
Messungen 'PunktmessungX.tdms' benennen'''


# unterordner = 'Testdurchlauf9_5s_amp1'

# dataimport = [tdms.read(unterordner+'\Punktmessung1.tdms').groups()[0].channels(),
#               tdms.read(unterordner+'\Punktmessung2.tdms').groups()[0].channels(),
#               tdms.read(unterordner+'\Punktmessung3.tdms').groups()[0].channels(),
#               tdms.read(unterordner+'\Punktmessung4.tdms').groups()[0].channels(),
#               tdms.read(unterordner+'\Punktmessung5.tdms').groups()[0].channels(),
#               tdms.read(unterordner+'\Punktmessung6.tdms').groups()[0].channels()]

#%% Read data

print('Choose .tdms-Directory.')
root = tk.Tk()
root.withdraw()
unterordner = filedialog.askdirectory(title='Choose .tdms-Directory.') #Öffnet Dialogfenster zum auswählen der Datei

dataimport = []

for (dirpath, dirs, files) in os.walk(unterordner): #Schleife über Ordner und eventuelle Unterverzeichnisse

    for file in files:
        if os.path.splitext(file)[-1] == '.tdms':
            dataimport.append( 
                tdms.read(unterordner+'/{}'.format(file)).groups()[0].channels()
                )


'''Sweeplänge einstellen, Sample dt aus dataimport, variablen definieren'''
channelproperties = dataimport[0][0].properties
choosesensor = ''
sweeplen = 50000
sensitivity = 0
samplespersec = 10000
sampledt = round(list(dataimport[0][0].properties.items())[5][1],7)
resolution = samplespersec/sweeplen
channelid = 0
verstärkung = 100

'''Anregung (F) Kanal 0; '''
count = np.arange(0,6,1)
F =  np.empty(6, dtype=list)
Fpeaks = np.empty(6, dtype=list)
d2Fdx2 = np.empty(6, dtype=list)
Fshort = np.empty(6,dtype=list)
for i in count:
    F[i] = dataimport[i][0].data[0:100000]*verstärkung/100
    d2Fdx2[i] = np.gradient(np.gradient(F[i]))
    Fpeaks[i] = sp.find_peaks(np.diff(d2Fdx2[i]),np.max(np.diff(d2Fdx2[i]))*0.9)[0][0]
    Fshort[i] = F[i][Fpeaks[i]:Fpeaks[i]+sweeplen]

'''Anregung FFT [0Hz - 250Hz]'''
Ffft = np.empty(6, dtype=list)
for i in count:
    Ffft[i] = np.fft.fft(Fshort[i], norm = 'ortho')/(np.pi)
    Ffft[i] = Ffft[i][0:1250]

'''Amplitude berechnen und wählen des Sensors'''
setamp = round(max(F[0]),1)
if setamp > 1.99:
    choosesensor = 'L'
    sensitivity = 0.4 
    channelid = 1
else: 
    choosesensor = 'W' 
    sensitivity = -8
    channelid = 2

'''Messreihe; detrended; aufgeteilt in Druckseite und Saugseite wegen Richtung der Messung'''
X = np.empty(6, dtype=list)
Xshort = np.empty(6, dtype=list)
for i in np.arange(0,3,1):
    X[i] =  sp.detrend(dataimport[i][channelid].data)/(-1*sensitivity)
    Xshort[i] = X[i][Fpeaks[i]:(Fpeaks[i]+sweeplen)]
for i in np.arange(3,6,1):
    X[i] =  sp.detrend(dataimport[i][channelid].data)/(sensitivity)
    Xshort[i] = X[i][Fpeaks[i]:(Fpeaks[i]+sweeplen)]

'''Butterworth (Second Order Sections) Filter für Messungen'''
sosX = sp.butter(30, 250, output='sos', btype='low', fs=samplespersec)
for i in count:
        Xshort[i] = sp.sosfiltfilt(sosX, Xshort[i])

'''KO Transformation graphisch geschätzt'''
ko_cor = [1/np.cos(np.radians(-1.146)),1/np.cos(np.radians(-16.7)),1/np.cos(np.radians(-25.1)),1/np.cos(np.radians(-1.146)),1/np.cos(np.radians(-6.84)),1/np.cos(np.radians(-5))]
for i in count:
    Xshort[i] = Xshort[i]*ko_cor[i]

'''Fast Time Fourier Transform der Messung;
Bilden der Übertragungsfunktionen HxF;
Ermitteln der Eigenfrequenzen'''
Xfft = np.empty(6, dtype=list)
HxFeigvalues = np.empty(6, dtype=list)
HxFeigfrequencies = np.empty(6,dtype=list)
HxF = np.empty(6, dtype=list)
HxFeigvalue1 = np.empty(6, dtype=float)
HxFeigvalue2 = np.empty(6, dtype=float)
sosHxF = sp.butter(4, 65, output='sos', btype='low', fs=500)
HxF_abs_filt = np.empty(6, dtype=list)
for i in count:
    Xfft[i] = np.fft.rfft(Xshort[i], norm = 'ortho')/(2*np.pi)
    Xfft[i] = Xfft[i][0:1250]
    HxF[i] = (Xfft[i]/Ffft[i])[0:1250]
    if choosesensor == 'L':
        HxF_abs_filt[i] = sp.sosfiltfilt(sosHxF, np.abs(HxF[i]))
        HxFeigvalues[i] = sp.find_peaks(HxF_abs_filt[i][30:1200], np.max(HxF_abs_filt[i][30:1200])*0.2, distance = 70)[0]+30
        HxFeigvalue1[i] = HxFeigvalues[i][0]
    else:
        HxFeigvalues[i] = sp.find_peaks(np.abs(HxF[i][30:1200]), np.max(np.abs(HxF[i][30:1200]))*0.2, distance = 70)[0]+30
        HxFeigvalue1[i] = HxFeigvalues[i][0]
    if len(HxFeigvalues[i]) > 1:
        HxFeigvalue2[i] = HxFeigvalues[i][1]
    else:
        HxFeigvalue2[i] = 0
for i in count:
    if HxFeigvalue2[i] == 0:
        if choosesensor == 'L':
            HxFeigvalue2[i] = sp.find_peaks(np.abs(HxF_abs_filt[i][int(HxFeigvalue2[5]-50):int(HxFeigvalue2[5]+50)]), np.max(np.abs(HxF_abs_filt[i][int(HxFeigvalue2[5]-50):int(HxFeigvalue2[5]+50)]))*0.2, distance = 50)[0][0]+(HxFeigvalue2[5]-50)
            HxFeigvalues[i] = np.append(int(HxFeigvalues[i]), int(HxFeigvalue2[i]))
            HxFeigfrequencies[i] = HxFeigvalues[i]/(sweeplen*sampledt)
        else:
            HxFeigvalue2[i] = sp.find_peaks(np.abs(HxF[i][int(HxFeigvalue2[5]-50):int(HxFeigvalue2[5]+50)]), np.max(np.abs(HxF[i][int(HxFeigvalue2[5]-50):int(HxFeigvalue2[5]+50)]))*0.2, distance = 50)[0][0]+(HxFeigvalue2[5]-50)
            HxFeigvalues[i] = np.append(int(HxFeigvalues[i]), int(HxFeigvalue2[i]))
            HxFeigfrequencies[i] = HxFeigvalues[i]/(sweeplen*sampledt)
    else:
        HxFeigvalue2[i] = HxFeigvalue2[i]
        HxFeigfrequencies[i] = HxFeigvalues[i]/(sweeplen*sampledt)
HxFeigvaluecombined_1 = int(np.round(np.average(HxFeigvalue1[HxFeigvalue1 != 0]),0))
HxFeigfreqcombined_1 = np.round(HxFeigvaluecombined_1/(sweeplen*sampledt),1)
if np.sum(HxFeigvalue2) == 0:
    HxFeigvaluecombined_2 = 0
    HxFeigfreqcombined_2 = 0
else:
    HxFeigvaluecombined_2 = int(np.round(np.average(HxFeigvalue2[HxFeigvalue2 != 0]),0))
    HxFeigfreqcombined_2 = np.round(HxFeigvaluecombined_2/(sweeplen*sampledt), 1)
fftfreq = np.fft.fftfreq(sweeplen, d = sampledt)

'''Berechnen der Dämpfungsmaße der einzelnen Eigenfrequenzen für jeden Punkt'''
df1_1 = np.empty(6, dtype=int)
df2_1 = np.empty(6, dtype=int)
df1_2 = np.empty(6, dtype=int)
df2_2 = np.empty(6, dtype=int)
HxF1_1 = np.empty(6, dtype=list)
HxF2_1 = np.empty(6, dtype=list)
HxF1_2 = np.empty(6, dtype=list)
HxF2_2 = np.empty(6, dtype=list)    
D1 = np.empty(6, dtype=list)
D2 = np.empty(6, dtype=list)
HxFmaxes1 = np.zeros(6)
HxFmaxes2 = np.zeros(6)
HxFmaxes1sqrt = np.zeros(6)
HxFmaxes2sqrt = np.zeros(6)
for i in count:
    HxFmaxes1[i] = np.abs(HxF[i][HxFeigvalues[i][0]])
    HxFmaxes1sqrt[i] = HxFmaxes1[i]/np.sqrt(2)
    if choosesensor == 'L':
        HxFmaxes2[i] = np.abs(HxF_abs_filt[i][HxFeigvalues[i][1]])
        HxFmaxes2sqrt[i] = HxFmaxes2[i]/np.sqrt(2)
    else:
        HxFmaxes2[i] = np.abs(HxF[i][HxFeigvalues[i][1]])
        HxFmaxes2sqrt[i] = HxFmaxes2[i]/np.sqrt(2)
        continue  
for i in count:
    df1_1[i] = (np.abs(np.abs(HxF[i][30:HxFeigvalues[i][0]])-(np.abs(HxFmaxes1[i])/np.sqrt(2)))).argmin()+30
    df2_1[i] = (np.abs(np.abs(HxF[i][HxFeigvalues[i][0]:HxFeigvalues[i][0]+30])-(np.abs(HxFmaxes1[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][0]
    HxF1_1[i] = np.abs(HxF[i][df1_1[i]])
    HxF2_1[i] = np.abs(HxF[i][df2_1[i]])
    if choosesensor == 'L':       
        df1_2[i] = (np.abs(np.abs(HxF_abs_filt[i][HxFeigvalues[i][1]-20:HxFeigvalues[i][1]])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][1]-20
        df2_2[i] = (np.abs(np.abs(HxF_abs_filt[i][HxFeigvalues[i][1]:HxFeigvalues[i][1]+20])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][1]   
        HxF1_2[i] = np.abs(HxF_abs_filt[i][df1_2[i]])
        HxF2_2[i] = np.abs(HxF_abs_filt[i][df2_2[i]])
    else:
        df1_2[i] = (np.abs(np.abs(HxF[i][HxFeigvalues[i][1]-20:HxFeigvalues[i][1]])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][1]-20
        df2_2[i] = (np.abs(np.abs(HxF[i][HxFeigvalues[i][1]:HxFeigvalues[i][1]+20])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][1]   
        HxF1_2[i] = np.abs(HxF[i][df1_2[i]])
        HxF2_2[i] = np.abs(HxF[i][df2_2[i]])
for i in count:
    D1[i] = ((df2_1[i]/(sweeplen*sampledt)*2*np.pi - df1_1[i]/(sweeplen*sampledt)*2*np.pi)/((2*HxFeigvalues[i][0])/(sweeplen*sampledt)*2*np.pi))
    if len(HxFeigvalues[i]) > 1:
        D2[i] = ((df2_2[i]/(sweeplen*sampledt)*2*np.pi - df1_2[i]/(sweeplen*sampledt)*2*np.pi)/((2*HxFeigvalues[i][1])/(sweeplen*sampledt)*2*np.pi))
    else:
        D2[i] = 0

'''Interpolation zur Verfeinerung der Dämpfungsmaße'''    
x_new1l = np.empty(6, dtype=float)      
x_new1r = np.empty(6, dtype=float)      
D1_int = np.empty(6, dtype=float)
x_new2l = np.empty(6, dtype=float) 
x_new2r = np.empty(6, dtype=float)      
D2_int = np.empty(6, dtype=float)
for i in count:
    x_new1l[i] = np.interp(HxFmaxes1sqrt[i], [np.abs(HxF[i][df1_1[i]-2]), np.abs(HxF[i][df1_1[i]]), np.abs(HxF[i][df1_1[i]+2])], [df1_1[i]-2, df1_1[i], df1_1[i]+2])
    x_new1r[i] = np.interp(HxFmaxes1sqrt[i], [np.abs(HxF[i][df2_1[i]+2]), np.abs(HxF[i][df2_1[i]]), np.abs(HxF[i][df2_1[i]-2])], [df2_1[i]+2, df2_1[i], df2_1[i]-2])
    D1_int[i] = (x_new1r[i]/(sweeplen*sampledt)*2*np.pi-x_new1l[i]/(sweeplen*sampledt)*2*np.pi)/((2*HxFeigvalues[i][0])/(sweeplen*sampledt)*2*np.pi)
    if choosesensor == 'L':
        x_new2l[i] = np.interp(HxFmaxes2sqrt[i], [HxF_abs_filt[i][df1_2[i]-2], np.abs(HxF_abs_filt[i][df1_2[i]]), np.abs(HxF_abs_filt[i][df1_2[i]+2])], [df1_2[i]-2, df1_2[i], df1_2[i]+2])
        x_new2r[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df2_2[i]+2]), np.abs(HxF[i][df2_2[i]]), np.abs(HxF[i][df2_2[i]-2])], [df2_2[i]+2, df2_2[i], df2_2[i]-2])
        D2_int[i] = (x_new2r[i]/(sweeplen*sampledt)*2*np.pi-x_new2l[i]/(sweeplen*sampledt)*2*np.pi)/((2*HxFeigvalues[i][1])/(sweeplen*sampledt)*2*np.pi)
    else:
        x_new2l[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df1_2[i]-2]), np.abs(HxF[i][df1_2[i]]), np.abs(HxF[i][df1_2[i]+2])], [df1_2[i]-2, df1_2[i], df1_2[i]+2])
        x_new2r[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df2_2[i]+2]), np.abs(HxF[i][df2_2[i]]), np.abs(HxF[i][df2_2[i]-2])], [df2_2[i]+2, df2_2[i], df2_2[i]-2])
        D2_int[i] = (x_new2r[i]/(sweeplen*sampledt)*2*np.pi-x_new2l[i]/(sweeplen*sampledt)*2*np.pi)/((2*HxFeigvalues[i][1])/(sweeplen*sampledt)*2*np.pi)
    

'''Durchschnittliche Dämpfungsgrade berechnen
Eigenvektor bilden (physikalisch r_1 und r_2 sowie regelungstechnisch phi_1 und phi_2)'''
HxFcomplexval1 = np.empty(6, dtype=complex)
HxFcomplexval2 = np.empty(6, dtype=complex)
HxFangle1 = np.empty(6)
HxFangle2 = np.empty(6)
for i in count:
    HxFcomplexval1[i] = HxF[i][int(HxFeigvalue1[i])]
    HxFcomplexval2[i] = HxF[i][int(HxFeigvalue2[i])]
    HxFangle1[i] = np.angle(HxFcomplexval1[i], deg = True)
    HxFangle2[i] = np.angle(HxFcomplexval2[i], deg = True)

mode_1_dir = np.sign(HxFcomplexval1.imag)
mode_2_dir = np.sign(HxFcomplexval2.imag)
phi_1 = np.empty(6)
phi_2 = np.empty(6)
for i in count:
    phi_1[i] = HxFmaxes1[i]*mode_1_dir[i]
    phi_2[i] = HxFmaxes2[i]*mode_2_dir[i]
D1_avg = np.average(D1_int)
if np.sum(D2) == 0:
    D2_avg = 0
else:
    D2_avg = np.average(D2_int[D2_int != 0])
r_1 = np.empty(1, dtype=float)
r_1[0] = np.sqrt(np.abs((HxFcomplexval1[0])*(2*np.pi*HxFeigfrequencies[0][0])**2))*(D1_int[0]**4+4*D1_int[0]**2)**(1/4)
for i in np.arange(1,6,1):
    r_1 = np.append(r_1, np.abs((HxFcomplexval1[i])*(2*np.pi*HxFeigfrequencies[i][0])**2)/r_1[0]*np.sqrt(D1_int[i]**4+4*D1_int[i]**2))

r_2 = np.empty(1, dtype=float)
r_2[0] = np.sqrt(np.abs((HxFcomplexval2[0])*(2*np.pi*HxFeigfrequencies[0][1])**2))*(D2_int[0]**4+4*D2_int[0]**2)**(1/4)
for i in np.arange(1,6,1):
    r_2 = np.append(r_2, np.abs((HxFcomplexval2[i])*(2*np.pi*HxFeigfrequencies[i][1])**2)/r_1[0]*np.sqrt(D2_int[i]**4+4*D2_int[i]**2))
r_1 = r_1*mode_1_dir
r_2 = r_2*mode_2_dir
r_1_norm = r_1/np.max(np.abs(r_1))
r_2_norm = r_2/np.max(np.abs(r_2))


#%% PLOTS

plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.rc('axes', labelsize=14)

'''P - VARIABLE ZUR AUSWAHL DES ZU PLOTTENDEN PUNKTES'''
p = 0

'''PLOT 0: ANREGUNG UND MESSUNG IN PUNKT P'''

'''plot0, ad = plt.subplots(2,1)
plot0.set_figwidth(12)
plot0.set_figheight(8)
ad[0].set_ylabel('Amplitude in V')
ad[0].set_yscale("linear")
ad[0].plot(F[0]*100)
plot0.tight_layout()
ad[1].set_ylabel('Amplitude in mm')
ad[1].set_xlabel('Stützstelle')
ad[1].set_yscale("linear")
ad[1].plot(X[0])
#'''


'''PLOT 1: AMPLITUDEN- UND PHASENFREQUENZGANG DER ÜBERTRAGUNGSFUNKTION
Rote Punkte zeigen die Auswertung des gewählten Punktes,
Grüne Punkte zeigen die Auswertung der gemittelten Werte
Speicherfunktion ist auskommentiert'''

'''plot1, ax = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios': [3,1]})
plot1.set_figwidth(12)
plot1.set_figheight(8)
ax[0].set_ylabel('Amplitude in mm/100V')
ax[1].set_xlabel('Frequenz in Hz')
ax[0].set_yscale("log")
ax[0].set_title('Punktmessung '+ str(p+1))
plot1.tight_layout()
ax[0].grid(which = 'major', axis = 'y', linewidth = 0.25, linestyle='--')
if len(HxFeigvalues[p]) > 1: 
    ax[0].plot(fftfreq[HxFeigvalues[p][0]], np.abs(HxF[p])[HxFeigvalues[p][0]], 'r+', label = 'Punktmessung')
    ax[0].plot(fftfreq[HxFeigvalues[p][1]], np.abs(HxF[p])[HxFeigvalues[p][1]], 'r+')
else:
    ax[0].plot(fftfreq[HxFeigvalues[p][0]], np.abs(HxF[p])[HxFeigvalues[p][0]], 'r+', label = 'Punktmessung')
ax[0].plot(fftfreq[HxFeigvaluecombined_1], np.abs(HxF[p][HxFeigvaluecombined_1]), 'g*', label = 'Mittelwert')
ax[0].plot(fftfreq[HxFeigvaluecombined_2], np.abs(HxF[p][HxFeigvaluecombined_2]), 'g*')
ax[0].legend()
ax[1].set_ylabel('Phase')
ax[1].grid(which = 'major', axis = 'y', linewidth = 0.25, linestyle='--')
ax[1].set_yticks(np.arange(-180, 181, 90))
ax[1].plot(fftfreq[10:1000], np.angle(HxF[p][10:1000], deg=True))
ax[0].plot(fftfreq[10:1000], np.abs(HxF[p][10:1000]))
#plt.savefig(unterordner+'\Punktmessung '+str(p+1)+' bei '+str(setamp*100)+'V.svg', format = 'svg', dpi = 2000)
#'''

'''Punkte in [mm] ausgehend von unterer Einspannung'''
p1 = [40, 12.2]
p2 = [60.5, 6.75]
p3 = [75, 0.85]
p4 = [75, -6.15]
p5 = [49, -4.2]
p6 = [26, -2]

px = [p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]]
py = [p1[1], p2[1], p3[1], p4[1], p5[1], p6[1]]


'''PLOT 2: MODEN 1 UND 2 AUF SCHAUFEL 
Graphische Darstellung mit Hilfe der Pfeile für Richtung und Amplitude
Grünes Kreuz markiert den gewählten Nullpunkt'''

# plot2, ay = plt.subplots(2,1)
# imgschaufel = plt.imread('Schaufel.JPG')
# plot2.set_figwidth(12)
# plot2.set_figheight(8)
# ay[0].imshow(imgschaufel, extent = [82.1,-69.9, -15.8, int(152*451/1692)-15.8])
# ay[0].plot(0,0,'g+')
# ay[0].plot(px,py, 'r*')
# ay[0].set_title('Mode 1')
# ay[1].set_xlabel('x-Koordinate in mm')
# ay[1].set_ylabel('y-Koordinate in mm')
# plot2.tight_layout()
# plotscale = 10
# for i in count:
#     ay[0].arrow(px[i], py[i], 0, 1*plotscale*phi_1[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')

# ay[1].imshow(imgschaufel, extent = [82.1,-69.9, -15.8, int(152*451/1692)-15.8])
# ay[1].plot(0,0,'g+')
# ay[1].plot(px,py, 'r*')
# ay[1].set_title('Mode 2')
# for i in count: 
#     ay[1].arrow(px[i], py[i], 0, -1*plotscale*2*phi_2[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')
# #plt.savefig(unterordner+'\Moden.svg', format = 'svg', dpi = 2000)
# #'''


'''PLOT 3: VERGLEICH DER AMPLITUDENFREQUENZGÄNGE ALLER MESSPUNKTE'''

'''plot3, az = plt.subplots()
plot3.set_figwidth(12)
plot3.set_figheight(8)
az.set_ylabel('Komplexe Amplitude')
az.set_xlabel('Frequenz in Hz')
az.set_yscale("log")
plot3.tight_layout()
az.grid(which = 'major', axis = 'y', linewidth = 0.5, linestyle='--')
for i in count:
    az.plot(fftfreq[0:1250], np.abs(HxF[i]), label = 'Punkt '+str(i+1))
az.legend()
#plt.savefig(unterordner+'\Amplitudenfrequenzgänge.svg', format = 'svg', dpi = 2000)'''


'''PLOT 4: ANREGUNG, AUSLENKUNG UND ÜBERTRAGUNGSFUNKTION EINES MESSPUNKTES'''

'''plot4, ab = plt.subplots(3,1)
plot4.set_figwidth(12)
plot4.set_figheight(8)
ab[0].set_ylabel('Anregung in V')
ab[1].set_ylabel('Messung in V')
ab[2].set_ylabel('Komplexe Amplitude')
ab[0].set_xlabel('Zeit in s')
ab[1].set_xlabel('Zeit in s')
ab[2].set_xlabel('Frequenz in Hz')
ab[2].set_yscale('log')
plot4.tight_layout()
ab[0].plot(np.arange(0,50000,1)/10000,Fshort[p])
ab[1].plot(np.arange(0,50000,1)/10000,Xshort[p])
ab[2].plot(fftfreq[:1250],np.abs(HxF[p][:1250]))
#plt.savefig(unterordner+'\Punkt '+str(p+1)+' Anregung_Messung_Übertragungsfunktion_1.svg', format = 'svg', dpi = 2000)'''


'''PLOT 5: EINZELVERGLEICH ALLER ÜBERTRAGUNGSFUNKTIONEN'''

'''plot5, ac = plt.subplots(2,3, sharey=True)
plot5.set_figwidth(12)
plot5.set_figheight(8)
plot5.tight_layout()
for i in np.arange(0,3,1):
    ac[0,i].plot(fftfreq[HxFeigvalues[i][0]], np.abs(HxF[i])[HxFeigvalues[i][0]], 'r+', label = 'Punktmessung')
    ac[0,i].plot(fftfreq[HxFeigvalues[i][1]], np.abs(HxF[i])[HxFeigvalues[i][1]], 'r+')
    ac[0,i].plot(fftfreq[HxFeigvaluecombined_1], np.abs(HxF[i][HxFeigvaluecombined_1]), 'g*', label = 'Mittelwert')
    ac[0,i].set_yscale('log')
    ac[0,i].plot(fftfreq[0:1250], np.abs(HxF[i]))
    ac[0,i].set_title('Punkt '+str(i+1))
    ac[0,1].set_ylabel('Komplexe Amplitude')
for i in np.arange(3,6,1):
    ac[1,i-3].plot(fftfreq[HxFeigvalues[i][0]], np.abs(HxF[i])[HxFeigvalues[i][0]], 'r+', label = 'Punktmessung')
    ac[1,i-3].plot(fftfreq[HxFeigvalues[i][1]], np.abs(HxF[i])[HxFeigvalues[i][1]], 'r+')
    ac[1,i-3].plot(fftfreq[HxFeigvaluecombined_1], np.abs(HxF[i][HxFeigvaluecombined_1]), 'g*', label = 'Mittelwert')
    ac[1,i-3].set_yscale('log')
    ac[1,i-3].plot(fftfreq[0:1250], np.abs(HxF[i]))
    ac[1,i-3].set_title('Punkt '+str(i+1))
    ac[1,1].set_xlabel('Frequenz [Hz]')
#plt.savefig(unterordner+'\Übersicht_Amplitudenfrequenzgang.svg', format = 'svg', dpi = 2000)'''


'''PLOT 6: KONTROLLE DES ERMITTELTEN DÄMPFUNGSMAßES'''

'''plot61, D_intplot = plt.subplots()
plt.ylabel('Komplexe Amplitude')
plt.xlabel('Frequenz in Hz')
plt.title('Berechnung des Dämpfungsgrades')
plt.yscale('log')
plt.plot(fftfreq[5:1000] ,np.abs(HxF[p])[5:1000])
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new1l[p],4)*10000/(sweeplen*sampledt))], HxFmaxes1sqrt[p], 'r+')
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new1r[p],4)*10000/(sweeplen*sampledt))], HxFmaxes1sqrt[p], 'r+')
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new2l[p],4)*10000/(sweeplen*sampledt))], HxFmaxes2sqrt[p], 'r+')
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new2r[p],4)*10000/(sweeplen*sampledt))], HxFmaxes2sqrt[p], 'r+')

plot62, Dplot = plt.subplots()
plt.ylabel('Komplexe Amplitude')
plt.xlabel('Frequenz in Hz')
plt.title('Berechnung des Dämpfungsgrades')
plt.yscale('log')
plt.plot(fftfreq[5:1000] ,np.abs(HxF[p])[5:1000])
plt.plot(fftfreq[df1_1[p]], np.abs(HxF1_1[p]), 'r+')
plt.plot(fftfreq[df2_1[p]], np.abs(HxF2_1[p]), 'r+')
plt.plot(fftfreq[df1_2[p]], np.abs(HxF1_2[p]), 'r+')
plt.plot(fftfreq[df2_2[p]], np.abs(HxF2_2[p]), 'r+')
#'''


'''Savefiles für tabellarische Auswertung'''
# from tabulate import tabulate as tb
# import pandas as pd


# '''datasave = open(unterordner+'\Messung Stats.txt', 'w')
# datasave.write('Erste Mode bei '+str(HxFeigfreqcombined_1)+'Hz (max Abweichung '+str((np.max(HxFeigvalue1)-np.min(HxFeigvalue1))/5)+'Hz)\n'+'Modale Daempfung '+str(round(D1_avg,5))+' \n'+'Zweite Mode bei '+str(+HxFeigfreqcombined_2)+'Hz (max Abweichung '+str((np.max(HxFeigvalue2)-np.min(HxFeigvalue2[HxFeigvalue2 != 0]))/5)+'Hz)\n'+'Modale Daempfung '+str(round(D2_avg,5))+'\n'+'Messauflösung '+str(resolution)+'Hz')
# datasave.close()'''
# table = [['Messpunkt', 1,2,3,4,5,6], 
#          ['-','-','-','-','-','-','-'],
#          ['Erste Eigenfrequenz',HxFeigfrequencies[0][0],HxFeigfrequencies[1][0],HxFeigfrequencies[2][0],HxFeigfrequencies[3][0],HxFeigfrequencies[4][0],HxFeigfrequencies[5][0]],
#          ['Dämpfungsmaß 1te Mode', np.round(D1_int[0],5),np.round(D1_int[1],5),np.round(D1_int[2],5),np.round(D1_int[3],5),np.round(D1_int[4],5),np.round(D1_int[5],5)],
#          ['Auslenkung 1te Mode', np.round(phi_1[0],5), np.round(phi_1[1],5), np.round(phi_1[2],5), np.round(phi_1[3],5), np.round(phi_1[4],5), np.round(phi_1[5],5)],
#          ['-','-','-','-','-','-','-'],
#          ['Zweite Eigenfrequenz',HxFeigfrequencies[0][1],HxFeigfrequencies[1][1],HxFeigfrequencies[2][1],HxFeigfrequencies[3][1],HxFeigfrequencies[4][1],HxFeigfrequencies[5][1]],
#          ['Dämpfungsmaß 2te Mode', np.round(D2_int[0],5),np.round(D2_int[1],5),np.round(D2_int[2],5),np.round(D2_int[3],5),np.round(D2_int[4],5),np.round(D2_int[5],5)],
#          ['Auslenkung 2te Mode', np.round(phi_2[0],5), np.round(phi_2[1],5), np.round(phi_2[2],5), np.round(phi_2[3],5), np.round(phi_2[4],5), np.round(phi_2[5],5)],
#         ]
# print(tb(table))
# print('Erste Mode bei '+str(HxFeigfreqcombined_1)+'Hz\n'+'Modale Daempfung '+str(round(D1_avg,5))+' \n'+'Zweite Mode bei '+str(+HxFeigfreqcombined_2)+'Hz\n'+'Modale Daempfung '+str(round(D2_avg,5))) 

# dataframe = pd.DataFrame(table)
# dataframe.to_excel(unterordner+'/auswertung_'+str(int(setamp*100))+'V.xlsx', index=False, header=False)
# #'''

'''plot7, hxfplot = plt.subplots()
plt.xlabel('Frequenz in Hz')
plt.ylabel('Amplitude in mm/100V')
if choosesensor == 'W':
    plt.plot(fftfreq[0:1250], np.abs(HxF[0]), color='red')
    plt.plot(fftfreq[0:1250], np.abs(HxF[1]), color='red')
    plt.plot(fftfreq[0:1250], np.abs(HxF[2]), color='red')
    plt.plot(fftfreq[0:1250], np.abs(HxF[3]), color='green')
    plt.plot(fftfreq[0:1250], np.abs(HxF[4]), color='green')
    plt.plot(fftfreq[0:1250], np.abs(HxF[5]), color='green')
else:
    plt.plot(fftfreq[0:1250], np.abs(HxF_abs_filt[0]), color='red')
    plt.plot(fftfreq[0:1250], np.abs(HxF_abs_filt[1]), color='red')
    plt.plot(fftfreq[0:1250], np.abs(HxF_abs_filt[2]), color='red')
    plt.plot(fftfreq[0:1250], np.abs(HxF_abs_filt[3]), color='green')
    plt.plot(fftfreq[0:1250], np.abs(HxF_abs_filt[4]), color='green')
    plt.plot(fftfreq[0:1250], np.abs(HxF_abs_filt[5]), color='green')
plt.yscale('log')
plot7.tight_layout
plot7.set_figwidth(10)
plot7.set_figheight(6)
plt.savefig(unterordner+'\Übersicht_Amplitudenfrequenzgang.jpeg', format = 'jpeg', dpi = 150)
#'''

