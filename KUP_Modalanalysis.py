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

#%% Read data

# Asks for tdms-directory
print('Choose .tdms-Directory.')
root = tk.Tk()
root.withdraw()
sub_dic = filedialog.askdirectory(title='Choose .tdms-Directory.') #Öffnet Dialogfenster zum auswählen der Datei

# Empty list to store all tdms-files
data_import = []

# Walks through all files in sub directory
for (dirpath, dirs, files) in os.walk(sub_dic): #Schleife über Ordner und eventuelle Unterverzeichnisse

    for file in files:
        if os.path.splitext(file)[-1] == '.tdms':
            data_import.append( tdms.read(sub_dic+'/{}'.format(file)).groups()[0].channels() )

#%% Data preprocessing

# Extracting channel properties from first tdms-file
channel_prop = data_import[0][0].properties
channel_header = [ key for key in channel_prop.keys() ]

# Constants of sweep
sweep_len = 50000
sensitivity = 0
sample_freq = 10000
resolution = sample_freq/sweep_len
gain = 100

# Sample time increment from channel properties, rounded to 7th decimal
sampledt = round(channel_prop['wf_increment'], 7)

'''Anregung (F) Kanal 0 '''

# Number of tdms-files to be evaluated
num_of_files = 6#np.arange(0, len(data_import))
# Initialize simulation F
F = np.empty(num_of_files, dtype=list)
# Indices of peaks of stimulation
F_peaks = np.empty(num_of_files, dtype=list)
# Hessian of stimulation
d2Fdx2 = np.empty(num_of_files, dtype=list)
# ???
F_short = np.empty(num_of_files, dtype=list)
# FFT of stimulation
F_fft = np.empty(num_of_files, dtype=list)

# Loop over all files
for i in range(0, num_of_files):
    # Get data from data_import
    F[i] = data_import[i][0].data[0:100000]*gain/100
    # Calculate hessian
    d2Fdx2[i] = np.gradient(np.gradient(F[i]))
    # Finds peak, height: required height for identification as peak, here 90% of max
    # Only first peak is taken into account
    F_peaks[i] = sp.find_peaks( np.diff(d2Fdx2[i]), 
                                height = np.max(np.diff(d2Fdx2[i]))*0.9 )[0][0]
    # Stimulation is shortened to start at first peak and end at length of sweep
    F_short[i] = F[i][F_peaks[i]:F_peaks[i]+sweep_len]

    # FFT of stimulation for frequencies 0Hz - 250Hz
    F_fft[i] = np.fft.fft(F_short[i], norm = 'ortho')/(np.pi)
    F_fft[i] = F_fft[i][0:1250]

# Calculate amplitude and decide on sensor to be used
# Get amplitude from max of stimulation
set_amp = round( max( F[0] ), 1 )

# Decides on sensor to evaluate depending on amplitude
if set_amp > 1.99:
    choosesensor = 'L'
    sensitivity = 0.4 
    channel_id = 1
else: 
    choosesensor = 'W' 
    sensitivity = -8
    channel_id = 2

# Removing linear trend from measurements and separating data into SS and DS
# Iinitialize
X = np.empty(6, dtype=list)
X_short = np.empty(6, dtype=list)

# Loop over all files
for i in range(0, num_of_files):
    # Separates datasets into Saugseite and Druckseite due to direction of measurements
    if i < int(num_of_files/2):
        X[i] =  sp.detrend(data_import[i][channel_id].data)/(-1*sensitivity)
    else:
        X[i] =  sp.detrend(data_import[i][channel_id].data)/(sensitivity)
    # Shortening measurements in analogy to stimulation
    X_short[i] = X[i][F_peaks[i]:(F_peaks[i]+sweep_len)]

# Butterworth-Filter for measurements
N = 30 # Order of filter
Wn = 250 # Critical frequency
# output-type: sos, second-order-sections for general filter purposes
# fs = sampling frequency of digital system
butter_X = sp.butter(N, Wn, output = 'sos', btype = 'low', fs = sample_freq)
for i in range(0, num_of_files):
        X_short[i] = sp.sosfiltfilt(butter_X, X_short[i])

# KO-Transformation, graphically estimated ???
ko_cor = [ 1/np.cos(np.radians(-1.146)), 
          1/np.cos(np.radians(-16.7)), 
          1/np.cos(np.radians(-25.1)), 
          1/np.cos(np.radians(-1.146)), 
          1/np.cos(np.radians(-6.84)), 
          1/np.cos(np.radians(-5)) ]
# Transforms over all files
for i in range(0, num_of_files):
    X_short[i] = X_short[i]*ko_cor[i]


#%% FFT of measurements and calculation of Transfer Function, Determination of Eigen-Frequencies

X_fft = np.empty(6, dtype=list)
HxF_eigval = np.empty(6, dtype=list)
HxF_eigfreq = np.empty(6,dtype=list)
HxF = np.empty(6, dtype=list)
HxF_eigval1 = np.empty(6, dtype=float)
HxF_eigval2 = np.empty(6, dtype=float)
butter_HxF = sp.butter(4, 65, output='sos', btype='low', fs=500)
HxF_abs_filt = np.empty(6, dtype=list)


###############################################################################
# Loop over all files
for i in range(0, num_of_files):
    # DFFT of shortened measurements
    X_fft[i] = np.fft.rfft(X_short[i], norm = 'ortho')/(2*np.pi)
    X_fft[i] = X_fft[i][0:1250]
    
    # Transfer function from output/input
    HxF[i] = (X_fft[i]/F_fft[i])[0:1250]
    
    if choosesensor == 'L':
        HxF_abs_filt[i] = sp.sosfiltfilt(butter_HxF, np.abs(HxF[i]))
        HxF_eigval[i] = sp.find_peaks(HxF_abs_filt[i][30:1200], np.max(HxF_abs_filt[i][30:1200])*0.2, distance = 70)[0]+30
        HxF_eigval1[i] = HxF_eigval[i][0]
    else:
        HxF_eigval[i] = sp.find_peaks(np.abs(HxF[i][30:1200]), np.max(np.abs(HxF[i][30:1200]))*0.2, distance = 70)[0]+30
        HxF_eigval1[i] = HxF_eigval[i][0]
    if len(HxF_eigval[i]) > 1:
        HxF_eigval2[i] = HxF_eigval[i][1]
    else:
        HxF_eigval2[i] = 0
        
# Loop over all files
for i in range(0, num_of_files):
    if HxF_eigval2[i] == 0:
        if choosesensor == 'L':
            HxF_eigval2[i] = sp.find_peaks(np.abs(HxF_abs_filt[i][int(HxF_eigval2[5]-50):int(HxF_eigval2[5]+50)]), np.max(np.abs(HxF_abs_filt[i][int(HxF_eigval2[5]-50):int(HxF_eigval2[5]+50)]))*0.2, distance = 50)[0][0]+(HxF_eigval2[5]-50)
            HxF_eigval[i] = np.append(int(HxF_eigval[i]), int(HxF_eigval2[i]))
            HxF_eigfreq[i] = HxF_eigval[i]/(sweep_len*sampledt)
        else:
            HxF_eigval2[i] = sp.find_peaks(np.abs(HxF[i][int(HxF_eigval2[5]-50):int(HxF_eigval2[5]+50)]), np.max(np.abs(HxF[i][int(HxF_eigval2[5]-50):int(HxF_eigval2[5]+50)]))*0.2, distance = 50)[0][0]+(HxF_eigval2[5]-50)
            HxF_eigval[i] = np.append(int(HxF_eigval[i]), int(HxF_eigval2[i]))
            HxF_eigfreq[i] = HxF_eigval[i]/(sweep_len*sampledt)
    else:
        HxF_eigval2[i] = HxF_eigval2[i]
        HxF_eigfreq[i] = HxF_eigval[i]/(sweep_len*sampledt)

# Mean of all first eigen values
HxF_eigval1_combined = int( np.round( np.average(HxF_eigval1[HxF_eigval1 != 0]), 0 ) )
# Calculating and rounding first eigen frequency
HxF_eigfreq1_combined = np.round( HxF_eigval1_combined/(sweep_len*sampledt), 1 )

if np.sum(HxF_eigval2) == 0:
    HxFeigvaluecombined_2 = 0
    HxFeigfreqcombined_2 = 0
else:
    HxFeigvaluecombined_2 = int(np.round(np.average(HxF_eigval2[HxF_eigval2 != 0]),0))
    HxFeigfreqcombined_2 = np.round(HxFeigvaluecombined_2/(sweep_len*sampledt), 1)
fft_freq = np.fft.fftfreq(sweep_len, d = sampledt)

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
for i in range(0, num_of_files):
    HxFmaxes1[i] = np.abs(HxF[i][HxF_eigval[i][0]])
    HxFmaxes1sqrt[i] = HxFmaxes1[i]/np.sqrt(2)
    if choosesensor == 'L':
        HxFmaxes2[i] = np.abs(HxF_abs_filt[i][HxF_eigval[i][1]])
        HxFmaxes2sqrt[i] = HxFmaxes2[i]/np.sqrt(2)
    else:
        HxFmaxes2[i] = np.abs(HxF[i][HxF_eigval[i][1]])
        HxFmaxes2sqrt[i] = HxFmaxes2[i]/np.sqrt(2)
        continue  
for i in range(0, num_of_files):
    df1_1[i] = (np.abs(np.abs(HxF[i][30:HxF_eigval[i][0]])-(np.abs(HxFmaxes1[i])/np.sqrt(2)))).argmin()+30
    df2_1[i] = (np.abs(np.abs(HxF[i][HxF_eigval[i][0]:HxF_eigval[i][0]+30])-(np.abs(HxFmaxes1[i])/np.sqrt(2)))).argmin()+HxF_eigval[i][0]
    HxF1_1[i] = np.abs(HxF[i][df1_1[i]])
    HxF2_1[i] = np.abs(HxF[i][df2_1[i]])
    if choosesensor == 'L':       
        df1_2[i] = (np.abs(np.abs(HxF_abs_filt[i][HxF_eigval[i][1]-20:HxF_eigval[i][1]])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxF_eigval[i][1]-20
        df2_2[i] = (np.abs(np.abs(HxF_abs_filt[i][HxF_eigval[i][1]:HxF_eigval[i][1]+20])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxF_eigval[i][1]   
        HxF1_2[i] = np.abs(HxF_abs_filt[i][df1_2[i]])
        HxF2_2[i] = np.abs(HxF_abs_filt[i][df2_2[i]])
    else:
        df1_2[i] = (np.abs(np.abs(HxF[i][HxF_eigval[i][1]-20:HxF_eigval[i][1]])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxF_eigval[i][1]-20
        df2_2[i] = (np.abs(np.abs(HxF[i][HxF_eigval[i][1]:HxF_eigval[i][1]+20])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxF_eigval[i][1]   
        HxF1_2[i] = np.abs(HxF[i][df1_2[i]])
        HxF2_2[i] = np.abs(HxF[i][df2_2[i]])
for i in range(0, num_of_files):
    D1[i] = ((df2_1[i]/(sweep_len*sampledt)*2*np.pi - df1_1[i]/(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[i][0])/(sweep_len*sampledt)*2*np.pi))
    if len(HxF_eigval[i]) > 1:
        D2[i] = ((df2_2[i]/(sweep_len*sampledt)*2*np.pi - df1_2[i]/(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[i][1])/(sweep_len*sampledt)*2*np.pi))
    else:
        D2[i] = 0

'''Interpolation zur Verfeinerung der Dämpfungsmaße'''    
x_new1l = np.empty(6, dtype=float)      
x_new1r = np.empty(6, dtype=float)      
D1_int = np.empty(6, dtype=float)
x_new2l = np.empty(6, dtype=float) 
x_new2r = np.empty(6, dtype=float)      
D2_int = np.empty(6, dtype=float)
for i in range(0, num_of_files):
    x_new1l[i] = np.interp(HxFmaxes1sqrt[i], [np.abs(HxF[i][df1_1[i]-2]), np.abs(HxF[i][df1_1[i]]), np.abs(HxF[i][df1_1[i]+2])], [df1_1[i]-2, df1_1[i], df1_1[i]+2])
    x_new1r[i] = np.interp(HxFmaxes1sqrt[i], [np.abs(HxF[i][df2_1[i]+2]), np.abs(HxF[i][df2_1[i]]), np.abs(HxF[i][df2_1[i]-2])], [df2_1[i]+2, df2_1[i], df2_1[i]-2])
    D1_int[i] = (x_new1r[i]/(sweep_len*sampledt)*2*np.pi-x_new1l[i]/(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[i][0])/(sweep_len*sampledt)*2*np.pi)
    if choosesensor == 'L':
        x_new2l[i] = np.interp(HxFmaxes2sqrt[i], [HxF_abs_filt[i][df1_2[i]-2], np.abs(HxF_abs_filt[i][df1_2[i]]), np.abs(HxF_abs_filt[i][df1_2[i]+2])], [df1_2[i]-2, df1_2[i], df1_2[i]+2])
        x_new2r[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df2_2[i]+2]), np.abs(HxF[i][df2_2[i]]), np.abs(HxF[i][df2_2[i]-2])], [df2_2[i]+2, df2_2[i], df2_2[i]-2])
        D2_int[i] = (x_new2r[i]/(sweep_len*sampledt)*2*np.pi-x_new2l[i]/(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[i][1])/(sweep_len*sampledt)*2*np.pi)
    else:
        x_new2l[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df1_2[i]-2]), np.abs(HxF[i][df1_2[i]]), np.abs(HxF[i][df1_2[i]+2])], [df1_2[i]-2, df1_2[i], df1_2[i]+2])
        x_new2r[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df2_2[i]+2]), np.abs(HxF[i][df2_2[i]]), np.abs(HxF[i][df2_2[i]-2])], [df2_2[i]+2, df2_2[i], df2_2[i]-2])
        D2_int[i] = (x_new2r[i]/(sweep_len*sampledt)*2*np.pi-x_new2l[i]/(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[i][1])/(sweep_len*sampledt)*2*np.pi)
    

'''Durchschnittliche Dämpfungsgrade berechnen
Eigenvektor bilden (physikalisch r_1 und r_2 sowie regelungstechnisch phi_1 und phi_2)'''
HxFcomplexval1 = np.empty(6, dtype=complex)
HxFcomplexval2 = np.empty(6, dtype=complex)
HxFangle1 = np.empty(6)
HxFangle2 = np.empty(6)
for i in range(0, num_of_files):
    HxFcomplexval1[i] = HxF[i][int(HxF_eigval1[i])]
    HxFcomplexval2[i] = HxF[i][int(HxF_eigval2[i])]
    HxFangle1[i] = np.angle(HxFcomplexval1[i], deg = True)
    HxFangle2[i] = np.angle(HxFcomplexval2[i], deg = True)

mode_1_dir = np.sign(HxFcomplexval1.imag)
mode_2_dir = np.sign(HxFcomplexval2.imag)
phi_1 = np.empty(6)
phi_2 = np.empty(6)
for i in range(0, num_of_files):
    phi_1[i] = HxFmaxes1[i]*mode_1_dir[i]
    phi_2[i] = HxFmaxes2[i]*mode_2_dir[i]
D1_avg = np.average(D1_int)
if np.sum(D2) == 0:
    D2_avg = 0
else:
    D2_avg = np.average(D2_int[D2_int != 0])
r_1 = np.empty(1, dtype=float)
r_1[0] = np.sqrt(np.abs((HxFcomplexval1[0])*(2*np.pi*HxF_eigfreq[0][0])**2))*(D1_int[0]**4+4*D1_int[0]**2)**(1/4)
for i in np.arange(1,6,1):
    r_1 = np.append(r_1, np.abs((HxFcomplexval1[i])*(2*np.pi*HxF_eigfreq[i][0])**2)/r_1[0]*np.sqrt(D1_int[i]**4+4*D1_int[i]**2))

r_2 = np.empty(1, dtype=float)
r_2[0] = np.sqrt(np.abs((HxFcomplexval2[0])*(2*np.pi*HxF_eigfreq[0][1])**2))*(D2_int[0]**4+4*D2_int[0]**2)**(1/4)
for i in np.arange(1,6,1):
    r_2 = np.append(r_2, np.abs((HxFcomplexval2[i])*(2*np.pi*HxF_eigfreq[i][1])**2)/r_1[0]*np.sqrt(D2_int[i]**4+4*D2_int[i]**2))
r_1 = r_1*mode_1_dir
r_2 = r_2*mode_2_dir
r_1_norm = r_1/np.max(np.abs(r_1))
r_2_norm = r_2/np.max(np.abs(r_2))

###############################################################################

#%% PLOTS

plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.rc('axes', labelsize=14)

# Declare x-axis as time
def x_axis(y_vals, sample_freq):
 return np.arange(0, len(y_vals)/sample_freq, step=1/sample_freq)

# Choose the point to be evaluated
while True:    
    try:
        p = int(input('Which point out of the {} points should be evaluated\n?'.format(num_of_files)))
        if 0 < p-1 < num_of_files :
            break
    except (ValueError, TypeError):
        pass

#%% Plot 0: Excitation and measurements in point P

# plot_0, ax_0 = plt.subplots(2,1)
# plot_0.set_figwidth(12)
# plot_0.set_figheight(8)
# # Excitation/Stimulation
# ax_0[0].set_ylabel('Amplitude in V')
# ax_0[0].set_yscale("linear")
# ax_0[0].plot(x_axis(F[0], sample_freq), F[0]*100)
# # Measurements in P
# ax_0[1].set_ylabel('Amplitude in mm')
# ax_0[1].set_xlabel('Stützstelle')
# ax_0[1].set_yscale("linear")
# ax_0[1].plot(x_axis(X[0], sample_freq), X[0])
# plot_0.tight_layout()

#%% Plot 1: Amplitude- and Phaseresponse of Transfer-Function

# # Red dots: Evaluation of p
# # Green dots: Evaluation of averaged values

# plot_1, ax_1 = plt.subplots( 2, 1, sharex=True, gridspec_kw={'height_ratios': [3,1]} )
# plot_1.set_figwidth(12)
# plot_1.set_figheight(8)
# ax_1[0].set_ylabel('Amplitude in mm/100V')
# ax_1[1].set_xlabel('Frequenz in Hz')
# ax_1[0].set_yscale("log")
# ax_1[0].set_title('Punktmessung' + str(p+1) + ' - Übertragungsfunktion')
# plot_1.tight_layout()
# ax_1[0].grid(which = 'major', axis = 'y', linewidth = 0.25, linestyle='--')
# if len(HxF_eigval[p]) > 1: 
#     ax_1[0].plot(fft_freq[HxF_eigval[p][0]], np.abs(HxF[p])[HxF_eigval[p][0]], 'r+', label = 'Punktmessung')
#     ax_1[0].plot(fft_freq[HxF_eigval[p][1]], np.abs(HxF[p])[HxF_eigval[p][1]], 'r+')
# else:
#     ax_1[0].plot(fft_freq[HxF_eigval[p][0]], np.abs(HxF[p])[HxF_eigval[p][0]], 'r+', label = 'Punktmessung')
# ax_1[0].plot(fft_freq[HxF_eigval1_combined], np.abs(HxF[p][HxF_eigval1_combined]), 'g*', label = 'Mittelwert')
# ax_1[0].plot(fft_freq[HxFeigvaluecombined_2], np.abs(HxF[p][HxFeigvaluecombined_2]), 'g*')
# ax_1[0].legend()
# ax_1[1].set_ylabel('Phase')
# ax_1[1].grid(which = 'major', axis = 'y', linewidth = 0.25, linestyle='--')
# ax_1[1].set_yticks(np.arange(-180, 181, 90))
# ax_1[1].plot(fft_freq[10:1000], np.angle(HxF[p][10:1000], deg=True))
# ax_1[0].plot(fft_freq[10:1000], np.abs(HxF[p][10:1000]))
# # Save plot in measurement directory
# plt.savefig( sub_dic + str(p+1) + 'P'  + str(set_amp*100) + 'V.png', dpi = 2000)

#%% Plot 2: First and second mode plotted on blade
# ??? Punkte müssen neu vermessen werden

# Green dot: selected null point

# '''Punkte in [mm] ausgehend von unterer Einspannung'''
# p1 = [40, 12.2]
# p2 = [60.5, 6.75]
# p3 = [75, 0.85]
# p4 = [75, -6.15]
# p5 = [49, -4.2]
# p6 = [26, -2]

# px = [p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]]
# py = [p1[1], p2[1], p3[1], p4[1], p5[1], p6[1]]

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
# for i in num_of_files:
#     ay[0].arrow(px[i], py[i], 0, 1*plotscale*phi_1[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')

# ay[1].imshow(imgschaufel, extent = [82.1,-69.9, -15.8, int(152*451/1692)-15.8])
# ay[1].plot(0,0,'g+')
# ay[1].plot(px,py, 'r*')
# ay[1].set_title('Mode 2')
# for i in num_of_files: 
#     ay[1].arrow(px[i], py[i], 0, -1*plotscale*2*phi_2[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')
# #plt.savefig(sub_dic+'\Modes.png', dpi = 2000)

#%% Plot 3: Comparison Bode-Diagrams all Measurement-Points

# plot3, az = plt.subplots()
# plot3.set_figwidth(12)
# plot3.set_figheight(8)
# az.set_ylabel('Komplexe Amplitude')
# az.set_xlabel('Frequenz in Hz')
# az.set_yscale("log")
# plot3.tight_layout()
# az.grid(which = 'major', axis = 'y', linewidth = 0.5, linestyle='--')
# for i in num_of_files:
#     az.plot(fft_freq[0:1250], np.abs(HxF[i]), label = 'Punkt '+str(i+1))
# az.legend()
# #plt.savefig(sub_dic+'\Amplitudenfrequenzgänge.svg', format = 'svg', dpi = 2000)

#%% Plot 4: Excitation, Deflection and Transfer Function of one measurement point

plot_4, ax_4 = plt.subplots(3,1)
plot_4.set_figwidth(12)
plot_4.set_figheight(8)
ax_4[0].set_ylabel('Anregung in V')
ax_4[1].set_ylabel('Messung in V')
ax_4[2].set_ylabel('Komplexe Amplitude')
ax_4[0].set_xlabel('Zeit in s')
ax_4[1].set_xlabel('Zeit in s')
ax_4[2].set_xlabel('Frequenz in Hz')
ax_4[2].set_yscale('log')
plot_4.tight_layout()
ax_4[0].plot(np.arange(0,50000,1)/10000,F_short[p])
ax_4[1].plot(np.arange(0,50000,1)/10000,X_short[p])
ax_4[2].plot(fft_freq[:1250],np.abs(HxF[p][:1250]))
#plt.savefig(sub_dic+'\Punkt '+str(p+1)+' Anregung_Messung_Übertragungsfunktion_1.svg', format = 'svg', dpi = 2000)'''

#%% Plot 5: Single comparison of all Transfer Function
'''PLOT 5: EINZELVERGLEICH ALLER ÜBERTRAGUNGSFUNKTIONEN'''

# plot5, ac = plt.subplots(2,3, sharey=True)
# plot5.set_figwidth(12)
# plot5.set_figheight(8)
# plot5.tight_layout()
# for i in np.arange(0,3,1):
#     ac[0,i].plot(fft_freq[HxF_eigval[i][0]], np.abs(HxF[i])[HxF_eigval[i][0]], 'r+', label = 'Punktmessung')
#     ac[0,i].plot(fft_freq[HxF_eigval[i][1]], np.abs(HxF[i])[HxF_eigval[i][1]], 'r+')
#     ac[0,i].plot(fft_freq[HxF_eigval1_combined], np.abs(HxF[i][HxF_eigval1_combined]), 'g*', label = 'Mittelwert')
#     ac[0,i].set_yscale('log')
#     ac[0,i].plot(fft_freq[0:1250], np.abs(HxF[i]))
#     ac[0,i].set_title('Punkt '+str(i+1))
#     ac[0,1].set_ylabel('Komplexe Amplitude')
# for i in np.arange(3,6,1):
#     ac[1,i-3].plot(fft_freq[HxF_eigval[i][0]], np.abs(HxF[i])[HxF_eigval[i][0]], 'r+', label = 'Punktmessung')
#     ac[1,i-3].plot(fft_freq[HxF_eigval[i][1]], np.abs(HxF[i])[HxF_eigval[i][1]], 'r+')
#     ac[1,i-3].plot(fft_freq[HxF_eigval1_combined], np.abs(HxF[i][HxF_eigval1_combined]), 'g*', label = 'Mittelwert')
#     ac[1,i-3].set_yscale('log')
#     ac[1,i-3].plot(fft_freq[0:1250], np.abs(HxF[i]))
#     ac[1,i-3].set_title('Punkt '+str(i+1))
#     ac[1,1].set_xlabel('Frequenz [Hz]')
# #plt.savefig(sub_dic+'\Übersicht_Amplitudenfrequenzgang.svg', format = 'svg', dpi = 2000)


'''PLOT 6: KONTROLLE DES ERMITTELTEN DÄMPFUNGSMAßES'''

# plot61, D_intplot = plt.subplots()
# plt.ylabel('Komplexe Amplitude')
# plt.xlabel('Frequenz in Hz')
# plt.title('Berechnung des Dämpfungsgrades')
# plt.yscale('log')
# plt.plot(fft_freq[5:1000] ,np.abs(HxF[p])[5:1000])
# plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new1l[p],4)*10000/(sweep_len*sampledt))], HxFmaxes1sqrt[p], 'r+')
# plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new1r[p],4)*10000/(sweep_len*sampledt))], HxFmaxes1sqrt[p], 'r+')
# plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new2l[p],4)*10000/(sweep_len*sampledt))], HxFmaxes2sqrt[p], 'r+')
# plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new2r[p],4)*10000/(sweep_len*sampledt))], HxFmaxes2sqrt[p], 'r+')

# plot62, Dplot = plt.subplots()
# plt.ylabel('Komplexe Amplitude')
# plt.xlabel('Frequenz in Hz')
# plt.title('Berechnung des Dämpfungsgrades')
# plt.yscale('log')
# plt.plot(fft_freq[5:1000] ,np.abs(HxF[p])[5:1000])
# plt.plot(fft_freq[df1_1[p]], np.abs(HxF1_1[p]), 'r+')
# plt.plot(fft_freq[df2_1[p]], np.abs(HxF2_1[p]), 'r+')
# plt.plot(fft_freq[df1_2[p]], np.abs(HxF1_2[p]), 'r+')
# plt.plot(fft_freq[df2_2[p]], np.abs(HxF2_2[p]), 'r+')



'''Savefiles für tabellarische Auswertung'''
# from tabulate import tabulate as tb
# import pandas as pd


# '''datasave = open(sub_dic+'\Messung Stats.txt', 'w')
# datasave.write('Erste Mode bei '+str(HxF_eigfreq1_combined)+'Hz (max Abweichung '+str((np.max(HxF_eigval1)-np.min(HxF_eigval1))/5)+'Hz)\n'+'Modale Daempfung '+str(round(D1_avg,5))+' \n'+'Zweite Mode bei '+str(+HxFeigfreqcombined_2)+'Hz (max Abweichung '+str((np.max(HxF_eigval2)-np.min(HxF_eigval2[HxF_eigval2 != 0]))/5)+'Hz)\n'+'Modale Daempfung '+str(round(D2_avg,5))+'\n'+'Messauflösung '+str(resolution)+'Hz')
# datasave.close()'''
# table = [['Messpunkt', 1,2,3,4,5,6], 
#           ['-','-','-','-','-','-','-'],
#           ['Erste Eigenfrequenz',HxF_eigfreq[0][0],HxF_eigfreq[1][0],HxF_eigfreq[2][0],HxF_eigfreq[3][0],HxF_eigfreq[4][0],HxF_eigfreq[5][0]],
#           ['Dämpfungsmaß 1te Mode', np.round(D1_int[0],5),np.round(D1_int[1],5),np.round(D1_int[2],5),np.round(D1_int[3],5),np.round(D1_int[4],5),np.round(D1_int[5],5)],
#           ['Auslenkung 1te Mode', np.round(phi_1[0],5), np.round(phi_1[1],5), np.round(phi_1[2],5), np.round(phi_1[3],5), np.round(phi_1[4],5), np.round(phi_1[5],5)],
#           ['-','-','-','-','-','-','-'],
#           ['Zweite Eigenfrequenz',HxF_eigfreq[0][1],HxF_eigfreq[1][1],HxF_eigfreq[2][1],HxF_eigfreq[3][1],HxF_eigfreq[4][1],HxF_eigfreq[5][1]],
#           ['Dämpfungsmaß 2te Mode', np.round(D2_int[0],5),np.round(D2_int[1],5),np.round(D2_int[2],5),np.round(D2_int[3],5),np.round(D2_int[4],5),np.round(D2_int[5],5)],
#           ['Auslenkung 2te Mode', np.round(phi_2[0],5), np.round(phi_2[1],5), np.round(phi_2[2],5), np.round(phi_2[3],5), np.round(phi_2[4],5), np.round(phi_2[5],5)],
#         ]
# print(tb(table))
# print('Erste Mode bei '+str(HxF_eigfreq1_combined)+'Hz\n'+'Modale Daempfung '+str(round(D1_avg,5))+' \n'+'Zweite Mode bei '+str(+HxFeigfreqcombined_2)+'Hz\n'+'Modale Daempfung '+str(round(D2_avg,5))) 

# dataframe = pd.DataFrame(table)
# dataframe.to_excel(sub_dic+'/auswertung_'+str(int(set_amp*100))+'V.xlsx', index=False, header=False)


# plot7, hxfplot = plt.subplots()
# plt.xlabel('Frequenz in Hz')
# plt.ylabel('Amplitude in mm/100V')
# if choosesensor == 'W':
#     plt.plot(fft_freq[0:1250], np.abs(HxF[0]), color='red')
#     plt.plot(fft_freq[0:1250], np.abs(HxF[1]), color='red')
#     plt.plot(fft_freq[0:1250], np.abs(HxF[2]), color='red')
#     plt.plot(fft_freq[0:1250], np.abs(HxF[3]), color='green')
#     plt.plot(fft_freq[0:1250], np.abs(HxF[4]), color='green')
#     plt.plot(fft_freq[0:1250], np.abs(HxF[5]), color='green')
# else:
#     plt.plot(fft_freq[0:1250], np.abs(HxF_abs_filt[0]), color='red')
#     plt.plot(fft_freq[0:1250], np.abs(HxF_abs_filt[1]), color='red')
#     plt.plot(fft_freq[0:1250], np.abs(HxF_abs_filt[2]), color='red')
#     plt.plot(fft_freq[0:1250], np.abs(HxF_abs_filt[3]), color='green')
#     plt.plot(fft_freq[0:1250], np.abs(HxF_abs_filt[4]), color='green')
#     plt.plot(fft_freq[0:1250], np.abs(HxF_abs_filt[5]), color='green')
# plt.yscale('log')
# plot7.tight_layout
# plot7.set_figwidth(10)
# plot7.set_figheight(6)
# # plt.savefig(sub_dic+'\Übersicht_Amplitudenfrequenzgang.jpeg', format = 'jpeg', dpi = 150)


