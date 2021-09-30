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

plt.close('all')

#%% Read data

# Asks for tdms-directory
print('Choose .tdms-Directory.\n')
root = tk.Tk()
root.withdraw()
sub_dic = filedialog.askdirectory(title='Choose .tdms-Directory.') #Öffnet Dialogfenster zum auswählen der Datei
print('\'{}\' was chosen.'.format(os.path.split(sub_dic)[-1]))

# For debugging
# sub_dic = 'C:/Users/FWass/TUB_KUP/kup_frequenzgang/measures/kleineAluschaufel_DS_sweep'

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

# Read in channel names
channel_names = [ data_import[0][idx].name for idx in range(len(data_import[0])) ]

# Sample time increment from channel properties, rounded to 7th decimal
sampledt = round(channel_prop['wf_increment'], 7)

# Constants of sweep
sweep_len = 20000
sample_freq = int(1/sampledt)
resolution = sample_freq/sweep_len
gain = 100
sensitivity = 0

# Number of tdms-files to be evaluated
num_of_files = int(len(data_import))
# num_of_files = 15

# Find minimal number of periods contained in every .tdms-file
data_size = []
for file in range(0, num_of_files):    
    data_size.append( data_import[file][0].data.size )
# Integer number of periods times frequency results in common signal length
data_size = int(min(data_size)/sample_freq)*sample_freq

# Rescale data to data_size
for file in range(0, num_of_files):    
    data_import[file] = [ data_import[file][tdmschannel][0:data_size] for tdmschannel in range(0, np.shape(data_import[0])[0]) ]

#%% Anregungssignal verarbeiten
# Anregung F: Kanal 0

# Initialize stimulation F
F = []
# Indices of peaks of stimulation
F_peaks = []
# Hessian of stimulation
d2Fdx2 = []
# Anregung auf eine Sweeplänge skaliert
F_short = []
# FFT of stimulation
F_fft = []

# num_of_files = 3

# Loop over all files
for file in range(0, num_of_files):
    # Get data from data_import
    # F.append( data_import[file][1]*gain/100 )
    F.append( data_import[file][channel_names.index('Input')]*gain/100 )
    # Calculate hessian
    d2Fdx2.append( np.gradient(np.gradient(F[file])) )
    
    # New method: Finds peaks through binary mask and filter
    data = abs(d2Fdx2[file])
    # Threshhold for peak identification
    peak_threshold= 0.6*max(data)
    # Binary mask for threshhold
    data_masked = data > peak_threshold
    # Transform bool to int
    data_masked = data_masked.astype(int)
    # Returns indices of peaks
    peak_idx = list( filter(lambda x: data_masked[x] == 1, range(len(data_masked))) )
    # Append to list of peak indices for each file
    F_peaks.append(peak_idx)
        
    # Check if pseudo peak
    for peak_num, peak in enumerate(F_peaks[file]):        
        # Tolerance in which other potential peaks are searched and compared
        x_tol = 20
        # Minimal distance between two peaks
        min_dist = 0.5*sweep_len
        # Check if any higher peak in tolerated x direction
        if any( abs(d2Fdx2[file][peak - x_tol:peak + x_tol]) > abs(d2Fdx2[file][peak]) ):
            # Removes pseudo peaks
            F_peaks[file].remove(peak)
        # Check if minimal distance between peaks is violated
        if peak_num < len(F_peaks[file])-1:
            if np.diff(F_peaks[file])[peak_num] < min_dist:
                # Removes pseudo peaks
                F_peaks[file].remove(F_peaks[file][peak_num + 1])
    
    # # Plot peaks and original signal, Debugging
    # fig = plt.figure()
    # plt.plot(abs(d2Fdx2[file]))
    # plt.plot(F_peaks[file], abs(d2Fdx2[file])[F_peaks[file]], '+r')
               
    # Reshape F to rows -> one sweep length and average over it
    # Shortening of stimulation signal
    F_short.append( np.mean( F[file].reshape( ( int(data_size/sweep_len), sweep_len ) ), 0 ) )
    
    # FFT of stimulation for frequencies 0Hz - 250Hz
    F_fft.append( np.fft.fft(F_short[file], norm = 'ortho')[0:1250]/(np.pi) )

# Calculate amplitude and decide on sensor to be used
# Get amplitude from max of stimulation
set_amp = round( max( F[0] ), 1 )

# =============================================================================
# # Decides on sensor to evaluate depending on amplitude
# # 'L' - Laser-Sensor
# if set_amp > 1.99:
#     choosesensor = 'L'
#     sensitivity = 0.4 
#     channel_id = channel_names.index('Laser-Triangulation')
# # 'W' - Wirbelstrom-Sensor
# else: 
#     choosesensor = 'W' 
#     sensitivity = -8
#     channel_id = channel_names.index('Wirbelstrom')
# =============================================================================
# Uncomment the above when 'Wirbelstrom' is used, actually never the case
choosesensor = 'L'
channel_id = channel_names.index('Laser-Triangulation')
sensitivity = 0.4

# Removing linear trend from measurements and separating data into SS and DS
# Iinitialize
X = np.ones([num_of_files, data_size])
X_short = np.ones([num_of_files, sweep_len])

# Loop over all files
for file in range(0, num_of_files):
    # Separates datasets into Saugseite and Druckseite due to direction of measurements
    if file < int(num_of_files/2):
        X[file] =  sp.detrend(data_import[file][channel_id].data)/(-1*sensitivity)
    else:
        X[file] = sp.detrend(data_import[file][channel_id].data)/(sensitivity)
        
    # Shortening measurements in analogy to stimulation by reshaping + averaging
    # X_short[file] = X[file][F_peaks[file]:(F_peaks[file]+sweep_len)]
    X_short[file] = np.mean( X[file].reshape( ( int(data_size/sweep_len), sweep_len ) ), 0 )
    
# Butterworth-Filter for measurements
N = 30 # Order of filter
Wn = 250 # Critical frequency
# output-type: sos, second-order-sections for general filter purposes
# fs = sampling frequency of digital system
butter_X = sp.butter(N, Wn, output = 'sos', btype = 'low', fs = sample_freq)
for file in range(0, num_of_files):
        X_short[file] = sp.sosfiltfilt(butter_X, X_short[file])

#%% Koordinaten-Transformation, graphically estimated ??? -> NX?

# =============================================================================
# ko_cor = [ 1/np.cos(np.radians(-1.146)), 
#           1/np.cos(np.radians(-16.7)), 
#           1/np.cos(np.radians(-25.1)), 
#           1/np.cos(np.radians(-1.146)), 
#           1/np.cos(np.radians(-6.84)), 
#           1/np.cos(np.radians(-5)) ]
# # Transforms over all files
# for i in range(0, num_of_files):
#     X_short[i] = X_short[i]*ko_cor[i]
# =============================================================================


#%% FFT of measurements and calculation of Transfer Function, Determination of Eigen-Frequencies

X_fft = [] 
HxF = [] 

HxF_abs_filt = [] 
HxF_eigval = [] 
HxF_eigval1 = [] 
HxF_eigval2 = []

butter_HxF = sp.butter(4, 65, output='sos', btype='low', fs=500)

HxF_eigfreq = [] #np.empty([num_of_files, data_size])#,dtype=list)

# Loop over all files -> Woher die ganzen statischen Variablen??
for file in range(0, num_of_files):
    # DFFT of shortened measurements
    X_fft.append( np.fft.rfft(X_short[file], norm = 'ortho')[0:1250]/(2*np.pi) )
    
    # Transfer function from output/input
    HxF.append( X_fft[file]/F_fft[file] )
    
    # Filter if laser-data is used
    if choosesensor == 'L':
        HxF_abs_filt.append( sp.sosfiltfilt(butter_HxF, np.abs(HxF[file])) )
        HxF_eigval.append( sp.find_peaks(HxF_abs_filt[file][30:1200], np.max(HxF_abs_filt[file][30:1200])*0.2, distance = 70)[0]+30 )
        HxF_eigval1.append( HxF_eigval[file][0] )
    # Else if Wirbelstrom, no filtering needed
    else:
        HxF_eigval.append( sp.find_peaks(np.abs(HxF[file][30:1200]), np.max(np.abs(HxF[file][30:1200]))*0.2, distance = 70)[0]+30 )
        HxF_eigval1.append( HxF_eigval[file][0] )
    # If more then one eigenvalue, store second eigenvalue too
    if len(HxF_eigval[file]) > 1:
        HxF_eigval2.append( HxF_eigval[file][1] )
    # Else: second eigenvalue is zero
    else:
        HxF_eigval2.append( 0 )

# Loop over all files
for file in range(0, num_of_files):
    # To find second eigenvalue if not found before
    if HxF_eigval2[file] == 0:
        # For laser-Sensor
        if choosesensor == 'L':
            HxF_eigval2[file] = sp.find_peaks(np.abs(HxF_abs_filt[file][int(np.median(HxF_eigval2)-50):int(np.median(HxF_eigval2)+50)]),
                                            np.max(np.abs(HxF_abs_filt[file][int(np.median(HxF_eigval2)-50):int(np.median(HxF_eigval2)+50)]))*0.2,
                                            distance = 50)[0][0]+(np.median(HxF_eigval2)-50)
            HxF_eigval[file] = np.append(int(HxF_eigval[file]), int(HxF_eigval2[file]))
            HxF_eigfreq.append( HxF_eigval[file]/(sweep_len*sampledt) )
        # For Wirbelstrom
        else:
            HxF_eigval2[file] = sp.find_peaks(np.abs(HxF[file][int(120-50):int(120+50)]),
                                np.max(np.abs(HxF[file][int(120-50):int(120+50)]))*0.2,
                                distance = 50)[0][0]+(120-50)
            HxF_eigval[file] = np.append(int(HxF_eigval[file]), int(HxF_eigval2[file]))
            HxF_eigfreq.append( HxF_eigval[file]/(sweep_len*sampledt) )
    else:
        HxF_eigval2[file] = HxF_eigval2[file]
        # HxF_eigfreq[i] = HxF_eigval[i]/(sweep_len*sampledt)
        HxF_eigfreq.append( HxF_eigval[file]/(sweep_len*sampledt) )

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


#%% Calculate multiple k of second eigenvalue relative to first
k = []

for file in range(0, num_of_files):
    if HxF_eigval2[file]:
        k.append(HxF_eigval2[file]/HxF_eigval1[file])
    else:
        pass

#%% Berechnen der Dämpfungsmaße der einzelnen Eigenfrequenzen für jeden Punkt

HxFmaxes1 = [] 
HxFmaxes2 = [] 
HxFmaxes1sqrt = [] 
HxFmaxes2sqrt = [] 

df1_1 = [] 
df2_1 = []
df1_2 = [] 
df2_2 = [] 
HxF1_1 = [] 
HxF2_1 = [] 
HxF1_2 = [] 
HxF2_2 = [] 
D1 = [] 
D2 = [] 


for file in range(0, num_of_files):
    HxFmaxes1.append( np.abs(HxF[file][HxF_eigval[file][0]]) )
    HxFmaxes1sqrt.append( HxFmaxes1[file]/np.sqrt(2) )
    if choosesensor == 'L':
        HxFmaxes2.append( np.abs(HxF_abs_filt[file][HxF_eigval[file][1]]) )
        HxFmaxes2sqrt.append( HxFmaxes2[file]/np.sqrt(2) )
    else:
        HxFmaxes2.append( np.abs(HxF[file][HxF_eigval[file][1]]) )
        HxFmaxes2sqrt.append( HxFmaxes2[file]/np.sqrt(2) )
        continue  
for file in range(0, num_of_files):
    df1_1.append( (np.abs(np.abs(HxF[file] [30:HxF_eigval[file] [0]])-(np.abs(HxFmaxes1[file] )/np.sqrt(2)))).argmin()+30 )
    df2_1.append(  (np.abs(np.abs(HxF[file] [HxF_eigval[file] [0]:HxF_eigval[file] [0]+30])-(np.abs(HxFmaxes1[file] )/np.sqrt(2)))).argmin()+HxF_eigval[file] [0] )
    HxF1_1.append( np.abs(HxF[file] [df1_1[file] ]) )
    HxF2_1.append( np.abs(HxF[file] [df2_1[file] ]) )
    
    if choosesensor == 'L':       
        
        df1_2.append( (np.abs(np.abs(HxF_abs_filt[file] [HxF_eigval[file] [1]-20:HxF_eigval[file] [1]])-(np.abs(HxFmaxes2[file] )/np.sqrt(2)))).argmin()+HxF_eigval[file] [1]-20 )
        df2_2.append( (np.abs(np.abs(HxF_abs_filt[file] [HxF_eigval[file] [1]:HxF_eigval[file] [1]+20])-(np.abs(HxFmaxes2[file] )/np.sqrt(2)))).argmin()+HxF_eigval[file] [1] )
        HxF1_2.append( np.abs(HxF_abs_filt[file] [df1_2[file] ]) )
        HxF2_2.append( np.abs(HxF_abs_filt[file] [df2_2[file] ]) )
        
    else:
        
        df1_2.append( (np.abs(np.abs(HxF[file] [HxF_eigval[file] [1]-20:HxF_eigval[file] [1]])-(np.abs(HxFmaxes2[file] )/np.sqrt(2)))).argmin()+HxF_eigval[file] [1]-20 )
        df2_2.append( (np.abs(np.abs(HxF[file] [HxF_eigval[file] [1]:HxF_eigval[file] [1]+20])-(np.abs(HxFmaxes2[file] )/np.sqrt(2)))).argmin()+HxF_eigval[file] [1] )
        HxF1_2.append( np.abs(HxF[file] [df1_2[file] ]) )
        HxF2_2.append( np.abs(HxF[file] [df2_2[file] ]) )
        
for file in range(0, num_of_files):
    D1.append( ((df2_1[file] /(sweep_len*sampledt)*2*np.pi - df1_1[file] /(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[file] [0])/(sweep_len*sampledt)*2*np.pi)) )

    if len(HxF_eigval[file] ) > 1:
        D2.append( ((df2_2[file] /(sweep_len*sampledt)*2*np.pi - df1_2[file] /(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[file] [1])/(sweep_len*sampledt)*2*np.pi)) )
    else:
        D2.append( 0 )


#%% Interpolation zur Verfeinerung der Dämpfungsmaße

x_new1l = []  
x_new1r = [] 
D1_int = [] 
x_new2l = [] 
x_new2r = [] 
D2_int = []

for file in range(0, num_of_files):
    
    x_new1l.append( np.interp(HxFmaxes1sqrt[file] , [np.abs(HxF[file] [df1_1[file] -2]), np.abs(HxF[file] [df1_1[file] ]), np.abs(HxF[file] [df1_1[file] +2])], [df1_1[file] -2, df1_1[file] , df1_1[file] +2]) )
    x_new1r.append( np.interp(HxFmaxes1sqrt[file] , [np.abs(HxF[file] [df2_1[file] +2]), np.abs(HxF[file] [df2_1[file] ]), np.abs(HxF[file] [df2_1[file] -2])], [df2_1[file] +2, df2_1[file] , df2_1[file] -2]) )
    D1_int.append( (x_new1r[file] /(sweep_len*sampledt)*2*np.pi-x_new1l[file] /(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[file] [0])/(sweep_len*sampledt)*2*np.pi) )
        
    if choosesensor == 'L':
        
        x_new2l.append( np.interp(HxFmaxes2sqrt[file] , [HxF_abs_filt[file] [df1_2[file] -2], np.abs(HxF_abs_filt[file] [df1_2[file] ]), np.abs(HxF_abs_filt[file] [df1_2[file] +2])], [df1_2[file] -2, df1_2[file] , df1_2[file] +2]) )
        x_new2r.append( np.interp(HxFmaxes2sqrt[file] , [np.abs(HxF[file] [df2_2[file] +2]), np.abs(HxF[file] [df2_2[file] ]), np.abs(HxF[file] [df2_2[file] -2])], [df2_2[file] +2, df2_2[file] , df2_2[file] -2]) )
        D2_int.append( (x_new2r[file] /(sweep_len*sampledt)*2*np.pi-x_new2l[file] /(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[file] [1])/(sweep_len*sampledt)*2*np.pi) )
        
    else:
        
        x_new2l.append( np.interp(HxFmaxes2sqrt[file] , [np.abs(HxF[file] [df1_2[file] -2]), np.abs(HxF[file] [df1_2[file] ]), np.abs(HxF[file] [df1_2[file] +2])], [df1_2[file] -2, df1_2[file] , df1_2[file] +2]) )
        x_new2r.append( np.interp(HxFmaxes2sqrt[file] , [np.abs(HxF[file] [df2_2[file] +2]), np.abs(HxF[file] [df2_2[file] ]), np.abs(HxF[file] [df2_2[file] -2])], [df2_2[file] +2, df2_2[file] , df2_2[file] -2]) )
        D2_int.append( (x_new2r[file] /(sweep_len*sampledt)*2*np.pi-x_new2l[file] /(sweep_len*sampledt)*2*np.pi)/((2*HxF_eigval[file] [1])/(sweep_len*sampledt)*2*np.pi) )
    

#%% Durchschnittliche Dämpfungsgrade berechnen
# Eigenvektor bilden (physikalisch r_1 und r_2 sowie regelungstechnisch phi_1 und phi_2)

HxFcomplexval1 = [] 
HxFcomplexval2 = [] 
HxFangle1 = [] 
HxFangle2 = [] 

for file in range(0, num_of_files):

    HxFcomplexval1.append( HxF[file] [int(HxF_eigval1[file] )] )
    HxFcomplexval2.append( HxF[file] [int(HxF_eigval2[file] )] )
    HxFangle1.append( np.angle(HxFcomplexval1[file] , deg = True) )
    HxFangle2.append( np.angle(HxFcomplexval2[file] , deg = True) )

mode_1_dir = np.sign(np.array( HxFcomplexval1 ).imag)
mode_2_dir = np.sign(np.array( HxFcomplexval2 ).imag)

phi_1 = [] 
phi_2 = [] 

for file in range(0, num_of_files):
    phi_1.append( HxFmaxes1[file] *mode_1_dir[file]  )
    phi_2.append( HxFmaxes2[file] *mode_2_dir[file]  )
    
D1_avg = np.average(D1_int)

if np.sum(D2) == 0:
    D2_avg = 0
else:
    D2_avg = np.average(D2_int[D2_int != 0])
    
r_1 = np.empty(1, dtype=float)
r_1[0] = np.sqrt(np.abs((HxFcomplexval1[0])*(2*np.pi*HxF_eigfreq[0][0])**2))*(D1_int[0]**4+4*D1_int[0]**2)**(1/4)

for file in np.arange(1,num_of_files,1):
    r_1 = np.append(r_1, np.abs((HxFcomplexval1[file] )*(2*np.pi*HxF_eigfreq[file] [0])**2)/r_1[0]*np.sqrt(D1_int[file] **4+4*D1_int[file] **2))

r_2 = np.empty(1, dtype=float)
r_2[0] = np.sqrt(np.abs((HxFcomplexval2[0])*(2*np.pi*HxF_eigfreq[0][1])**2))*(D2_int[0]**4+4*D2_int[0]**2)**(1/4)

for file in np.arange(1,num_of_files,1):
    r_2 = np.append(r_2, np.abs((HxFcomplexval2[file] )*(2*np.pi*HxF_eigfreq[file] [1])**2)/r_1[0]*np.sqrt(D2_int[file] **4+4*D2_int[file] **2))

r_1 = r_1*mode_1_dir
r_2 = r_2*mode_2_dir
r_1_norm = r_1/np.max(np.abs(r_1))
r_2_norm = r_2/np.max(np.abs(r_2))

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
        p = int(input('Which point out of the {} points should be evaluated?\n'.format(num_of_files)))
        if 0 < p-1 < num_of_files :
            break
    except (ValueError, TypeError):
        print('You\'ve entered an invalid input.')
        pass

#%% Plot 0: Excitation and measurements in point P -> Check!

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

#%% Plot 1: Amplitude- and Phaseresponse of Transfer-Function -> Check!

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
# #ax_1[0].plot(fft_freq[10:1000], HxF_abs_filt[5][10:1000])
# # Save plot in measurement directory
# # plt.savefig( sub_dic + str(p+1) + 'P'  + str(set_amp*100) + 'V.png', dpi = 2000)

#%% Plot 2: First and second mode plotted on blade -> TO DO !!!
# ??? Punkte müssen neu vermessen werden

# Green dot: selected null point

# Punkte in [mm] ausgehend von unterer Einspannung
meas_pos = np.genfromtxt( 
                        os.getcwd() + '/Daten_Schaufel_Motoren_Saugseite.txt',
                        dtype=float,
                        delimiter='\t'
                        )
# Read positions and correct offset
px = 82.1 - meas_pos[:, 0]
py = -5 + meas_pos[:, 1]
p_ang = meas_pos[:, 2]


plot2, ay = plt.subplots(2,1)
imgschaufel = plt.imread('Schaufel.JPG')
plot2.set_figwidth(12)
plot2.set_figheight(8)
ay[0].imshow(imgschaufel, extent = [82.1,-69.9, -15.8, int(152*451/1692)-15.8])
ay[0].plot(0,0,'g+')
ay[0].plot(px,py, 'r*')
ay[0].set_title('Mode 1')
ay[1].set_xlabel('x-Koordinate in mm')
ay[1].set_ylabel('y-Koordinate in mm')
plot2.tight_layout()
plotscale = 10
for i in range(0, num_of_files):
    ay[0].arrow(px[i], py[i], 0, 1*plotscale*phi_1[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')

ay[1].imshow(imgschaufel, extent = [82.1,-69.9, -15.8, int(152*451/1692)-15.8])
ay[1].plot(0,0,'g+')
ay[1].plot(px,py, 'r*')
ay[1].set_title('Mode 2')
for i in range(0, num_of_files): 
    ay[1].arrow(px[i], py[i], 0, -1*plotscale*2*phi_2[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')
#plt.savefig(sub_dic+'\Modes.png', dpi = 2000)

#%% Plot 3: Comparison Bode-Diagrams all Measurement-Points -> Check!

# plot3, az = plt.subplots()
# plot3.set_figwidth(12)
# plot3.set_figheight(8)
# az.set_ylabel('Komplexe Amplitude')
# az.set_xlabel('Frequenz in Hz')
# az.set_yscale("log")
# plot3.tight_layout()
# az.grid(which = 'major', axis = 'y', linewidth = 0.5, linestyle='--')
# for i in range(0, num_of_files):
#     az.plot(fft_freq[0:1250], np.abs(HxF[i]), label = 'Punkt '+str(i+1))
# az.legend()
# #plt.savefig(sub_dic+'\Amplitudenfrequenzgänge.svg', format = 'svg', dpi = 2000)

#%% Plot 4: Excitation, Deflection and Transfer Function of one measurement point -> Check!

# plot_4, ax_4 = plt.subplots(3,1)
# plot_4.set_figwidth(12)
# plot_4.set_figheight(8)
# ax_4[0].set_ylabel('Anregung in V')
# ax_4[1].set_ylabel('Messung in V')
# ax_4[2].set_ylabel('Komplexe Amplitude')
# ax_4[0].set_xlabel('Zeit in s')
# ax_4[1].set_xlabel('Zeit in s')
# ax_4[2].set_xlabel('Frequenz in Hz')
# ax_4[2].set_yscale('log')
# plot_4.tight_layout()
# # ax_4[0].plot(np.arange(0,50000,1)/10000,F_short[p])
# # ax_4[1].plot(np.arange(0,50000,1)/10000,X_short[p])
# # ax_4[2].plot(fft_freq[:1250],np.abs(HxF[p][:1250]))
# ax_4[0].plot(np.arange(0,sweep_len,1)/sweep_len,F_short[p])
# ax_4[1].plot(np.arange(0,sweep_len,1)/sweep_len,X_short[p])
# ax_4[2].plot(fft_freq[:1250],np.abs(HxF[p]))
# #plt.savefig(sub_dic+'\Punkt '+str(p+1)+' Anregung_Messung_Übertragungsfunktion_1.svg', format = 'svg', dpi = 2000)'''

#%% Plot 5: Single comparison of all Transfer Function -> Check!

# plot5, ac = plt.subplots(2,3, sharey=True)
# plot5.set_figwidth(12)
# plot5.set_figheight(8)
# plot5.tight_layout()
# # Loop over Saug- and Druckseite! Range needs to be updated 
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

#%% PLOT 6: KONTROLLE DES ERMITTELTEN DÄMPFUNGSMAßES -> Check!

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

#%% Plot 7 -> Check!

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


#%% Savefiles für tabellarische Auswertung''' -> TO DO !!!

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


