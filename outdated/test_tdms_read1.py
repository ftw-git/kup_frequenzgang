# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 12:16:30 2017

@author: Tobias
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal as si
import tkinter as tk
#from fft_t import fft_fit

from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from nptdms import TdmsFile
from tkinter import filedialog

plt.close('all')

#%% main
root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename()
##
#file = open("D:/Messdaten/Messung_Drossel_KuP/20170815/Msg11_P_30Hz_v_ohne.tdms", 'r+b')

tdms_file = TdmsFile(file_path)

channels = tdms_file.group_channels(tdms_file.groups()[0])

c = len(channels)

# Iterate over all items in the properties dictionary and print them
for chan in channels:
    for name, value in chan.properties.items():
         print("{0}: {1}".format(name, value))

T = np.array([1.,1.,1.,1.])
#for i in range(1,c):     
#    T[i]=channels[i].property('wf_increment')
#    print(T)
    
f = 1/T
N = int(channels[1].data.size)
#time = channels[1].time_track()

plt.figure()
plt.xlabel('Time'), plt.ylabel('Amplitude'), plt.grid()
for i in range(0,c):
    plt.plot(time,channels[i].data, label = channels[i].channel)
plt.legend(loc = 'best')    

freq = np.ones([c,int(N/2)])
SIGNAL = np.ones([c,int(N/2)])+1j*np.ones([c,int(N/2)])
H = np.ones([c,int(N)])

for i in range(0,c):
#    if i==2:
#        s = 0.1021/0.0997
#    elif i==3:
#        s = 0.0997/0.1021 
#    else:
#        s = 1
    (freq[i],SIGNAL[i]) = fft_fit(time,channels[i].data,channels[i],N)
    plt.figure()
    plt.plot(freq[i],abs(SIGNAL[i]),'b-');
    plt.xlabel('Frequenz [Hz] (Auflösung: %3.3f Hz)' %(freq[i][8]-freq[i][7]));
    plt.ylabel('Absolutwert');
    plt.title('Amplitudenfrequenzspektrum ' + channels[i].channel,fontsize=18);
    plt.axis([0, 500, 0, np.max(np.abs(SIGNAL[i]))])
    plt.grid()

   
H1 = SIGNAL[2]/SIGNAL[0]
#H2 = SIGNAL[3]/SIGNAL[1]

plt.figure()
#plt.plot(freq[3],abs(H1))
plt.semilogx(freq[2],10*np.log10(abs(H1)))
plt.grid()

#plt.figure()
#plt.semilogx(freq[3],10*np.log10(abs(H2)))
#plt.grid()

#------------------------------------------------------------------------------
# Waterfall Diagramm
#------------------------------------------------------------------------------
r =int(np.floor(time[-1])-5)
d = int(np.floor(3/T[1])) - int(np.floor(2/T[1]))
z = 300
ri = int(N/z)
h = np.ones([ri,int(d/2)])*0
fn = np.ones(ri)*0

for l in range(0,int(ri)-1):
    print('l = %02i'%l)
    (fq, h[l,:]) = fft_fit(time[int(np.floor(2/T[1]))+z*l:int(np.floor(3/T[1]))+z*l], channels[2].data[int(np.floor(2/T[1]))+z*l:int(np.floor(3/T[1]))+z*l],channels[2],d)
    ki = find(np.diff(channels[3].data[int(np.floor(2/T[1]))+z*l:int(np.floor(3/T[1]))+z*l])>0.05)
    fn[l] = np.mean(1/np.diff(time[ki+int(np.floor(2/T[1]))+z*l])[find(1/np.diff(time[ki+int(np.floor(2/T[1]))+z*l])<100)])
#    print(fn[l])

FQ,FN = np.meshgrid(fq[0:1000],fn)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot a basic wireframe.
ax.plot_wireframe(FQ, FN, h, rstride=00, cstride=10)
plt.show()

n = int(find(fn==np.nanmax(fn)))

fn_neu = np.linspace(fn[0],np.nanmax(fn),n)

#hi = sp.interpolate.griddata((fq, fq1), h, (FQ, FQ1), method='cubic')
f2 = plt.figure()
CS = plt.contour(fq[0:150],fn_neu, abs(h[0:n,0:150]),15,linewidths=0, colors='k')
CS = plt.contourf(fq[0:150], fn_neu, abs(h[0:n,0:150]), 15)
plt.colorbar(), #plt.xlabel(r'x/c'), plt.ylabel(r'y/t')
plt.title('waterfall', fontsize=15)
plt.tight_layout()
plt.show()  

