# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:18:27 2020

@author: Philipp Tuntsch
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp

'''Messdaten importieren'''
P1 = [
      np.genfromtxt('HammerTest_Component1_1_+X_1_1.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_2_1.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_3_1.csv', delimiter = ',', skip_header=2)
      ]

P2 = [
      np.genfromtxt('HammerTest_Component1_1_+X_1_2.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_2_2.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_3_2.csv', delimiter = ',', skip_header=2)
      ]

P3 = [
      np.genfromtxt('HammerTest_Component1_1_+X_1_3.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_2_3.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_3_3.csv', delimiter = ',', skip_header=2)
      ]

P4 = [
      np.genfromtxt('HammerTest_Component1_1_+X_1_4.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_2_4.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_3_4.csv', delimiter = ',', skip_header=2)
      ]

P5 = [
      np.genfromtxt('HammerTest_Component1_1_+X_1_5.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_2_5.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_3_5.csv', delimiter = ',', skip_header=2)
      ]

P6 = [
      np.genfromtxt('HammerTest_Component1_1_+X_1_6.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_2_6.csv', delimiter = ',', skip_header=2),
      np.genfromtxt('HammerTest_Component1_1_+X_3_6.csv', delimiter = ',', skip_header=2)
      ]

C = np.genfromtxt("HammerTest_Coherence.csv", delimiter = ',', skip_header=2)
FRF = np.genfromtxt('HammerTest_FRF.csv', delimiter = ',', skip_header=2)

####
frf1 = FRF[:,0]+1j*FRF[:,1]
frf2 = FRF[:,2]+1J*FRF[:,3]
frf3 = FRF[:,4]+1J*FRF[:,5]
frf4 = FRF[:,6]+1J*FRF[:,7]
frf5 = FRF[:,8]+1J*FRF[:,9]
frf6 = FRF[:,10]+1J*FRF[:,11]
####

ko_cor = [1/np.cos(np.radians(-1.146)),1/np.cos(np.radians(-16.7)),1/np.cos(np.radians(-25.1)),1/np.cos(np.radians(-1.146)),1/np.cos(np.radians(-6.84)),1/np.cos(np.radians(-5))]

F1 = [P1[0][:,0], P1[1][:,0], P1[2][:,0]]
x1 = [sp.detrend(P1[0][:,1])/1000*ko_cor[0], sp.detrend(P1[1][:,1])/1000*ko_cor[0], sp.detrend(P1[2][:,1])/1000*ko_cor[0]]

F2 = [P2[0][:,0], P2[1][:,0], P2[2][:,0]]
x2 = [sp.detrend(P2[0][:,1])/1000*ko_cor[1], sp.detrend(P2[1][:,1])/1000*ko_cor[1], sp.detrend(P2[2][:,1])/1000*ko_cor[1]]

F3 = [P3[0][:,0], P3[1][:,0], P3[2][:,0]]
x3 = [sp.detrend(P3[0][:,1])/1000*ko_cor[2], sp.detrend(P3[1][:,1])/1000*ko_cor[2], sp.detrend(P3[2][:,1])/1000*ko_cor[2]]

F4 = [P4[0][:,0], P4[1][:,0], P4[2][:,0]]
x4 = [sp.detrend(P4[0][:,1])/1000*ko_cor[3], sp.detrend(P4[1][:,1])/1000*ko_cor[3], sp.detrend(P4[2][:,1])/1000*ko_cor[3]]

F5 = [P5[0][:,0], P5[1][:,0], P5[2][:,0]]
x5 = [sp.detrend(P5[0][:,1])/1000*ko_cor[4], sp.detrend(P5[1][:,1])/1000*ko_cor[4], sp.detrend(P5[2][:,1])/1000*ko_cor[4]]

F6 = [P6[0][:,0], P6[1][:,0], P6[2][:,0]]
x6 = [sp.detrend(P6[0][:,1])/1000*ko_cor[5], sp.detrend(P6[1][:,1])/1000*ko_cor[5], sp.detrend(P6[2][:,1])/1000*ko_cor[5]]

'''FFT von Wirbelstromsensor und Hammerschlag'''
F1fft = np.transpose(np.array([np.fft.fft(F1[0], norm = 'ortho')/np.pi, np.fft.fft(F1[1], norm = 'ortho')/np.pi, np.fft.fft(F1[2], norm = 'ortho')/np.pi]))
x1fft = np.transpose(np.array([np.fft.fft(x1[0], norm = 'ortho')/np.pi, np.fft.fft(x1[1], norm = 'ortho')/np.pi, np.fft.fft(x1[2], norm = 'ortho')/np.pi]))
HxF_1 = x1fft/F1fft
F2fft = np.transpose(np.array([np.fft.fft(F2[0], norm = 'ortho')/np.pi, np.fft.fft(F2[1], norm = 'ortho')/np.pi, np.fft.fft(F2[2], norm = 'ortho')/np.pi]))
x2fft = np.transpose(np.array([np.fft.fft(x2[0], norm = 'ortho')/np.pi, np.fft.fft(x2[1], norm = 'ortho')/np.pi, np.fft.fft(x2[2], norm = 'ortho')/np.pi]))
HxF_2 = x2fft/F2fft
F3fft = np.transpose(np.array([np.fft.fft(F3[0], norm = 'ortho')/np.pi, np.fft.fft(F3[1], norm = 'ortho')/np.pi, np.fft.fft(F3[2], norm = 'ortho')/np.pi]))
x3fft = np.transpose(np.array([np.fft.fft(x3[0], norm = 'ortho')/np.pi, np.fft.fft(x3[1], norm = 'ortho')/np.pi, np.fft.fft(x3[2], norm = 'ortho')/np.pi]))
HxF_3 = x3fft/F3fft
F4fft = np.transpose(np.array([np.fft.fft(F4[0], norm = 'ortho')/np.pi, np.fft.fft(F4[1], norm = 'ortho')/np.pi, np.fft.fft(F4[2], norm = 'ortho')/np.pi]))
x4fft = np.transpose(np.array([np.fft.fft(x4[0], norm = 'ortho')/np.pi, np.fft.fft(x4[1], norm = 'ortho')/np.pi, np.fft.fft(x4[2], norm = 'ortho')/np.pi]))
HxF_4 = x4fft/F4fft
F5fft = np.transpose(np.array([np.fft.fft(F5[0], norm = 'ortho')/np.pi, np.fft.fft(F5[1], norm = 'ortho')/np.pi, np.fft.fft(F5[2], norm = 'ortho')/np.pi]))
x5fft = np.transpose(np.array([np.fft.fft(x5[0], norm = 'ortho')/np.pi, np.fft.fft(x5[1], norm = 'ortho')/np.pi, np.fft.fft(x5[2], norm = 'ortho')/np.pi]))
HxF_5 = x5fft/F5fft
F6fft = np.transpose(np.array([np.fft.fft(F6[0], norm = 'ortho')/np.pi, np.fft.fft(F6[1], norm = 'ortho')/np.pi, np.fft.fft(F6[2], norm = 'ortho')/np.pi]))
x6fft = np.transpose(np.array([np.fft.fft(x6[0], norm = 'ortho')/np.pi, np.fft.fft(x6[1], norm = 'ortho')/np.pi, np.fft.fft(x6[2], norm = 'ortho')/np.pi]))
HxF_6 = x6fft/F6fft
fftfreq = np.fft.fftfreq(len(HxF_1), d = 1/1024)

'''Bilden des Mittels aus den 3 durchgeführten Schlägen'''
HxF = [np.mean(HxF_1[0:1000], axis=1), np.mean(HxF_2[0:1000], axis=1), np.mean(HxF_3[0:1000],axis=1), np.mean(HxF_4[0:1000],axis=1), np.mean(HxF_5[0:1000],axis=1), np.mean(HxF_6[0:1000],axis=1)]

'''Ermitteln der Eigenfrequenzen'''
sweeplen = 4096
sampledt = 1/1024
count = np.arange(0,6,1)
HxFeigvalues = np.empty(6, dtype=list)
HxFeigfrequencies = np.empty(6,dtype=list)
HxFeigvalue1 = np.empty(6, dtype=float)
HxFeigvalue2 = np.empty(6, dtype=float)
for i in count:
    HxFeigvalues[i] = sp.find_peaks(np.abs(HxF[i][30:1000]), np.max(np.abs(HxF[i][30:1000]))*0.2, distance = 70)[0]+30
    HxFeigvalue1[i] = HxFeigvalues[i][0]
    if len(HxFeigvalues[i]) > 1:
        HxFeigvalue2[i] = HxFeigvalues[i][1]
    else:
        HxFeigvalue2[i] = 0
for i in count:
    if HxFeigvalue2[i] == 0:
        HxFeigvalue2[i] = sp.find_peaks(np.abs(HxF[i][int(HxFeigvalue2[5]-50):int(HxFeigvalue2[5]+50)]), np.max(np.abs(HxF[i][int(HxFeigvalue2[5]-50):int(HxFeigvalue2[5]+50)]))*0.5, distance = 50)[0][0]+(HxFeigvalue2[5]-50)
        HxFeigvalues[i] = np.append(int(HxFeigvalues[i]), int(HxFeigvalue2[i]))
        HxFeigfrequencies[i] = HxFeigvalues[i]/(sweeplen*sampledt)
    else:
        HxFeigvalue2[i] = HxFeigvalue2[i]
        HxFeigfrequencies[i] = HxFeigvalues[i]/(sweeplen*sampledt)
HxFeigvaluecombined_1 = int(np.round(np.average(HxFeigvalue1[HxFeigvalue1 != 0]),0))
HxFeigfreqcombined_1 = np.round(HxFeigvaluecombined_1/(sweeplen*sampledt),2)
if np.sum(HxFeigvalue2) == 0:
    HxFeigvaluecombined_2 = 0
    HxFeigfreqcombined_2 = 0
else:
    HxFeigvaluecombined_2 = int(np.round(np.average(HxFeigvalue2[HxFeigvalue2 != 0]),0))
    HxFeigfreqcombined_2 = np.round(HxFeigvaluecombined_2/(sweeplen*sampledt), 2)

'''Ermitteln der Dämpfungsmaße'''    
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
    if len(HxFeigvalues[i]) > 1:
        HxFmaxes2[i] = np.abs(HxF[i][HxFeigvalues[i][1]])
        HxFmaxes2sqrt[i] = HxFmaxes2[i]/np.sqrt(2)
    else:
        continue  
for i in count:
    df1_1[i] = (np.abs(np.abs(HxF[i][30:HxFeigvalues[i][0]])-(np.abs(HxFmaxes1[i])/np.sqrt(2)))).argmin()+30
    df2_1[i] = (np.abs(np.abs(HxF[i][HxFeigvalues[i][0]:HxFeigvalues[i][0]+30])-(np.abs(HxFmaxes1[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][0]
    HxF1_1[i] = np.abs(HxF[i][df1_1[i]])
    HxF2_1[i] = np.abs(HxF[i][df2_1[i]])
    if len(HxFeigvalues[i]) > 1:
        df1_2[i] = (np.abs(np.abs(HxF[i][HxFeigvalues[i][1]-100:HxFeigvalues[i][1]])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][1]-100
        df2_2[i] = (np.abs(np.abs(HxF[i][HxFeigvalues[i][1]:HxFeigvalues[i][1]+25])-(np.abs(HxFmaxes2[i])/np.sqrt(2)))).argmin()+HxFeigvalues[i][1]   
        HxF1_2[i] = np.abs(HxF[i][df1_2[i]])
        HxF2_2[i] = np.abs(HxF[i][df2_2[i]])
    else:
        df1_2[i] = 0
        df2_2[i] = 0
for i in count:
    D1[i] = ((df2_1[i] - df1_1[i])/(2*HxFeigvalues[i][0]))
    if len(HxFeigvalues[i]) > 1:
        D2[i] = ((df2_2[i] - df1_2[i])/(2*HxFeigvalues[i][1]))
    else:
        D2[i] = 0

'''Verfeinerung der Dämpfungsmaße durch Interpolation'''
x_new1l = np.empty(6, dtype=float)      
x_new1r = np.empty(6, dtype=float)      
D1_int = np.empty(6, dtype=float)
x_new2l = np.empty(6, dtype=float) 
x_new2r = np.empty(6, dtype=float)      
D2_int = np.empty(6, dtype=float)
for i in count:
    x_new1l[i] = np.interp(HxFmaxes1sqrt[i], [np.abs(HxF[i][df1_1[i]-1]), np.abs(HxF[i][df1_1[i]]), np.abs(HxF[i][df1_1[i]+1])], [df1_1[i]-1, df1_1[i], df1_1[i]+1])
    x_new1r[i] = np.interp(HxFmaxes1sqrt[i], [np.abs(HxF[i][df2_1[i]+1]), np.abs(HxF[i][df2_1[i]]), np.abs(HxF[i][df2_1[i]-1])], [df2_1[i]+1, df2_1[i], df2_1[i]-1])
    x_new2l[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df1_2[i]-1]), np.abs(HxF[i][df1_2[i]]), np.abs(HxF[i][df1_2[i]+1])], [df1_2[i]-1, df1_2[i], df1_2[i]+1])
    x_new2r[i] = np.interp(HxFmaxes2sqrt[i], [np.abs(HxF[i][df2_2[i]+1]), np.abs(HxF[i][df2_2[i]]), np.abs(HxF[i][df2_2[i]-1])], [df2_2[i]+1, df2_2[i], df2_2[i]-1])
    D1_int[i] = (x_new1r[i]-x_new1l[i])/(2*HxFeigvalues[i][0])
    D2_int[i] = (x_new2r[i]-x_new2l[i])/(2*HxFeigvalues[i][1])

'''Bilden der Auslenkungsvektoren'''
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
plotdir = np.array([1,1,1,-1,-1,-1])
r_1 = r_1*mode_1_dir*plotdir
r_2 = r_2*mode_2_dir*plotdir
r_1_norm = r_1/np.max(np.abs(r_1))
r_2_norm = r_2/np.max(np.abs(r_2))


'''PLOTS'''

'''KOHÄRENZ'''
'''plot001, aa = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3,1]})
aa[0].set_yscale('log')
aa[0].plot(np.arange(1.25,250,0.25), np.abs(HxF_4[:,0][5:1000]))
aa[0].plot(np.arange(1.25,250,0.25), np.abs(HxF_4[:,1][5:1000]))
aa[0].plot(np.arange(1.25,250,0.25), np.abs(HxF_4[:,2][5:1000]))
aa[0].set_ylabel('Amplitude in mm/N')
aa[1].plot(np.arange(1.25,250,0.25), C[:,0][5:1000])
aa[1].set_ylabel('Kohärenz')
aa[1].set_xlabel('Frequenz in Hz')
#'''

'''Punkte in [mm] ausgehend von unterer Einspannung'''
p1 = [28, 14.2]
p2 = [60.5, 6.75]
p3 = [75, 0.85]
p4 = [75, -6.15]
p5 = [49, -4.2]
p6 = [26, -2]

px = [p1[0], p2[0], p3[0], p4[0], p5[0], p6[0]]
py = [p1[1], p2[1], p3[1], p4[1], p5[1], p6[1]]

'''
Plot der Moden 1 und 2 
Graphische Darstellung mit Hilfe der Pfeile für Richtung und Amplitude
Grünes Kreuz markiert den gewählten Nullpunkt
'''
'''plot2, ay = plt.subplots(2,1)
plt.rc('xtick', labelsize=14) 
plt.rc('ytick', labelsize=14) 
plt.rc('axes', labelsize=16)
imgschaufel = plt.imread('Schaufel.JPG')
plot2.set_figwidth(12)
plot2.set_figheight(8)
ay[0].imshow(imgschaufel, extent = [82.1, -69.9, -15.8, int(152*451/1692)-15.8])
ay[0].plot(0,0,'g+')
ay[0].plot(px,py, 'r*')
ay[0].set_title('Mode 1')
ay[1].set_xlabel('x-Koordinate in mm')
ay[1].set_ylabel('y-Koordinate in mm')
plot2.tight_layout()
plotscale = 0.03
for i in count:
    ay[0].arrow(px[i], py[i], 0, plotscale*r_1[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')

ay[1].imshow(imgschaufel, extent = [82.1, -69.9, -15.8, int(152*451/1692)-15.8])
ay[1].plot(0,0,'g+')
ay[1].plot(px,py, 'r*')
ay[1].set_title('Mode 2')
for i in count: 
    ay[1].arrow(px[i], py[i], 0, plotscale*r_2[i],width = 0.1, head_width = 1, head_length = 2, ls = '-', ec = 'b', fc = 'b')
#plt.savefig(unterordner+'\Moden.svg', format = 'svg', dpi = 2000)'''


'''KONTROLLE DES ERMITTELTEN DÄMPFUNGSMAßES'''
p = 0
'''plot, D_intplot = plt.subplots()
plt.ylabel('Komplexe Amplitude')
plt.xlabel('Frequenz in Hz')
plt.title('Berechnung des Dämpfungsgrades')
plt.yscale('log')
plt.plot(fftfreq[5:1000] ,np.abs(HxF[p])[5:1000])
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new1l[p],4)*10000/(sweeplen*sampledt))], HxFmaxes1sqrt[p], 'r+')
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new1r[p],4)*10000/(sweeplen*sampledt))], HxFmaxes1sqrt[p], 'r+')
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new2l[p],4)*10000/(sweeplen*sampledt))], HxFmaxes2sqrt[p], 'r+')
plt.plot(np.arange(0,1000,0.0001)[int(np.round(x_new2r[p],4)*10000/(sweeplen*sampledt))], HxFmaxes2sqrt[p], 'r+')

plot, Dplot = plt.subplots()
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


'''VERGLEICH ALLER AMPLITUDENFREQUENZGÄNGE, DRUCK UND SAUGSEITE'''
'''plot7, hxfplot = plt.subplots()
plt.xlabel('Frequenz in Hz')
plt.ylabel('Amplitude in mm/N')
plt.plot(fftfreq[5:1000], np.abs(HxF[0][5:1000]), color='red')
plt.plot(fftfreq[5:1000], np.abs(HxF[1][5:1000]), color='red')
plt.plot(fftfreq[5:1000], np.abs(HxF[2][5:1000]), color='red')
plt.plot(fftfreq[5:1000], np.abs(HxF[3][5:1000]), color='green')
plt.plot(fftfreq[5:1000], np.abs(HxF[4][5:1000]), color='green')
plt.plot(fftfreq[5:1000], np.abs(HxF[5][5:1000]), color='green')
plt.yscale('log')
plot7.tight_layout
plot7.set_figwidth(10)
plot7.set_figheight(6)
#'''

'''Savefiles für tabellarische Auswertung'''
from tabulate import tabulate as tb
import pandas as pd


'''datasave = open(unterordner+'\Messung Stats.txt', 'w')
datasave.write('Erste Mode bei '+str(HxFeigfreqcombined_1)+'Hz (max Abweichung '+str((np.max(HxFeigvalue1)-np.min(HxFeigvalue1))/5)+'Hz)\n'+'Modale Daempfung '+str(round(D1_avg,5))+' \n'+'Zweite Mode bei '+str(+HxFeigfreqcombined_2)+'Hz (max Abweichung '+str((np.max(HxFeigvalue2)-np.min(HxFeigvalue2[HxFeigvalue2 != 0]))/5)+'Hz)\n'+'Modale Daempfung '+str(round(D2_avg,5))+'\n'+'Messauflösung '+str(resolution)+'Hz')
datasave.close()'''
table = [['Messpunkt', 1,2,3,4,5,6], 
         ['-','-','-','-','-','-','-'],
         ['Erste Eigenfrequenz',HxFeigfrequencies[0][0],HxFeigfrequencies[1][0],HxFeigfrequencies[2][0],HxFeigfrequencies[3][0],HxFeigfrequencies[4][0],HxFeigfrequencies[5][0]],
         ['Dämpfungsmaß 1te Mode', np.round(D1_int[0],5),np.round(D1_int[1],5),np.round(D1_int[2],5),np.round(D1_int[3],5),np.round(D1_int[4],5),np.round(D1_int[5],5)],
         ['Auslenkung 1te Mode', np.round(r_1_norm[0],5), np.round(r_1_norm[1],5), np.round(r_1_norm[2],5), np.round(r_1_norm[3],5), np.round(r_1_norm[4],5), np.round(r_1_norm[5],5)],
         ['-','-','-','-','-','-','-'],
         ['Zweite Eigenfrequenz',HxFeigfrequencies[0][1],HxFeigfrequencies[1][1],HxFeigfrequencies[2][1],HxFeigfrequencies[3][1],HxFeigfrequencies[4][1],HxFeigfrequencies[5][1]],
         ['Dämpfungsmaß 2te Mode', np.round(D2_int[0],5),np.round(D2_int[1],5),np.round(D2_int[2],5),np.round(D2_int[3],5),np.round(D2_int[4],5),np.round(D2_int[5],5)],
         ['Auslenkung 2te Mode', np.round(r_2_norm[0],5), np.round(r_2_norm[1],5), np.round(r_2_norm[2],5), np.round(r_2_norm[3],5), np.round(r_2_norm[4],5), np.round(r_2_norm[5],5)],
        ]
print(tb(table))
print('Erste Mode bei '+str(HxFeigfreqcombined_1)+'Hz\n'+'Modale Daempfung '+str(round(D1_avg,5))+' \n'+'Zweite Mode bei '+str(+HxFeigfreqcombined_2)+'Hz\n'+'Modale Daempfung '+str(round(D2_avg,5))) 

dataframe = pd.DataFrame(table)
dataframe.to_excel('auswertung_hammer.xlsx', index=False, header=False)
#'''
