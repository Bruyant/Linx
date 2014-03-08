#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
LibLinx.py
Librairie de traitement du signal du lock in numerique

Created by N. Bruyant 2014-03-08

"""
from scipy import zeros, cos, pi, sin, fft, real, mean, optimize,\
fromfile, floor, exp, argmax, sqrt, diff, float32
from scipy.integrate import trapz, cumtrapz
from scipy.fftpack import fftfreq


def FindRefFreq(data):    
    # On calcule les fréquences de modulation de chaque lock-in
    # On approche la valeur avec une FFT
    data_fourier = abs(real(fft(data[:,1])))
    freq = fftfreq(len(data), 1/sample_rate)
    
    posj = argmax(data_fourier[0:round(len(data_fourier)/2)+1])
    
    freq[i] = freq[posj]
    
    # On effectue une DFT pour affiner le calcul
    deltaf = (freq[2]-freq[1])/1000
    fourier_lenght = 16084
    
    F = [0]*3
    expon = -2j*pi*data[0:fourier_lenght, 0]
    
    F[1] = abs(trapz(data[0:fourier_lenght,1]*exp(expon*freq[i]), data[0:fourier_lenght, 0]))
    F[2] = abs(trapz(data[0:fourier_lenght,1]*exp(expon*(freq[i]+deltaf)), data[0:fourier_lenght, 0]))
    F[0] = abs(trapz(data[0:fourier_lenght,1]*exp(expon*(freq[i]-deltaf)), data[0:fourier_lenght, 0]))
    
    
    if F[2]>F[1]:
        essaimax = F[1]
        while abs(deltaf)>0.0002:
            F[2] = abs(trapz(data[0:fourier_lenght,1]*exp(expon*(freq[i]+deltaf)),data[0:fourier_lenght, 0]))
            if F[2]>essaimax:
                essaimax = F[2]
                freq[i] = freq[i]+deltaf
            else:
                deltaf = -deltaf/10
    elif F[0]>F[1]:
        deltaf = -deltaf
        essaimax = F[1]
        
        while abs(deltaf)>0.0002:
            F[0] = abs(trapz(data[0:fourier_lenght,1]*exp(expon*(freq[i]+deltaf)), data[0:fourier_lenght, 0]))
            if F[0]>essaimax:
                essaimax = F[0]
                freq[i] = freq[i]+deltaf
            else:
                deltaf = -deltaf/10
                
                
if __name__ == '__main__':
    sample_rate=500E3
    NSamples=1000000
    data = zeros((NSamples,5))
    for i in range(0,NSamples):
        data[i, 0] = i/sample_rate    
    data[:,1]=sin(data[:,0]*1E4)
    FindRefFreq(data)               