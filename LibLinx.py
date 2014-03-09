#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
LibLinx.py
Librairie de traitement du signal du lock in numerique

Created by N. Bruyant 2014-03-08

"""
from scipy import zeros, cos, pi, sin, fft, real, mean, optimize,\
fromfile, floor, exp, argmax, sqrt, diff, float32
from numpy import arange,linspace
from scipy.integrate import trapz, cumtrapz
from scipy.fftpack import fftfreq


def FindRefFreq(data,SampleRate):
    time = arange(0,len(data))/SampleRate  
    # On calcule les fréquences de modulation de chaque lock-in
    # On approche la valeur avec une FFT
    data_fourier = abs(real(fft(data)))
    freqs = fftfreq(len(data), 1/SampleRate)
    # find maximum frequency
    posj = argmax(data_fourier[0:round(len(data_fourier)/2)+1])
    freq = freqs[posj]
    
    # Determine the freq spacing of the FFT to take 1/1000 as a deltaf
    deltaf = (freqs[2]-freqs[1])/1000
    
    # Select only first part of the signal
    #TODO to be supressed from this function
    fourier_lenght = 16084
    
    F = zeros(3)
    #Calculate prefactor
    expon = -2j*pi*time[0:fourier_lenght]

    #Calculate DFT for 3 point around the approx maximum
    F[1] = abs(trapz(data[0:fourier_lenght]*exp(expon*freq), dx=1/SampleRate))
    F[2] = abs(trapz(data[0:fourier_lenght]*exp(expon*(freq+deltaf)), dx=1/SampleRate))
    F[0] = abs(trapz(data[0:fourier_lenght]*exp(expon*(freq-deltaf)), dx=1/SampleRate))
    
    
    if F[2]>F[1]:
        essaimax = F[1]
        while abs(deltaf)>0.0002:
            F[2] = abs(trapz(data[0:fourier_lenght]*exp(expon*(freq+deltaf)),dx=1/SampleRate))
            if F[2]>essaimax:
                essaimax = F[2]
                freq = freq+deltaf
            else:
                deltaf = -deltaf/10
    elif F[0]>F[1]:
        deltaf = -deltaf
        essaimax = F[1]
        
        while abs(deltaf)>0.0002:
            F[0] = abs(trapz(data[0:fourier_lenght]*exp(expon*(freq+deltaf)),dx=1/SampleRate))
            if F[0]>essaimax:
                essaimax = F[0]
                freq = freq+deltaf
            else:
                deltaf = -deltaf/10
    else: 
        print "You are lucky ... nothing to do"
        
    return freq            
                
if __name__ == '__main__':
    sample_rate=500E3
    NSamples=600000
    data = zeros((NSamples,5))
    data[:,0] = arange(0,NSamples)/sample_rate
    FreqsTest=linspace(1000,20000,num=10)
    res=zeros(len(FreqsTest))
    for index,FreqRef in enumerate(FreqsTest):    
        data[:,1]=sin(data[:,0]*FreqRef*2*pi)
        res[index]=FindRefFreq(data[:,1],sample_rate)
    print res              