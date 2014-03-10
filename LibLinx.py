#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
LibLinx.py
Librairie de traitement du signal du lock in numerique

Created by N. Bruyant 2014-03-08

"""
from scipy import zeros, cos, pi, sin, fft, real, mean, optimize,
fromfile, floor, exp, argmax, sqrt, diff, float32
from numpy import arange, linspace
from scipy.integrate import trapz, cumtrapz
from scipy.fftpack import fftfreq
import logging

#configure logging
logging.basicConfig(filename='LibLinx.log', level=logging.DEBUG)

def FindRefFreq(Data, SampleRate):
    """
     
    Parameters
    ==========
        data: np array
        SampleRate: sample rate in seconds-1

    Returns
    =======
        detected frequency in Hz with respect to samplerate
    """
    logging.info('Start reference frequency detection')

    #calculate the time vector
    #TODO Verify the coherence with the stored one
    time = arange(0, len(data))/SampleRate

    # On calcule les fréquences de modulation de chaque lock-in
    # On approche la valeur avec une FFT
    Data_fourier = abs(real(fft(Data)))
    freqs = fftfreq(len(Data), 1/SampleRate)

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
    F[1] = abs(trapz(Data[0:fourier_lenght]*exp(expon*freq),
                     dx=1/SampleRate))
    F[2] = abs(trapz(Data[0:fourier_lenght]*exp(expon*(freq+deltaf)),
                     dx=1/SampleRate))
    F[0] = abs(trapz(Data[0:fourier_lenght]*exp(expon*(freq-deltaf)),
                     dx=1/SampleRate))

    #optimize the frequency to maximize the DFT
    if F[2] > F[1]:
        essaimax = F[1]
        while abs(deltaf) > 0.0002:
            F[2] = abs(trapz(Data[0:fourier_lenght]*exp(expon*(freq+deltaf)),
                       dx=1/SampleRate))
            if F[2] > essaimax:
                essaimax = F[2]
                freq += deltaf
            else:
                deltaf = -deltaf/10
    elif F[0] > F[1]:
        deltaf = -deltaf
        essaimax = F[1]

        while abs(deltaf) > 0.0002:
            F[0] = abs(trapz(Data[0:fourier_lenght]*exp(expon*(freq+deltaf)),
                       dx=1/SampleRate))
            if F[0] > essaimax:
                essaimax = F[0]
                freq += deltaf
            else:
                deltaf = -deltaf/10
    else:
        logging.debug("You are lucky F[0]>F[1]>F[2] and deltaf is %f", deltaf)

    logging.info('Detected frequency is %f Hz' , freq)
    return freq


def data_analyze():
    ''' 

    Treat raw data
    '''
    self.sig_out = zeros((len(self.data[:, 2]), len(self.coldata)*2))
    phase_ref = zeros(len(self.data[:, 2]))
    antiphase_ref = zeros(len(self.data[:, 2]))
    self.amplitude = zeros(len(self.coldata))
    self.zp = int(self.lineEdit_21.text())
    
    for j in range(0, len(self.coldata)):
        if len(self.colref) == 1:
            posref = self.colref
            thefreq = self.freq[0]
        else:
            posref = self.colref[j]
            thefreq = self.freq[j]
        
        # Phase detector
        self.amplitude[j] = (max(self.data[1:10000, int(posref)])-min(self.data[1:10000, int(posref)]))/2
        phase_ref[:] = self.amplitude[j]*cos(2*pi*thefreq*self.data[:, 0])
        
        # Phase detector using sinusoidal fit
        fitfunc       = lambda p, x: self.amplitude[j]*cos(2*pi*thefreq*x+p[0])
        errfunc       = lambda p, x, y: fitfunc(p, x) - y
        p0            = [0]
        fit_le        = int(self.zp*round(self.sample_rate/thefreq))
        p1, success   = optimize.leastsq(errfunc, p0[:], args=(self.data[0:fit_le,0], self.data[0:fit_le, int(posref)]))
        self.phasemin = p1*180/pi
        
        phase_ref[:] = self.amplitude[j]*cos(2*pi*thefreq*self.data[:, 0]+self.phasemin[j]*pi/180+self.dephase*pi/180)
        self.sig_out[:, 2*j] = phase_ref[:]*self.data[:, self.coldata[j]]
        self.sig_out[:, 2*j] = self.data_lpass(self.sig_out[:, 2*j], self.cutoff, self.sample_rate)
        
        antiphase_ref[:] = -self.amplitude[j]*sin(2*pi*thefreq*self.data[:, 0]+self.phasemin[j]*pi/180+self.dephase*pi/180)
        self.sig_out[:, 2*j+1] = antiphase_ref[:]*self.data[:, self.coldata[j]]
        self.sig_out[:, 2*j+1] = self.data_lpass(self.sig_out[:, 2*j+1], self.cutoff, self.sample_rate)
           
                
if __name__ == '__main__':
    sample_rate = 500E3
    NSamples = 600000
    data = zeros((NSamples, 5))
    data[:, 0] = arange(0, NSamples)/sample_rate
    FreqsTest = linspace(1000, 20000, num=5)
    res = zeros(len(FreqsTest))
    for index, FreqRef in enumerate(FreqsTest):
        data[:, 1] = sin(data[:, 0]*FreqRef*2*pi)
        res[index] = FindRefFreq(data[:, 1], sample_rate)
    print res