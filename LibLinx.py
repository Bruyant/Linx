#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
LibLinx.py
Librairie de traitement du signal du lock in numerique

Created by N. Bruyant 2014-03-08

"""
from scipy import zeros, cos, pi, sin, fft, real, mean, optimize, fromfile, floor, exp, argmax, sqrt, diff, float32
from numpy import arange, linspace
from scipy.integrate import trapz, cumtrapz
from scipy.fftpack import fftfreq
from scipy.signal import resample

from scipy.signal import buttord, ellipord, cheb1ord, ellip, cheby1, filtfilt, butter
# from modified_filter import butter  # error free scipy butter

import logging

# configure logging
logging.basicConfig(filename='LibLinx.log', level=logging.DEBUG)

# setup tolerences
deltaf_Tol = 0.0002


def FindRefFreq(RefData, SampleRate):
    """

    Parameters
    ==========
        RefData: np array
        SampleRate: sample rate in seconds-1

    Returns
    =======
        detected frequency in Hz with respect to samplerate
    """
    logging.info('Start reference frequency detection')

    # calculate the time vector
    # TODO Verify the coherence with the stored one
    time = arange(0, len(RefData)) / SampleRate

    # On calcule les fréquences de modulation de chaque lock-in
    # On approche la valeur avec une FFT
    data_fourier = abs(real(fft(RefData)))
    freqs = fftfreq(len(RefData), 1 / SampleRate)

    # find maximum frequency
    posj = argmax(data_fourier[0:round(len(data_fourier) / 2) + 1])
    freq = freqs[posj]

    # Determine the freq spacing of the FFT to take 1/1000 as a deltaf
    deltaf = (freqs[2] - freqs[1]) / 1000

    # Select only first part of the signal
    # TODO calculate the good fourier length with respect to the detected freq
    # Downsampling to 2*Shannon freq
# FIXME
#    downsampling= round(len(data_fourier)/posj/4)
#    logging.debug("downsampling %i",downsampling)
#    # precision visée 1E-6    ->  posj*1E_6*fourier_lenght=1/2
#    fourier_lenght = min(len(Data),1/(2*posj*1E-9)*downsampling)
#    logging.debug("fourier_lenght %i", fourier_lenght)
    fourier_lenght = 36000

    #cut signal
    RefData = RefData[0:fourier_lenght]
    time = time[0:fourier_lenght]

    #downsample signal
#    Data = resample(Data,downsampling)
#    time = resample(time,downsampling)
#    SampleRate=SampleRate/downsampling

    #initate error vector
    F = zeros(3)
    # Calculate prefactor
    expon = -2j * pi * time

    # Calculate DFT for 3 point around the approx maximum
    F[1] = abs(trapz(RefData * exp(expon * freq),
                     dx=1 / SampleRate))
    F[2] = abs(trapz(RefData * exp(expon * (freq + deltaf)),
                     dx=1 / SampleRate))
    F[0] = abs(trapz(RefData * exp(expon * (freq - deltaf)),
                     dx=1 / SampleRate))

    # optimize the frequency to maximize the DFT
    if F[2] > F[1]:
        essaimax = F[1]
        while abs(deltaf) > deltaf_Tol:
            F[2] = abs(trapz(RefData * exp(expon * (freq + deltaf)),
                       dx=1 / SampleRate))
            if F[2] > essaimax:
                essaimax = F[2]
                freq += deltaf
            else:
                deltaf = -deltaf / 10
    elif F[0] > F[1]:
        deltaf = -deltaf
        essaimax = F[1]

        while abs(deltaf) > deltaf_Tol:
            F[0] = abs(trapz(RefData * exp(expon * (freq + deltaf)),
                       dx=1 / SampleRate))
            if F[0] > essaimax:
                essaimax = F[0]
                freq += deltaf
            else:
                deltaf = -deltaf / 10
    else:
        logging.debug("strange case where F[0]>F[1]>F[2] and deltaf is %f", deltaf)

    logging.info('Detected frequency is %f Hz +/- %f Hz', freq, deltaf)
    return freq


def data_analyze(data, RefData, thefreq, SampleRate, cutoff=None,
                 Phase_accuracy=2, dephase=0):
    '''
    detection of the lock-in signal with a reconstructed reference signal

    Parameters :
        data: 1D np array
        freq: carrier detected frequency
        SampleRate: sample rate in Hertz
        cutoff: Low pass frequency cutoff in Hertz
        Phase_accuracy = 2
        dephase = 0 manual dephasing correction

    Returns :
        sig_out: 2D np array demodulated signal in phase and out-of-phase
    '''
    datalen = len(data)
    sig_out = zeros((datalen,2))
    phase_ref = zeros(datalen)
    antiphase_ref = zeros(datalen)

    #time vector
    time = arange(0, len(data)) / SampleRate

    if cutoff is None:  # if no cutoff definded choose a good guess
        cutoff = thefreq / 2

    # Phase detector
    # calculate the amplitude of the reference signal
    amplitude = (max(RefData[1:10000]) - min(RefData[1:10000])) / 2

    # calculte the phase of the reference
    #phase_ref = amplitude * cos(2 * pi * thefreq * data[:, 0])

    # Phase detector using sinusoidal fit of the ref signal
    #------------------------------------------
    #define the fitting function
    fitfunc = lambda p, x: amplitude * cos(2 * pi * thefreq * x + p[0])
    errfunc = lambda p, x, y: fitfunc(p, x) - y

    #init the phase to 0
    p0 = [0]

    # define the fit legnth
    fit_le = int(Phase_accuracy * round(SampleRate / thefreq))

    # Calculate the phase p1
    p1, success = optimize.leastsq(errfunc, p0[:], args=(time[0:fit_le],
                                   RefData[0:fit_le]))
    phasemin = p1 * 180 / pi

    # re-calculate the reference signal
    phase_ref = amplitude * cos(2 * pi * thefreq * time + phasemin * pi / 180 +
                                dephase * pi / 180)

    # multiplication of the ref by the signal
    sig_out[:, 0] = phase_ref * data
    # lowpass filtering
    sig_out[:, 0] = data_lpass(sig_out[:, 0], cutoff, SampleRate)

    # calculate the reference signal dephased by pi/2
    antiphase_ref[:] = -amplitude * sin(2 * pi * thefreq * time +
                                        phasemin * pi / 180 + dephase * pi / 180)
    # multiplication of the dephased ref by the signal
    sig_out[:, 1] = antiphase_ref[:] * data
    # lowpass filtering
    sig_out[:, 1] = data_lpass(sig_out[:, 1], cutoff, SampleRate)

    return sig_out
    
def data_lpass( x, Wp, srate, Ws= 2, Rp= 2,Rs= 24):
    ''' Low-pass filter '''
    Wp = float(Wp*2/srate)
    (b, a)  =  butter(5, Wp, btype = 'low')
    y  =  filtfilt(b, a, x)
    
    return(y) 

if __name__ == '__main__':
    sample_rate = 250E3
    NSamples = 3000000
    data = zeros((NSamples, 2))
    data[:, 0] = arange(0, NSamples) / sample_rate
    FreqRef = 12345.678905
    RefData= sin(data[:, 0] * FreqRef * 2 * pi)
    thefreq=FindRefFreq(RefData, sample_rate)
    print("frequency error")
    print(thefreq-FreqRef)
    # test of demodultion
    #signal creation
    ModAmp=(1 + 0.1 * abs(sin(data[:, 0] * 20 * 2 * pi)))
    data[:, 1] = ModAmp * sin(data[:, 0] * FreqRef * 2 * pi+1)
    #plot(data[1:10000, 0], data[1:10000, 1])
    
    res=data_analyze(data[:,1], RefData, FreqRef, sample_rate)
    
    from pylab import plot
    idxmax=int(50E3)
    plot(res[1:idxmax, 0])
    plot(ModAmp[:idxmax]/2)
