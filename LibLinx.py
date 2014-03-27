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
import logging

# configure logging
logging.basicConfig(filename='LibLinx.log', level=logging.DEBUG)

# setup tolerences
deltaf_Tol = 0.0002


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

    # calculate the time vector
    # TODO Verify the coherence with the stored one
    time = arange(0, len(Data)) / SampleRate

    # On calcule les fréquences de modulation de chaque lock-in
    # On approche la valeur avec une FFT
    data_fourier = abs(real(fft(Data)))
    freqs = fftfreq(len(Data), 1 / SampleRate)

    # find maximum frequency
    posj = argmax(data_fourier[0:round(len(data_fourier) / 2) + 1])
    freq = freqs[posj]

    # Determine the freq spacing of the FFT to take 1/1000 as a deltaf
    deltaf = (freqs[2] - freqs[1]) / 1000

    # Select only first part of the signal
    # TODO calculate the good fourier length with respect to the detected freq
    # Downsampling to 2*Shannon freq
#FIXME
#    downsampling= round(len(data_fourier)/posj/4)
#    logging.debug("downsampling %i",downsampling)
#    # precision visée 1E-6    ->  posj*1E_6*fourier_lenght=1/2
#    fourier_lenght = min(len(Data),1/(2*posj*1E-9)*downsampling)
#    logging.debug("fourier_lenght %i", fourier_lenght)
    fourier_lenght = 16000

    #cut signal
    Data = Data[0:fourier_lenght]
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
    F[1] = abs(trapz(Data * exp(expon * freq),
                     dx=1 / SampleRate))
    F[2] = abs(trapz(Data * exp(expon * (freq + deltaf)),
                     dx=1 / SampleRate))
    F[0] = abs(trapz(Data * exp(expon * (freq - deltaf)),
                     dx=1 / SampleRate))

    # optimize the frequency to maximize the DFT
    if F[2] > F[1]:
        essaimax = F[1]
        while abs(deltaf) > deltaf_Tol:
            F[2] = abs(trapz(Data * exp(expon * (freq + deltaf)),
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
            F[0] = abs(trapz(Data * exp(expon * (freq + deltaf)),
                       dx=1 / SampleRate))
            if F[0] > essaimax:
                essaimax = F[0]
                freq += deltaf
            else:
                deltaf = -deltaf / 10
    else:
        logging.debug("strange case where F[0]>F[1]>F[2] and deltaf is %f", deltaf)

    logging.info('Detected frequency is %f Hz', freq)
    return freq


def data_analyze(data, freq):
    '''

    Treat raw data with respect to the carrier frequency
    '''
    sig_out = zeros((len(data[:, 2]), len(coldata) * 2))
    phase_ref = zeros(len(data[:, 2]))
    antiphase_ref = zeros(len(data[:, 2]))
    amplitude = zeros(len(coldata))
    zp = int(lineEdit_21.text())

    for j in range(0, len(coldata)):
        if len(colref) == 1:
            posref = colref
            thefreq = freq[0]
        else:
            posref = colref[j]
            thefreq = freq[j]

        # Phase detector
        amplitude[j] = (max(data[1:10000, int(posref)]) - min(data[1:10000, int(posref)])) / 2
        phase_ref[:] = amplitude[j] * cos(2 * pi * thefreq * data[:, 0])

        # Phase detector using sinusoidal fit
        fitfunc = lambda p, x: amplitude[j] * cos(2 * pi * thefreq * x + p[0])
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        p0 = [0]
        fit_le = int(zp * round(sample_rate / thefreq))
        p1, success = optimize.leastsq(errfunc, p0[:], args=(data[0:fit_le, 0], data[0:fit_le, int(posref)]))
        phasemin = p1 * 180 / pi

        phase_ref[:] = amplitude[j] * cos(2 * pi * thefreq * data[:, 0] + phasemin[j] * pi / 180 + dephase * pi / 180)
        sig_out[:, 2 * j] = phase_ref[:] * data[:, coldata[j]]
        sig_out[:, 2 * j] = data_lpass(sig_out[:, 2 * j], cutoff, sample_rate)

        antiphase_ref[:] = -amplitude[j] * sin(2 * pi * thefreq * data[:, 0] + phasemin[j] * pi / 180 + dephase * pi / 180)
        sig_out[:, 2 * j + 1] = antiphase_ref[:] * data[:, coldata[j]]
        sig_out[:, 2 * j + 1] = data_lpass(sig_out[:, 2 * j + 1], cutoff, sample_rate)


if __name__ == '__main__':
    sample_rate = 500E3
    NSamples = 3000000
    data = zeros((NSamples, 5))
    data[:, 0] = arange(0, NSamples) / sample_rate
    FreqRef = 1000.123
    data[:, 1] = sin(data[:, 0] * FreqRef * 2 * pi)
    print FindRefFreq(data[:, 1], sample_rate)-FreqRef

