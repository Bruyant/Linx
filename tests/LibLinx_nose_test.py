# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 21:58:26 2014

@author: nicolas
"""
from nose.tools import *
import logging

# libs to generate tests
from random import *
from numpy import arange, linspace, zeros, sin, pi, max
import itertools

# lib under test
from LibLinx import FindRefFreq

# Generals test parameters
_multiprocess_can_split_ = True
tolerance = 2E-6   # 0.01% -> i.e; 1e-4
deltaf_Tol = 0.0002


# test bench init and closing
def setup():
    logging.basicConfig(filename='LibLinx_test.log', level=logging.DEBUG)
    print "Setup Done!"


def teardown():
    print "TEAR DOWN!"


#  Helpers functions
def GenerateRefSignal(FreqRef=54321, sample_rate=500E3, NSamples=1000000):
    data = GenerateRefSignalWithPhase(FreqRef, sample_rate, NSamples, Phase=0)
    return data


#  Helpers functions
def GenerateRefSignalWithPhase(FreqRef=54321, sample_rate=500E3,
                               NSamples=1000000, Phase=0):
    data = zeros((NSamples, 5))
    data[:, 0] = arange(0, NSamples) / sample_rate
    data[:, 1] = sin(data[:, 0] * FreqRef * 2 * pi + Phase)
    return data


def RelativeError(FreqRef, Res):
    return abs(max(Res - FreqRef)) / max(FreqRef)


def AbsError(FreqRef, Res):
    return abs(max(Res - FreqRef))


def check_FindRefFreq(FreqRef=54321, sample_rate=500E3,
                      NSamples=1000000):
    data = GenerateRefSignal(FreqRef, sample_rate, NSamples)

    # do the test
    res = FindRefFreq(data[:, 1], sample_rate)
    logging.debug('FreqRef={0},sample_rate={1},NSamples={2}'.format(
        FreqRef, sample_rate, NSamples))
    logging.debug('FreqDelta={0},FreqDelta/Freq={1}'.format(
        FreqRef - res, RelativeError(FreqRef, res)))
    assert AbsError(FreqRef, res) < deltaf_Tol or RelativeError(FreqRef, res) < tolerance


def check_FindRefFreqPhase(FreqRef=54321, sample_rate=500E3,
                           NSamples=1000000, Phase=0):
    data = GenerateRefSignalWithPhase(FreqRef, sample_rate, NSamples, Phase)
    # do the test
    res = FindRefFreq(data[:, 1], sample_rate)
    logging.debug('FreqRef={0},sample_rate={1},NSamples={2}, Phase={3}'.format(
        FreqRef, sample_rate, NSamples, Phase))
    logging.debug('FreqDelta={0},FreqDelta/Freq={1}'.format(
        FreqRef - res, RelativeError(FreqRef, res)))
    assert AbsError(FreqRef, res) < deltaf_Tol or RelativeError(FreqRef, res) < tolerance

# test functions


def test_FindRefFreqMultiple():
    FreqRef = arange(12345.678901, 100E3, 10000)
    sample_rate = (500E3, 250E3, 100E3, 50e3)
    NSamples = (2000000, 1000000, 100000, 50000)
    caseslists = [FreqRef,
                  sample_rate,
                  NSamples]
    # generate cross product
    cases = list(itertools.product(*caseslists))
    logging.debug(cases)
    # clean non physical tests
    goodcases = []
    for FreqRef, sample_rate, NSamples in cases:
        if FreqRef < sample_rate / 2:
            goodcases.append((FreqRef, sample_rate, NSamples))
    # start tests
    for FreqRef, sample_rate, NSamples in goodcases:
        yield check_FindRefFreq, FreqRef, sample_rate, NSamples


def test_FindRefFreqRandom():
    # add some randomness to improve coverage
    FreqRef = (12345 + random() * 50E3 for r in xrange(20))
    sample_rate = (250E3 + random() * 250E3 for r in xrange(5))
    NSamples = (100000, 90000)
    caseslists = [FreqRef,
                  sample_rate,
                  NSamples]
    # generate cross product
    cases = list(itertools.product(*caseslists))

    # clean non physical tests
    goodcases = []
    for FreqRef, sample_rate, NSamples in cases:
        if FreqRef < sample_rate / 2:
            goodcases.append((FreqRef, sample_rate, NSamples))
    # start tests
    for FreqRef, sample_rate, NSamples in goodcases:
        yield check_FindRefFreq, FreqRef, sample_rate, NSamples


def test_FindRefFreqPhase():
    # generate cross product
    FreqRef = 60000
    sample_rate = 500E3
    NSamples = 20000
    Phases = arange(0, pi, pi / 20)
    # start tests
    for Phase in Phases:
        yield check_FindRefFreqPhase, FreqRef, sample_rate, NSamples, Phase
