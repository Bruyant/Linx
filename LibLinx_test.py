# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 21:58:26 2014

@author: nicolas
"""

import unittest
import logging

import random
from numpy import arange, linspace, zeros, sin, pi, max
from LibLinx import FindRefFreq

_multiprocess_can_split_ = True

tolerance = 0.001   # 0.01% -> i.e; 1e-4


def GenerateRefSignal(FreqRef=54321, sample_rate=500E3, NSamples=1000000):
    data = zeros((NSamples, 5))
    data[:, 0] = arange(0, NSamples)/sample_rate
    data[:, 1] = sin(data[:, 0]*FreqRef*2*pi)
    return data


def ComputeError(FreqRef, Res):
    return abs(max(Res-FreqRef))/max(FreqRef)*100


class TestFindRefFreq(unittest.TestCase):
    def setUp(self):
        # base setup code

        self.FreqRef = 54321
        self.sample_rate = 500E3
        self.NSamples = 1000000
        self.data = GenerateRefSignal(self.FreqRef,
                                    self.sample_rate,self.NSamples)
        logging.debug('FreqRef={0},sample_rate={1},NSamples={2}'.format(
                        self.FreqRef,self.sample_rate,self.NSamples))                         


class TestFindRefFreqBasic(TestFindRefFreq):
    def runTest(self):
        res=FindRefFreq(self.data[:,1],500E3)
        self.assertTrue(ComputeError(54321,res)<tolerance) 

class TestFindRefFreqSingleRandom(TestFindRefFreq):
    def runTest(self):
        #Single random test around 50KHz
        self.FreqRef=40E3+random.randint(0,20000)
        self.data[:,1]=sin(self.data[:,0]*self.FreqRef*2*pi)
        res=FindRefFreq(self.data[:,1],500E3)
        self.assertTrue(ComputeError(self.FreqRef,res)<tolerance)
        
class TestFindRefFreqSamplerateRandom(TestFindRefFreq):
    def runTest(self):
        #random samplerate
        self.sample_rate=200E3+random.randint(0,3E5)
        self.data=GenerateRefSignal(54321,self.sample_rate)
        self.res=FindRefFreq(self.data[:,1],self.sample_rate)
        print self.res,self.sample_rate
        self.assertTrue(ComputeError(54321,self.res)<tolerance)


class TestFindRefFreqLinspaceFreq(TestFindRefFreq):
    def runTest(self):
        self.sample_rate = 500E3
        self.data = GenerateRefSignal(54321, self.sample_rate)
        # define frequency to test
        FreqsTest = linspace(1000, 20000, num=5)
        res = zeros(len(FreqsTest))

        for index, FreqRef in enumerate(FreqsTest):
            self.data[:, 1] = sin(self.data[:, 0]*FreqRef*2*pi)
            res[index] = FindRefFreq(self.data[:, 1], self.sample_rate)

        self.assertTrue(ComputeError(FreqsTest, res) < tolerance)


class test_many_SamplesRates(TestFindRefFreq):
    def runTest(self):   
        for _ in range(200):
            TestFindRefFreqSamplerateRandom()


class test_many_Freqs(TestFindRefFreq):
    def runTest(self):
        for _ in range(200):
            TestFindRefFreqSingleRandom()

if __name__ == '__main__':
    #configure logging
    logging.basicConfig(filename='LibLinx_test.log', level=logging.DEBUG)
    unittest.main()