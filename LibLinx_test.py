# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 21:58:26 2014

@author: nicolas
"""

import unittest
import testtools
from numpy import arange,linspace,zeros,sin,pi

from LibLinx import FindRefFreq

def GenerateRefSignal(FreqRef=12345,sample_rate=500E3,NSamples=1000000):
    data = zeros((NSamples,5))
    data[:,0] = arange(0,NSamples)/sample_rate
    data[:,1]=sin(data[:,0]*FreqRef*2*pi)
    return data

class TestFindRefFreq(unittest.TestCase):
    def setUp(self):
        # base setup code
        self.data=GenerateRefSignal()

class TestFindRefFreqSingle(TestFindRefFreq):
    def runTest(self):
        #first test
        res=FindRefFreq(self.data[:,1],500E3)
        self.assertTrue(res-12345<0.01)
    
class TestFindRefFreq2Freq(TestFindRefFreq):
    def runTest(self):
        #mutiply by 2 the ref frequency
        self.data[:,1]=sin(self.data[:,0]*2*12345*2*pi)
        res=FindRefFreq(self.data[:,1],500E3)
        self.assertTrue(res-12345*2<0.01)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFindRefFreq)
    concurrent_suite = testtools.ConcurrentStreamTestSuite(lambda: ((case, None) for case in suite))
    concurrent_suite.run(testtools.StreamResult())