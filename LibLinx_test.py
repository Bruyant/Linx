# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 21:58:26 2014

@author: nicolas
"""

import unittest
from numpy import arange,linspace,zeros,sin,pi

def GenerateRefSignal(FreqRef):
    sample_rate=500E3
    NSamples=1000000
    data = zeros((NSamples,5))
    data[:,0] = arange(0,NSamples)/sample_rate
    data[:,1]=sin(data[:,0]*FreqRef*2*pi)
    return data

class SimpleWidgetTestCase(unittest.TestCase):
    def setUp(self):
        # base setup code
        self.n=1

class DefaultWidgetSizeTestCase(SimpleWidgetTestCase):
    def runTest(self):
        #first test
        self.assertEqual(self.n,1)
    
class WidgetResizeTestCase(SimpleWidgetTestCase):
    def runTest(self):
        #second test
        self.assertEqual(self.n,2)

if __name__ == '__main__':
    unittest.main()