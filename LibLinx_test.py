# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 21:58:26 2014

@author: nicolas
"""

import unittest

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