#!/usr/local/bin/python
'''
Created on Mar 03, 2016

@author: nathan
'''


import sys
sys.path.append('/Users/nathan/Documents/PYTHON')
from SaxsDisplay import SaxsPlot


if len(sys.argv) < 2:
    sys.exit('Useage: LogBin.py saxsfile_1.dat saxsfile_2.day... saxsfile_n.dat')
else:
    for file in sys.argv[1:]:
        try:
            if file[-4:] == '.dat':
                pass
            else:
                raise ValueError('File was not of type dat')
        except:
            sys.exit('Files must be of type dat')




job = SaxsPlot()
job.GetImages()
job.ParseDat()
job.LogBin()
job.OutputDat()
