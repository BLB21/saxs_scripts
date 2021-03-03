#!/usr/local/bin/python

import matplotlib.pyplot as plt
from os.path import isfile
#from matplotlib.backends.backend_pdf import PdfPages
import logging
import numpy

class autoRg():
    """Calculate Rg for a given dat file
    
    The autoRg class takes a dat file as input and calculates
    the Rg.

    
    """
    
    '''
    Constructor
    '''
    def __init__(self):
        ###start a log file
        self.logger = logging.getLogger('autoRg')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new autoRg job')
        ###set some parameters
        self.q = []
        self.i = []
        self.e = []

    def getDatFile(self, datfile=None):
        ############################################
        #USER SUPPLIED ARGS: USE THESE AS FILE LIST#
        ############################################
        if isfile(datfile) and datfile[-4:] == '.dat':
            self.parseDatFile(datfile)
        else:
            self.logger.error(datfile+' does not exist')

            
    def parseDatFile(self, datfile=None):
        self.logger.info('Reading and parsing dat file: '+str(datfile))
        qdata = []
        idata = []
        edata = []
                
        datafile = open(datfile, "r")
        saxsdata = datafile.readlines()
        for row in saxsdata:
            try:
                q = float(row.split()[0])
                i = float(row.split()[1])
                e = float(row.split()[2])
                if q > 0:
                    self.q.append(q)
                    self.i.append(i)
                    self.e.append(e)
            except:pass
        self.logger.info('Parsed with '+str(len(self.q))+' data points.')

    def plotData(self, x=None, y=None):
        if x == None:
            x = self.q
        if y == None:
            y = self.i
        fit, ax = plt.subplots()
        #ax.set_yscale("log", nonposy='clip')
        line, = ax.plot(x, y)
        plt.show()

    def gradientScan(self):
        window_size = 15
        gradients = []
        indeces = []
        for index in range(1, 100):#, q in enumerate(self.q[:100]):
            x = numpy.array([ numpy.power(q, 2) for q in self.q[index:index+window_size]], dtype=numpy.float)
            y = numpy.array([ numpy.log(i) for i in self.i[index:index+window_size]], dtype=numpy.float)
            A = numpy.vstack([x, numpy.ones(len(x))]).T
            m, c = numpy.linalg.lstsq(A, y)[0]
            gradients.append(m)
            indeces.append(index)
        self.plotData(indeces, gradients)

        
if __name__ == '__main__':

    job = autoRg()
    job.getDatFile('/Users/nathan/Documents/SAXS/TEST_DATA/bsa_10mgml_average_average_svd.dat')        
    job.gradientScan()
