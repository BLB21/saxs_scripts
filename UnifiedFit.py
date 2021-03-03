#!/usr/local/bin/python
'''
#IMPORT LIBRARIES AND SET UP PATHS
'''
import sys
import logging
import numpy
from scipy.special import erf as erf
from os.path import isfile as isfile
import pylab
from PIL import Image
from subprocess import check_output
 
#import wx, sys, logging, wx.lib.dialogs, ssh, epics, urllib, ast, time, glob, re, os, subprocess, numpy, pylab, math, multiprocessing, getpass
from scipy import optimize
#from numpy import mat
#
#from os.path import isdir as isdir
#from subprocess import call, Popen, PIPE, STDOUT, check_output
#from SetUp import SetUp as SetUp
#from redisobj import RedisHashMap


class UnifiedFit():
    """Performs a unified fit on a dat file
    
    Does a unified fit on a dat file. Options to run datrg to get starting
    value of Rg and I(0). 
    
    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        try:
            self.options = dict(options)
        except:
            sys.exit('UnifiedFit class needs an input dictionary, none provided')

        ###start a log file
        self.logger = logging.getLogger('UnifiedFit')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new unified fit job')

    def RunAutoRg(self):
        if self.options['autorg']:
            command = 'autorg -f csv '+self.options['datfile']
            output = check_output(command, shell=True).rstrip()
            try:
                autorg_headers = output.split('\n')[0].split(',')
                autorg_data = output.split('\n')[-1].split(',')
                autorg_dict = {}
                for (index,header) in enumerate(autorg_headers):
                    autorg_dict[header] = autorg_data[index]
                self.options['rg'] = float(autorg_dict['Rg'])
                self.options['io'] = float(autorg_dict['I(0)'])
                self.logger.info('Ran autorg, with quality: '+autorg_dict['Quality'])
            except:
                self.logger.error('autorg failed! Proceeding with defaults')

    def ParseDat(self):
        if isfile(self.options['datfile']) and self.options['datfile'][-4:] == '.dat':
            self.logger.info('Parsing '+self.options['datfile'])
        else:
            self.logger.error(self.options['datfile']+' does not exist.')
            sys.exit()

        self.q = []
        self.i = []
                    
        datafile = open(self.options['datfile'], "r")
        saxsdata = datafile.readlines()
        for row in saxsdata:
            try:
                q = float(row.split()[0])
                i = float(row.split()[1])
                e = float(row.split()[2])
                if q > 0:
                    self.q.append(q)
                    self.i.append(i)
            except:pass
        self.logger.info('Parsed with '+str(len(self.q))+' data points.')


    def FitUnified(self):
        self.logger.info('Fitting unified equation')
        
        x = numpy.array(self.q)
        y = numpy.array(self.i)
        #G, Rg, B, P, Bkg
        def f(x, G, Rg, B, P, Bkg):
            return ((G*numpy.exp(-(x**2)*(Rg**2)/3))+(P*(((erf(x*B/6**0.5))**3)/x)**P)+Bkg)
        def resid_f(p, y, x):
            G, Rg, B, P, Bkg = p
            return y - f(x, G, Rg, B, P, Bkg)

#        def BandP(x, B, P):
#            return ((self.options['io']*numpy.exp(-(x**2)*(self.options['rg']**2)/3))+(P*((((erf(x*B/6**0.5)))**3)/x)**P)+self.options['bkg'])
#        def resid_bandp(p, y, x):
#            B, P = p
#            return y - BandP(x, B, P)


        G0, Rg0, B0, P0, Bkg0 = self.options['io'], self.options['rg'], 1.0E-07, self.options['porod'], 0.0

        [G, Rg, B, P, Bkg], flag = optimize.leastsq(resid_f, [G0, Rg0, B0, P0, Bkg0], args=(y, x))
#        [B, P], flag = optimize.leastsq(resid_bandp, [B0, P0], args=(y, x))

        self.logger.info('Rg Unified = '+str(Rg)+': Rg Guinier = '+str(self.options['rg']))
        self.logger.info('Porod Slope Unified = '+str(P)+': Starting value = '+str(self.options['porod']))
        self.logger.info('B = '+str(B))
        self.logger.info('I(0) Unified = '+str(G)+': I(0) Guinier = '+str(self.options['io']))
        self.logger.info('Background = '+str(Bkg)+': Starting value = '+str(Bkg0))

        self.options['b'] = B
        self.options['porod'] = P
        self.options['io'] = G
        self.options['bkg'] = Bkg
        self.options['rg'] = Rg

        ###CREATE A PYLAB PLOT OBJECT
        self.myfigure = pylab.figure()  
        self.myplot = self.myfigure.add_subplot(111)
        self.myplot.set_title('Deviation of beam centre X (blue) and Y (green) as a function of distance')
        self.myplot.set_xlabel('Distance (mm)')
        self.myplot.set_ylabel('Beam centre (pixels)')
        self.myplot.set_yscale('log', nonposy='clip')
        self.myplot.set_xscale('log')

        #plot the datapoints
        self.myplot.plot(x, y, 'ro')
        
        # plot the smooth model fit
        xc = numpy.linspace(x.min(), x.max(), 150)
        self.myplot.plot(xc, f(xc, self.options['io'], self.options['rg'], self.options['b'], self.options['porod'], self.options['bkg'] ))
        
        
        #save the plot to file
        output_image='unifiedfit.png'
        self.myfigure.savefig(output_image, format="png")
        img = Image.open('unifiedfit.png')
        img.show()

if __name__ == '__main__':

    '''
    parse command line options
    '''
    if len(sys.argv) < 2:
        sys.argv.append('-h')

    from optparse import OptionParser
    from optparse import OptionGroup
    
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-f", "--datfile", action="store", type="string", dest="datfile", help="The dat file you want to perform the fit on.")

    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-r", "--rg", action="store", type="float", dest="rg", default=20.0, help="A starting value for Rg. Default is 20. If option to run datrg is checked then this value will overwrite the starting value.")
    optional.add_option("-b", "--bkg", action="store", type="float", dest="bkg", default=0.0, help="A starting value for background level. Default is 0.")
    optional.add_option("-i", "--io", action="store", type="float", dest="io", default=0.01, help="A starting value for I(0). Default is 0.01. If option to run datrg is checked then this value will overwrite the starting value.")
    optional.add_option("-p", "--porod", action="store", type="float", dest="porod", default=3.5, help="A starting value for the Porod slope. Default is 3.5.")
    optional.add_option("-a", "--autorg", action="store_true", dest="autorg", default=False, help="Run datrg to get a starting value for Rg and I(0). Default is to not run it. This option acts as a switch and does not need a value.")
    optional.add_option("-g", "--graph", action="store_true", dest="graph", default=False, help="Graph the result, default is false, This option acts as a switch and does not need a value.")

    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    

    options = eval(str(options))
    job = UnifiedFit(options)
    job.RunAutoRg()
    job.ParseDat()
    job.FitUnified()
