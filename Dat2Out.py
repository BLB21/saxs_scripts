#!/anaconda/bin/python

'''
Created on Jul 24, 2019

@author: nathan
'''

import sys
import os
import logging
import numpy as np


class Dat2Out():
    """A dat file and a P(r) dat file from Scatter convert to an out file
    
    When you do the P(r) refine or 2file steps in Scatter and you don't
    have ATSAS it will output two dat files, one with the P(r) in it and
    the other with the fitted data (q, I_OBS(q), error, I_CALC(q)). Dat2Out
    reads in BOTH of these files and outputs an out file that should work in
    ATSAS dummy atom modelling programs or DENS.

    """
    
    '''
    Constructor
    '''
    def __init__(self):
        ###start a log file
        self.logger = logging.getLogger('Dat2Out')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a Dat2Out job') 

        ###Some parameters
        self.dat_pr_file = None
        self.dat_data_file = None
        self.dat_out_file = None
        self.pr_array = []
        self.dat_array = []
        self.output_text = []
        
    def setDatPr(self, dat_pr=None):
        if dat_pr:
            if os.path.isfile(dat_pr):
                try:
                    if dat_pr[-4:] == '.dat':
                        self.dat_pr_file = dat_pr
                        self.logger.info('Set dat Pr file')
                        return True
                    else:
                        self.logger.error('setDatPr needs a file of type .dat')
                        return False
                except:
                    self.logger.error('please include the .dat filename extensions')
                    return False
            else:
                self.logger.error('setDatPr file does not exist')
                return False
        else:
            self.logger.error('setDatPr function needs a filename as input')
            return False
        
    def setDatData(self, dat_data=None):
        if dat_data:
            if os.path.isfile(dat_data):
                try:
                    if dat_data[-4:] == '.dat':
                        self.dat_data_file = dat_data
                        self.logger.info('Set dat data file')
                        return True
                    else:
                        self.logger.error('setDatData needs a file of type .dat')
                        return False
                except:
                    self.logger.error('please include the .dat filename extensions')
                    return False
            else:
                self.logger.error('setDatData file does not exist')
                return False
        else:
            self.logger.error('setDatData function needs a filename as input')
            return False

    def parseDatPr(self):
        if self.dat_pr_file:
            with open(self.dat_pr_file, 'r') as f:
                d = f.readlines()
            for line in d:
                line = line.split()
                try:
                    r = float(line[0])
                    p = float(line[1])
                    e = float(line[2])
                    self.pr_array.append( (r,p,e) )
                except:
                    pass
            self.logger.info('Parsed Pr data with '+str(len(self.pr_array))+' points')

    def parseDatData(self):
        if self.dat_data_file:
            with open(self.dat_data_file, 'r') as f:
                d = f.readlines()
            for line in d:
                line = line.split()
                try:
                    q = float(line[0])
                    iobs = float(line[1])
                    e = float(line[2])
                    icalc = float(line[3])
                    self.dat_array.append( (q,iobs,e,icalc) )
                except:
                    pass
            self.logger.info('Parsed dat data with '+str(len(self.dat_array))+' points')

    def formatOutput(self):
        header = [
            '           ####    G N O M   ---   Version 4.5a revised 09/02/02     ####',
            '',
            ' '*50+datetime.now().strftime('%a %b %d %H:%M:%S %Y'),
            '           ===    Run No   1   ===',
            ' Run title:  REMARK 265',
            '',
            '',
'  Number of points omitted at the end:         310

          To optimize the performance, some input data points 
           have been joined. New number of points equals to 484

   *******   Working Directory :
   *******    Input file(s) : /Users/nathan/bsa_1_average_average_refi
           Condition P(rmin) = 0 is used. 
           Condition P(rmax) = 0 is used. 

          Experimental conditions for the first run
          -----------------------------------------

          Highest ALPHA is found to be   0.5409E+03

  The measure of inconsistency AN1 equals to    0.9389E+00
     Alpha    Discrp  Oscill  Stabil  Sysdev  Positv  Valcen    Total   

        
        
if __name__ == '__main__':
    job = Dat2Out()
    job.setDatPr('bsa_1_average_average_refined_pr.dat')
    job.setDatData('bsa_1_average_average_refined_sx.dat')
    job.parseDatPr()
    job.parseDatData()
    job.logger.info('FINISHED NORMALLY')

