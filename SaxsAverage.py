#!/anaconda/bin/python
'''
Created on Apr 07, 2014

@author: nathan
'''

import sys
import os
import logging
import numpy as np


if len(sys.argv) < 2:
    print 'Useage: SaxsAverage.py file1.dat file2.dat ... fileN.dat'
    print 'SaxsAdd.py will average together all pdb files in the list'
    print 'and write to standard out.'
    sys.exit()

filelist = sys.argv[1:]

class SaxsAdd():
    """Sum together a list of SAXS dat files
    
    SaxsAdd sums together each file in a list and writes the result
    to standard out. The output is a dat file that can be piped to a
    file.

    """
    
    '''
    Constructor
    '''
    def __init__(self, filelist):
        ###start a log file
        self.logger = logging.getLogger('SaxsAdd')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new SaxsAdd job.') 
        self.filelist = filelist
        self.imagedict = {}
        
    def TestFiles(self):
        for file in self.filelist:
            if not os.path.isfile(file):
                self.logger.error(file+' does not exist')
                sys.exit()

    def ParseFiles(self):
        if len(self.filelist) < 2:
            self.logger.error('There are less than 2 valid files, cannot sum')
            sys.exit()

        
        for file in self.filelist:
            if os.path.splitext(file)[-1] == '.dat':
                qdata = []
                idata = []
                edata = []
                datafile = open(file, "r")
                saxsdata = datafile.readlines()
                lo_q = []
                hi_q = []
                no_datapoints = len(saxsdata)
                for index, row in enumerate(saxsdata):
                    try:
                        qdata.append(float(row.split()[0]))
                        idata.append(float(row.split()[1]))
                        edata.append(float(row.split()[2]))
                        if no_datapoints*0.9 < index < no_datapoints*0.95:
                            hi_q.append( float(row.split()[1]) )
                        elif no_datapoints*0.015 < index < no_datapoints*0.05:
                            lo_q.append( float(row.split()[1]) )
                        else:
                            pass
                    except:pass
                self.imagedict[os.path.splitext(file)[0]] = {'Q': qdata, 'I': idata, 'E': edata, 'lo_q': sum(lo_q), 'hi_q': sum(hi_q), 'outlier': False}
                self.logger.info('Parsed '+str(os.path.splitext(file)[0])+' with '+str(len(qdata))+' data points.')
                testset = set(self.imagedict[self.imagedict.keys()[0]]['Q'])
                for file in self.imagedict.keys()[1:]:
                    targetset = set(self.imagedict[file]['Q'])
                    if not len(testset - targetset) == 0:
                        self.logger.error('The files do not all have the same Q range')
                        sys.exit()


    def rejectOutliers(self):
        #remove air shots
        target = max([self.imagedict[x]['hi_q'] for x in self.imagedict.keys()])
        for image in self.imagedict.keys():
            if self.imagedict[image]['hi_q'] < target*0.98:
                self.imagedict[image]['outlier'] = True
            else:
                self.imagedict[image]['outlier'] = False
        self.logger.info(str(sum([self.imagedict[x]['outlier'] for x in self.imagedict.keys()]))+' of '+str(len(self.imagedict.keys()))+' rejected because of air')
        #remove aggregation
        target = min([self.imagedict[x]['lo_q'] for x in self.imagedict.keys() if not self.imagedict[x]['outlier']])
        for image in self.imagedict.keys():
            if self.imagedict[image]['lo_q'] > target*1.03:
                self.imagedict[image]['outlier'] = True
        self.logger.info(str(sum([self.imagedict[x]['outlier'] for x in self.imagedict.keys()]))+' of '+str(len(self.imagedict.keys()))+' rejected for air and rad damage')
        
    def Average(self):
        qdata = []
        idata = []
        edata = []
        n = sum([1 for x in self.imagedict.items() if not x[1]['outlier']])
        self.logger.info('n = '+str(n))
        for index, q in enumerate(self.imagedict[self.imagedict.keys()[0]]['Q']):

            i = np.array([ x[1]['I'][index] for x in self.imagedict.items() if not x[1]['outlier']]).mean()
            e = np.sqrt(np.array([ np.power(x[1]['E'][index],2) for x in self.imagedict.items() if not x[1]['outlier']]).sum()/n)
            
            qdata.append(q)
            idata.append(i)
            edata.append(e)
        self.logger.info('Averaged '+str(len(qdata))+' points')
        self.imagedict['outfile'] = {'Q': qdata, 'I': idata, 'E': edata}

    def Output(self):
        #string = 'Sum of images: '
        #for file in filelist:
        #    string = string+file+', '
        #print string[:-2]
        print "%-14s %-14s %-8s" % ("Q(A-1)","I(au)","Error")
        for q in self.imagedict['outfile']['Q']:
            index = self.imagedict['outfile']['Q'].index(q)
            q = "{0:11.9f}".format(q)
            i = "{0:15.9f}".format(self.imagedict['outfile']['I'][index])
            e = "{0:15.9f}".format(self.imagedict['outfile']['E'][index])
            print q+i+e


if __name__ == '__main__':
    job = SaxsAdd(filelist)
    job.TestFiles()
    job.ParseFiles()
    job.rejectOutliers()
    job.Average()
    job.Output()


