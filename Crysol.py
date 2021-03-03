#!/usr/local/bin/python

'''
Created on Apr 18, 2016

@author: nathan
'''





import glob
import logging
import matplotlib
from optparse import OptionParser
from optparse import OptionGroup
import os
import re
from subprocess import check_output
import sys

#import __main__
#__main__.pymol_argv = ['pymol', '-qc']
#sys.path.append('/usr/local/Cellar/pymol/1.7.2.1/lib/python2.7/site-packages/')
#import pymol

class Crysol():
    """Run crysol and harvest the results for one or many pdb files
    
    This function of this class is to be able to run crysol on many
    pdb files at the same time and collate the results into a report
    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('Crysol')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new crysol job')        
        try:
            self._options = dict(options)
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()

        ###CHECK PDB FILE(S)
        self._options['filelist'] = []
        filelist = glob.glob(self._options['pdbfile'])

        for file in filelist:
            if os.path.isfile(file) and file[-4:] == '.pdb':
                self._options['filelist'].append(file)
            else:
                self.logger.error("The file "+str(file)+" was in your list but either does not exist or is not a pdb file. It will be left out.")
        if len(self._options['filelist']) < 1:
            self.logger.error('There are no valid pdb files, will exit')
            sys.exit()
        else:
            self.logger.info('Will run crysol on '+str(len(self._options['filelist']))+' pdb files')


        ###CHECK DAT FILE
        if self._options['datfile'] == None:
            pass
        else:
            if os.path.isfile(self._options['datfile']) and self._options['datfile'][-4:] == '.dat':
                self.logger.info('Will compare pdb model with '+str(self._options['datfile']))
            else:
                self.logger.error("The dat file you specified either does not exist or is not of type '.dat'")
                sys.exit()

    def RunCrysol(self):
        self.chi_output = {}
        for pdbfile in self._options['filelist']:
            command = 'crysol '
            if self._options['subtract']:
                command = command + '-cst '
            command = command +'-ns '+str(self._options['no_points'])+' '+str(pdbfile)
            if not self._options['datfile'] == None:
                command = command+' '+str(self._options['datfile'])
            
            filelist_before = os.listdir(os.getcwd())
            output = check_output(command, shell=True)
            filelist_after = os.listdir(os.getcwd())
            crysol_files = list(set(filelist_after) - set(filelist_before))
            if not self._options['quiet']:
                fitfile = [ m for m in crysol_files if m[-4:] == '.fit'][0]
                with open(fitfile) as datfile:
                    filedata = datfile.readlines()
                datdata = {'Q': [], 'DATA': [], 'FIT': []}
                for line in filedata:
                    line = line.rsplit()
                    try:
                        datdata['Q'].append(float(line[0]))
                        datdata['DATA'].append(float(line[1]))
                        datdata['FIT'].append(float(line[2]))
                    except:
                        pass

                
            if not self._options['keep']:
                self.logger.info('Deleting crysol output files')
                for file in crysol_files:
                    os.remove(file)
            else:
                self.logger.info('Will keep crysol output files')
            output = output.split('\n')
            try:
                index = [i for i, item in enumerate(output) if re.search('.*Chi.*', item)][-1]
                self.chi_output[pdbfile] = "%.3f" % (float(output[index].split(':')[-1]))
            except:
                self.chi_output[pdbfile] = "****"
            

        for key in self.chi_output:
            print str(key)+','+self.chi_output[key]

if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.argv.append('-h')
        
    '''
    parse command line options
    '''
    
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-p", "--pdbfile", action="store", type="string", dest="pdbfile", help="The pdb file or files you want to calculate SAXS from. Can glob files i.e. '*.pdb'")

    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-d", "--datfile", action="store", type="string", dest="datfile", default=None, help="The dat file you want to compare your models with. Takes only one dat file.")
    optional.add_option("-k", "--keep", action="store_true", dest="keep", default=False, help="Keep the output files. Default is false i.e. delete them.")
    optional.add_option("-s", "--subtract", action="store_true", dest="subtract", default=False, help="Subtract constant, if flag is present will subtract, default is not to subtract.")
    optional.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, help="Do not plot the fit. Default is to plot.")
    #optional.add_option("-f", "--foxs", action="store_true", dest="foxs", default=False, help="Use foxs to calculate amplitudes. Default is false i.e. use crysol.")
    optional.add_option("-n", "--no_points", action="store", type="int", dest="no_points", default=51, help="The number of points in the output fit (default: 51)")

    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    if not options.pdbfile:
        sys.argv = [sys.argv[0], '-h']
        (options, args) = parser.parse_args()


    options = eval(str(options))
    job = Crysol(options)
    job.RunCrysol()
