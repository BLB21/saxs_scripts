#!/usr/local/bin/python3
'''
Created on Aug 24, 2015

@author: nathan
'''

import re
import os
import sys
import logging
import numpy
from optparse import OptionParser
from optparse import OptionGroup
from rpy2.robjects.packages import importr
from subprocess import check_output

if len(sys.argv) < 2:
    sys.argv.append('-h')
    
class NMASaxs():
    """Make a set of pdb files by NMA and calculate fit to data
    
    The NMA tool takes a pdb file and uses the Bio3D package
    (http://thegrantlab.org/bio3d/tutorials/normal-mode-analysis)
    in the R statistics package (https://www.r-project.org/) to
    do normal mode analysis. The program then calls crysol or foxs
    to calculate chi fits to a supplied dat file either individually
    or as a set.

    """
    
    '''
    Constructor
    '''
    def __init__(self):
        ###start a log file
        self.logger = logging.getLogger('NMASaxs')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new NMASaxs job')        

        self.no_models = 10
        self.magnitude = 40
        self.output = 'nma_model.pdb'
        self.all_atom = True
        self.input_pdb = False
        self.allatom = False
        self.rtb = False
        self.mode = 7
        self.filelist_before = os.listdir(os.getcwd())
        self.unique_models = {}
        
    def SetInputPDB(self, input=False):
        if os.path.isfile(input) and input[-4:] == '.pdb':
            self.input_pdb = input
            self.logger.info('Input PDB file set to '+str(input))
        else:
            self.logger.error('Failed to set input pdb, does not exist or is not of type ".pdb"')

    def SetOutputPDB(self, output='nma_model.pdb'):
        if output[-4:] == '.pdb':
            self.output = output
            self.logger.info('Ouput PDB file set to '+output)
        else:
            self.logger.error('Failed to set output pdb, string must end in ".pdb"')

    def SetAllAtom(self, allatom=False):
        if type(allatom) == type(True):
            self.allatom = allatom
            if allatom:
                self.logger.info('Will create an all atom model')
            else:
                self.logger.info('Will create a C-alpha model')
        else:
            self.logger.error('Set allatom needs to be of type boolean, did not set')

    def SetRTB(self, rtb=False):
        if type(rtb) == type(True):
            self.rtb = rtb
            if rtb:
                self.logger.info('Will use RTB approximation. Does nothing for Calpha models')
            else:
                self.logger.info('Will not use RTB approximation.')
        else:
            self.logger.error('Set RTB needs to be of type boolean, did not set')

    def SetRange(self, mag=40):
        try:
            mag = int(mag)
            self.magnitude=mag
            self.logger.info('Range set to '+str(mag)+'A.')
        except:
            self.logger.error('Range must be an integer, unable to set')

    def SetMode(self, mode=7):
        try:
            mode = int(mode)
            self.mode=mode
            self.logger.info('Mode set to '+str(mode)+'.')
        except:
            self.logger.error('Mode must be an integer, unable to set')
            
    def SetNumberOfModels(self, models=10):
        try:
            models = int(models)
            self.no_models=models
            self.logger.info('Number of models set to '+str(models)+'.')
        except:
            self.logger.error('Number of models must be an integer, unable to set') 
       
    def CalculateNMA(self):
        if not self.input_pdb:
            self.logger.error('Cannot run CalculateNMA without an input PDB file')
            return False
        self.logger.info('Calculating NMA, this might take a few mins')
        bio3d = importr('bio3d')
        mypdb = bio3d.read_pdb(self.input_pdb)
        if self.allatom:
            self.logger.info('Running all atom nma')
            modes = bio3d.aanma(mypdb, outmodes='noh', rtb=self.rtb)
        else:
            modes = bio3d.nma(mypdb)
        trj = bio3d.mktrj(modes, mode=self.mode, mag=self.magnitude, step=float(self.magnitude)/float(self.no_models), file=self.output)
        self.logger.info('Trajecory models concatenated into file: '+self.output)

    def SplitOutputModel(self):
        #PARSE TRJ OUTPUT PDB FILE
        try:
            pdbdata = open(self.output, 'r').readlines()
        except:
            self.logger.error('Could not find the output file')
            return False
        all_models = {}
        model_number = 0
        for line in pdbdata:
            if line[0:5] == "MODEL":
                model_number += 1
                all_models[model_number] = []
            if line[0:4] == "ATOM":
                all_models[model_number].append(line)
        self.logger.info('Trajectory parsed into '+str(len(all_models.keys()))+' redundant models')
        self.unique_models[0] = all_models[1]
        window = int(round((len(all_models.keys())-2)/4))
        first_window = range(2,window+2)
        newno = 1
        for modelno in first_window:
            self.unique_models[newno] = all_models[modelno]
            newno += 1
        second_window = range((2*window)+3, (3*window)+3)
        newno = -1
        for modelno in second_window:
            self.unique_models[newno] = all_models[modelno]
            newno -= 1
        self.logger.info(str(len(self.unique_models))+' unique models')
        return True

    def RemoveTempFiles(self):
        self.logger.info('Deleting temporary files')
        filelist_after = os.listdir(os.getcwd())
        temp_files = list(set(filelist_after) - set(self.filelist_before))
        if len(temp_files) == 0:
            self.logger.info('No temporary files to delete')
        for file in temp_files:
            self.logger.info('Removing file: '+file)
            os.remove(file)

    def OutputMultipleFiles(self):
        self.logger.info('Outputting all models to multiple, single-model files')
        output_prefix = self.output[0:-4]+'_'
        print(output_prefix)
        for modelno in sorted(self.unique_models.keys()):
            with open(output_prefix+str(modelno).zfill(3)+'.pdb', 'w') as outfile:
                for line in self.unique_models[modelno]:
                    outfile.write(line)
                outfile.write('END\n')


    def OutputSingleFile(self):
        self.logger.info('Outputting all models to a single multi-model file')
        with open(self.output, 'w') as outfile:
            for modelno in sorted(self.unique_models.keys()):
                outfile.write('MODEL'+str(modelno).rjust(9)+'\n')
                for line in self.unique_models[modelno]:
                    outfile.write(line)
                outfile.write('ENDMDL\n')
            outfile.write('END\n')
                    
if __name__ == '__main__':
    '''
    Parse the command line options
    '''
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-p", "--pdb", action="store", type="string", dest="pdbfile", help="The pdb file you want to do the NMA on")

    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-o", "--out", action="store", type="string", dest="outfile", default="nma_model.pdb", help="The name of a pdb file you want to write out. (default is trajectory.pdb).")
    optional.add_option("-r", "--range", action="store", type="int", dest="range", default=40, help="The range of the movement along the major mode, (default 40)")
    optional.add_option("-m", "--mode", action="store", type="int", dest="mode", default=7, help="The mode you want to use, default 7, should be between 7 and maybe low teens.")
    optional.add_option("-n", "--number", action="store", type="int", dest="number", default=10, help="The number of models in the range, (default is 10)")
    optional.add_option("-c", "--clean", action="store_false", dest="clean", default=True, help="Remove the temporary files, this is a switch, no setting needed (default True)")
    optional.add_option("-a", "--allatom", action="store_true", dest="allatom", default=False, help="Calculate an all-atom model, alternative is a Calpha trace (default False)")
    optional.add_option("-f", "--fast", action="store_true", dest="rtb", default=False, help="Use RTB approximation. Speeds things up for all atom refinement. (default False)")
    optional.add_option("-t", "--trajectory", action="store_true", dest="trajectory", default=False, help="Output a single file with all models in it. (default False.)")
    optional.add_option("-e", "--eachmodel", action="store_true", dest="eachmodel", default=False, help="Output each model as a separate pdb file. (default False.)")            


    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    print(f'each model is {str(options.eachmodel)}')
    print(f'trajectory is {str(options.trajectory)}')
    '''
    Run the program
    '''
    job = NMASaxs()
    job.SetInputPDB(options.pdbfile)
    job.SetOutputPDB(options.outfile)
    job.SetMode(options.mode)
    job.SetNumberOfModels(options.number)
    job.SetRange(options.range)
    job.SetAllAtom(options.allatom)
    job.SetRTB(options.rtb)
    job.CalculateNMA()
    job.SplitOutputModel()    
    if options.clean:
        job.RemoveTempFiles()
    if options.eachmodel:
        job.OutputMultipleFiles()
    if options.trajectory:
        job.OutputSingleFile()
