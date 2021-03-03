#!/usr/local/bin/python
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
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('NMASaxs')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new NMASaxs job')        
        try:
            self._options = dict(options)
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()
        ###set the default value for step if it has not been set
        if self._options['step'] == 0:
            self._options['step'] = int(round(self._options['magnitude'] / 10))
        ###get a snapshot of the files before for use in cleanup afterwards
        self.filelist_before = os.listdir(os.getcwd())

    def CalculateNMA(self):
        self.logger.info('Running Normal Mode Analysis')
        self.output = 'mag'+str(self._options['magnitude'])+'.pdb'
        bio3d = importr('bio3d')
        mypdb = bio3d.read_pdb(self._options['pdbfile'])
        modes = bio3d.nma(mypdb)
        trj = bio3d.mktrj(modes, mode=self._options['mode'], mag=self._options['magnitude'], step=self._options['step'], file=self.output)
        self.logger.info('Trajecory models concatenated into file: '+self.output)
        #PARSE TRJ OUTPUT PDB FILE
        pdbdata = open(self.output, 'r').readlines()
        all_models = {}
        model_number = 0
        for line in pdbdata:
            if line[0:5] == "MODEL":
                model_number += 1
                all_models[model_number] = []
            if line[0:4] == "ATOM":
                all_models[model_number].append(line)
        self.logger.info('Trajectory parsed into '+str(len(all_models.keys()))+' redundant models')
        self.unique_models = {}
        self.unique_models[0] = all_models[1]
        window = (len(all_models.keys())-2)/4
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
        

    def CrysolOnSingles(self):
        if self._options['individual']:
            self.logger.info('Running Crysol on single models')
            self.unique_dats = {}
            for modelno in sorted(self.unique_models.keys()):
                outfile_name = 'MODEL_'+str(modelno).zfill(3)+'.pdb'
                with open(outfile_name, 'w') as outfile:
                    outfile.write('\n'.join(map(str.strip, self.unique_models[modelno])))

                command = 'crysol '+str(outfile_name)+' '+str(self._options['datfile'])
                self.logger.info('Running: '+command)
                output = check_output(command, shell=True)
                output = output.split('\n')

                try:
                    index = [i for i, item in enumerate(output) if re.search('.*Chi.*', item)][-1]
                    chi_string = "%.2f".rjust(5) % (float(output[index].split(':')[-1]))
                except:
                    self.logger.error('Could not resolve chi as an integer, defaulting to 50')
                    chi_string = "%.2f".rjust(5) % 50.0
                self.logger.info(str(modelno)+': '+chi_string)
                replacement_list = []
                for line in self.unique_models[modelno]:
                    replacement_list.append(line[0:60]+chi_string+line[66:])
                self.unique_models[modelno] = replacement_list

    def CrysolOnMultiples(self):
        if not self._options['individual']:
            self.logger.info('Running Crysol and averaging all models')

            ###Parse the input dat file
            datfile = open(self._options['datfile'], 'r').readlines()
            self.datdict = {'Q': [], 'I': [], 'E': []}
            for line in datfile:
                line = line.split()
                try:
                    self.datdict['Q'].append(float(line[0]))
                    self.datdict['I'].append(float(line[1]))
                    self.datdict['E'].append(float(line[2]))
                except:
                    pass
            ###Run Crysol on all the individuals and parse the outputs
            self.unique_dats = {}
            for modelno in sorted(self.unique_models.keys()):
                outfile_name = 'MODEL_'+str(modelno).zfill(3)+'.pdb'
                with open(outfile_name, 'w') as outfile:
                    outfile.write('\n'.join(map(str.strip, self.unique_models[modelno])))
    
                command = 'crysol '+str(outfile_name)+' '+str(self._options['datfile'])
                self.logger.info('Running: '+command)
                output = check_output(command, shell=True)
                output = output.split('\n')
                fitfile_name = [ line.split()[-1].rstrip() for line in output if re.match(re.compile('.*saved to file.*'), line) ][-1]
                fitdata = open(fitfile_name, 'r').readlines()
                qs = []
                obss = []
                mods = []
                self.unique_dats[outfile_name] = {'Q': [], 'OBS': [], 'MOD': [], 'ERR': []} 
                for line in fitdata:
                    line = line.split()
                    try:
                        if float(line[0]) in self.datdict['Q']:
                            index = self.datdict['Q'].index(float(line[0]))
                            self.unique_dats[outfile_name]['Q'].append(float(line[0]))
                            self.unique_dats[outfile_name]['OBS'].append(float(line[1]))
                            self.unique_dats[outfile_name]['MOD'].append(float(line[2]))
                            self.unique_dats[outfile_name]['ERR'].append(self.datdict['E'][index])
                    except:
                        pass
                self.logger.info('Parsed '+fitfile_name+' with '+str(len(self.unique_dats[outfile_name]['Q']))+' points')
    
            ###Average them all together
            self.averaged_dat = {'Q': [], 'OBS': [], 'MOD': [], 'ERR': []}
            for q in self.unique_dats[self.unique_dats.keys()[0]]['Q']:
                index = self.unique_dats[self.unique_dats.keys()[0]]['Q'].index(q)
                total_obs = []
                total_mod = []
                total_err = []
                for model in self.unique_dats.keys():
                    total_obs.append(self.unique_dats[model]['OBS'][index])
                    total_mod.append(self.unique_dats[model]['MOD'][index])
                    total_err.append(self.unique_dats[model]['ERR'][index])
                self.averaged_dat['Q'].append(q)
                self.averaged_dat['OBS'].append(numpy.mean(total_obs))
                self.averaged_dat['MOD'].append(numpy.mean(total_mod))
                self.averaged_dat['ERR'].append(numpy.mean(total_err))

            ###Calculate Final Chi
            total = 0
            #WORK OUT THE SCALE FACTOR
            firstterm = 0
            secondterm = 0
            for q in self.averaged_dat['Q']:
                index = self.averaged_dat['Q'].index(q)
                firstterm += ( ( self.averaged_dat['OBS'][index] * self.averaged_dat['MOD'][index] ) / self.averaged_dat['ERR'][index]**2 )
                secondterm += ( self.averaged_dat['MOD'][index]**2 / self.averaged_dat['ERR'][index]**2 )
            scalefactor = firstterm / secondterm
            self.logger.info('Scalefactor is '+str(scalefactor))

            for q in self.averaged_dat['Q']:
                index = self.averaged_dat['Q'].index(q)
                total += ( ( self.averaged_dat['OBS'][index] - scalefactor * self.averaged_dat['MOD'][index] ) / self.averaged_dat['ERR'][index] )
            final_chi = ( total )**2 / len(self.averaged_dat['Q'] )
            self.logger.info('Chi score for average is '+str(final_chi))
            print self.averaged_dat


                                         


    def RemoveTempFiles(self):
        filelist_after = os.listdir(os.getcwd())
        temp_files = list(set(filelist_after) - set(self.filelist_before))
        if self._options['clean']:
            self.logger.info('Deleting temporary files')
            for file in temp_files:
                os.remove(file)
        else:
            self.logger.info('Will leave temporary files')

    def OutputSingleFile(self):
        if self._options['trajectory']:
            self.logger.info('Outputting all models to one multi-model file')
            with open(self._options['outfile'], 'w') as outfile:
                for modelno in sorted(self.unique_models.keys()):
                    outfile.write('MODEL'+str(modelno).rjust(9)+'\n')
                    for line in self.unique_models[modelno]:
                        outfile.write(line)
                    outfile.write('ENDMDL\n')

    def OutputMultipleFiles(self):
        if not self._options['trajectory']:
            self.logger.info('Outputting all models to multiple, single-model files')
            for modelno in sorted(self.unique_models.keys()):
                with open('MODEL_'+str(modelno).zfill(3)+'.pdb', 'w') as outfile:
                    for line in self.unique_models[modelno]:
                        outfile.write(line)
                    outfile.write('END\n')
                              
if __name__ == '__main__':
    '''
    Parse the command line options
    '''
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-p", "--pdb", action="store", type="string", dest="pdbfile", help="The pdb file you want to do the NMA on")

    optional = OptionGroup(parser, "Optional Arguments")
    required.add_option("-o", "--out", action="store", type="string", dest="outfile", default="trajectory.pdb", help="The name of a pdb file you want to write out. (default is trajectory.pdb).")
    optional.add_option("-d", "--dat", action="store", type="string", dest="datfile", help="The dat file you want to do compare the NMA with")
    optional.add_option("-m", "--magnitude", action="store", type="int", dest="magnitude", default=40, help="The magnitude of the movement along the major mode, (default 40)")
    optional.add_option("-n", "--normalmode", action="store", type="int", dest="mode", default=7, help="The mode you want to use, default 7, should be between 7 and maybe low teens.")
    optional.add_option("-s", "--step", action="store", type="int", dest="step", default=0, help="The step size, should be much less than magnitude, (default is magnitude/10)")
    optional.add_option("-c", "--clean", action="store_false", dest="clean", default=True, help="Remove the temporary files, this is a switch, no setting needed (default True)")
    optional.add_option("-i", "--individual", action="store_true", dest="individual", default=False, help="Run crysol on each structure separately and report chi, default is to treat the files as a set and average")
    optional.add_option("-t", "--trajectory", action="store_false", dest="trajectory", default=False, help="Save output pdb files as a single multimodel file. (default is False).")


    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    options = eval(str(options))

    '''
    Run the program
    '''
    job = NMASaxs(options)
    job.CalculateNMA()
    job.CrysolOnSingles()
    job.CrysolOnMultiples()
    job.RemoveTempFiles()    
    job.OutputSingleFile()
    job.OutputMultipleFiles()
