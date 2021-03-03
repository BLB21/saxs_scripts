#!/opt/anaconda3/bin/python

'''
Created on Apr 03, 2014

@author: nathan
'''


import sys
import os
import shutil
import random
import subprocess
from multiprocessing import Process as Process
from time import sleep as sleep
from time import time as time
import glob
import logging
from optparse import OptionParser
from optparse import OptionGroup
from subprocess import check_output

import __main__
__main__.pymol_argv = ['pymol', '-qc']
sys.path.append('/usr/local/Cellar/pymol/1.8.2.1/libexec/lib/python2.7/site-packages/')
import pymol

class DummyAtom():
    """Run a dummy atom modelling job on a given out file
    
    The DummyAtom class contains functions for running a dummy atom
    modelling job on a given 'out' file. You can specify the number
    of replicates or the program that you wish to use to do the
    modelling. The program will wait for the modelling jobs to finish
    and then run damaver after creating the appropriate directory
    structure.
    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('DummyAtom')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new dummy atom modelling job')        
        try:
            self._options = dict(options)
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()

        ###CHECK OUT FILE
        if os.path.isfile(self._options['file']) and self._options['file'][-4:] == '.out':
            logname = os.path.splitext(os.path.basename(self._options['file']))[0]
        else:
            self.logger.error("The file you specified either does not exist or is not of type '.out'")
            sys.exit()

        ###MAKE LOG FILE NAME IF NOT SPECIFIED
        if self._options['log'] == None:
            self._options['log'] = logname

        ###CHECK MODELLING PROGRAM IS SUPPORTED
        self.dammin = False
        self.dammif = False
        self.gasbor = False
        if self._options['program'] == 'dammif':
            self.dammif = True
        elif self._options['program'] == 'dammin':
            self.dammin = True
        elif self._options['program'] == 'gasbor':
            self.gasbor = True
        else:
            self.logger.error("Job must be either, 'dammin', 'dammif' or 'gasbor'.")
            sys.exit()

        ###IF GASBOR IS REQUIRED CHECK THAT NO. RESIDUES HAS BEEN SPECIFIED
        if self.gasbor and self._options['residues'] == 0:
            self.logger.error('To run gasbor you MUST provide the number of residues in the asymmertric particle')
            sys.exit()
        if self.gasbor and self._options['residues'] > 7999:
            self.logger.error('Gasbor cannot handle 8000 or more residues')
            sys.exit()

        ###CHECK SYMMETRY IS SUPPORTED
        supported = []
        for number in range(1,19+1):
            supported.append('P'+str(number))
        for number in range(2,12+1):
            supported.append('P'+str(number)+'2')
        supported.append('P23')
        supported.append('P432')
        supported.append('PICO')

        if not self._options['symmetry'] in supported:
            self.logger.error('Symmetry supported are: Point groups P1, ..., P19, Pn2 (n = 2, ..., 12), P23, P432 or PICO (icosahedral)')
            sys.exit()

        ###MAKE THE PROCESSING DIRECTORY
        index = 1
        dirname = self._options['program'].upper()
        while os.path.isdir(dirname):
            dirname = self._options['program'].upper()+'_'+str(index)
            index += 1
        self.logger.info('Making directory: '+dirname)
        os.mkdir(dirname)
        shutil.copyfile(self._options['file'], dirname+'/'+self._options['file'])


        ###CHECK THE PDB FILE TO SUPERIMPOSE
        if self._options['pdb'] != None:
            if os.path.isfile(self._options['pdb']) and self._options['pdb'][-4:] == '.pdb':
                self.logger.info('Will superimpose damfilt.pdb onto '+os.path.basename(self._options['pdb']))
                shutil.copyfile(self._options['pdb'], dirname+'/'+os.path.basename(self._options['pdb']))
            else:
                self.logger.error('The pdb file you specified either does not exist or is not of type .pdb')

        self.damaver_run = 0
        self.supcomb_run = 0
        os.chdir(dirname)
        
    def RunDammin(self):
        self.job_tracker = {}
        self.logger.info('Starting '+str(self._options['number'])+' '+self._options['program']+' jobs')
        index = 1
        while index < self._options['number'] + 1:
            self.job_tracker[index] = {}
            self.job_tracker[index]['command'] = 'dammin '+self._options['file']+\
                                                 ' /lo '+self._options['log'][:8]+'_'+str(index)+'.log'+\
                                                 ' /sy '+self._options['symmetry']+\
                                                 ' /SD '+str(random.randint(0,24000))
            def f():
                subprocess.call([self.job_tracker[index]['command']], shell=True, stdout=open(os.devnull, 'wb'))


            self.job_tracker[index]['process'] = Process(target=f)
            index += 1

        ###START ALL THE JOBS
        for index in self.job_tracker.keys():
            self.logger.info('Start job '+str(index))
            self.job_tracker[index]['process'].start()
            self.job_tracker[index]['starttime'] = time()

        ###WAIT FOR JOBS TO FINISH
        self.logger.info('Waiting for all jobs to finish')
        while True in {self.job_tracker[index]['process'].is_alive() for index in self.job_tracker.keys()}:
            sleep(5)
        self.logger.info('All processes have finished')

        ###CHECK EXIT STATUS
        for index in self.job_tracker.keys():
            if self.job_tracker[index]['process'].exitcode != 0:
                self.logger.error('Job number '+str(index)+' had an exit code of '+str(self.job_tracker[index]['process'].exitcode))


    def RunDammif(self):
        self.job_tracker = {}
        self.logger.info('Starting '+str(self._options['number'])+' '+self._options['program']+' jobs')
        index = 1
        while index < self._options['number'] + 1:
            self.job_tracker[index] = {}
            self.job_tracker[index]['command'] = 'dammif -p '+self._options['log'][:8]+'_'+str(index)+\
                                                 ' -m slow'+\
                                                 ' -s '+self._options['symmetry']+\
                                                 ' -q '+self._options['file']
            def f():
                subprocess.call([self.job_tracker[index]['command']], shell=True, stdout=open(os.devnull, 'wb'))

            self.job_tracker[index]['process'] = Process(target=f)
            index += 1
            

        ###START ALL THE JOBS
        for index in self.job_tracker.keys():
            self.logger.info('Start job '+str(index))
            self.job_tracker[index]['process'].start()
            self.job_tracker[index]['starttime'] = time()
            sleep(5)#YOU CAN'T SPECIFY A RANDOM SEED FOR DAMMIF SO SLEEP 5 GIVES RANDOM START


        ###WAIT FOR JOBS TO FINISH
        self.logger.info('Waiting for all jobs to finish')
        while True in {self.job_tracker[index]['process'].is_alive() for index in self.job_tracker.keys()}:
            sleep(5)
        self.logger.info('All processes have finished')

        ###CHECK EXIT STATUS
        for index in self.job_tracker.keys():
            if self.job_tracker[index]['process'].exitcode != 0:
                self.logger.error('Job number '+str(index)+' had an exit code of '+str(self.job_tracker[index]['process'].exitcode))

            

    def RunGasbor(self):
        self.job_tracker = {}
        self.logger.info('Starting '+str(self._options['number'])+' '+self._options['program']+' jobs')
        index = 1
        while index < self._options['number'] + 1:
            self.job_tracker[index] = {}
            self.job_tracker[index]['command'] = 'gasborp '+self._options['file']+\
                                                 ' '+str(self._options['residues'])+\
                                                 ' /lo '+self._options['log'][:4]+'_'+str(index)+'.log'+\
                                                 ' /sy '+self._options['symmetry']+\
                                                 ' /SD '+str(random.randint(0,24000))
            def f():
                subprocess.call([self.job_tracker[index]['command']], shell=True, stdout=open(os.devnull, 'wb'))

            self.job_tracker[index]['process'] = Process(target=f)
            index += 1

        ###START ALL THE JOBS
        for index in self.job_tracker.keys():
            self.logger.info('Start job '+str(index))
            self.job_tracker[index]['process'].start()
            self.job_tracker[index]['starttime'] = time()

        ###WAIT FOR JOBS TO FINISH
        self.logger.info('Waiting for all jobs to finish')
        while True in {self.job_tracker[index]['process'].is_alive() for index in self.job_tracker.keys()}:
            sleep(5)
        self.logger.info('All processes have finished')

        ###CHECK EXIT STATUS
        for index in self.job_tracker.keys():
            if self.job_tracker[index]['process'].exitcode != 0:
                self.logger.error('Job number '+str(index)+' had an exit code of '+str(self.job_tracker[index]['process'].exitcode))


    def RunDamaver(self):
        self.logger.info('Running Damaver')
            
        os.mkdir('DAMAVER')
        for file in glob.glob('*-1.pdb'):
            shutil.copyfile(file, 'DAMAVER/'+file)
        os.chdir('DAMAVER')
        
        command = 'damaver -a *.pdb'

        def f():
            subprocess.call([command], shell=True, stdout=open(os.devnull, 'wb'))

        p = Process(target=f)

        ###START THE JOB
        p.start()

        ###WAIT FOR COMPLETION AND CHECK EXIT CODE
        while p.is_alive():
            sleep(5)
        if p.exitcode != 0:
            self.logger.error('Damaver finished with a non-zero exit code, may be an error')
        else:
            self.damaver_run = 1
            self.logger.info('Damaver finished cleanly')

        self.logger.info('ALL DONE!')

    def SuperImpose(self):
        if self.damaver_run == 0 or self._options['pdb'] == None:
            return
        self.logger.info('Superimposing damfilt.pdb onto '+os.path.basename(self._options['pdb']))


        command = 'supcomb '+os.path.basename(self._options['pdb'])+' damfilt.pdb'

        def f():
            subprocess.call([command], shell=True, stdout=open(os.devnull, 'wb'))

        p = Process(target=f)

        ###START THE JOB
        p.start()

        ###WAIT FOR COMPLETION AND CHECK EXIT CODE
        while p.is_alive():
            sleep(5)
        if p.exitcode != 0:
            self.logger.error('Supcomb finished with a non-zero exit code, may be an error')
        else:
            self.supcomb_run = 1
            self.logger.info('Supcomb finished cleanly, superimposed molecule saved as damfiltr.pdb')
            
    def MakeImage(self):
        if self._options['image'] and self.supcomb_run:
            self.logger.info('Making an image')
        else:
            self.logger.info('Will not make an image')
            return

        pymol.finish_launching()
        pymol.cmd.bg_color(color='white')
        pymol.cmd.load('damfiltr.pdb', 'damfilt')
        pymol.cmd.orient('damfilt')
        pymol.cmd.disable('all')
        pymol.cmd.enable('damfilt')
        pymol.cmd.set('solvent_radius', 2)
        pymol.cmd.set('transparency', 0.5)
        pymol.cmd.set('surface_quality',2)
        pymol.cmd.alter('damfilt', 'vdw=3')
        pymol.cmd.hide('nonbonded')
        pymol.cmd.show('surface', 'damfilt')
        if self._options['pdb'] != None:
            pymol.cmd.set('transparency', 0.5)
            pymol.cmd.load(os.path.basename(self._options['pdb']), 'highres')
            pymol.cmd.enable('highres')
            pymol.cmd.show('cartoon', 'highres')
        pymol.cmd.ray(1500,1000,-1)
        pymol.cmd.png("damfilt_0deg.png", 1500, 1000, 300, 0)
        pymol.cmd.rotate('y', 90)
        pymol.cmd.ray(1500,1000,-1)
        pymol.cmd.png("damfilt_90deg.png", 1500, 1000, 300, 0)
        pymol.cmd.quit()

if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.argv.append('-h')
        
    '''
    parse command line options
    '''
    
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-f", "--file", action="store", type="string", dest="file", help="The out file you want to use for dummy atom modelling.")
    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-n", "--number", action="store", type="int", dest="number", default=10, help="The number of repeats you want to do (default = 10)")
    optional.add_option("-j", "--job", action="store", type="string", dest="program", default='dammif', help="The program you want to use, either 'dammin', 'dammif' or 'gasbor' (default is 'dammif')")
    optional.add_option("-l", "--log", action="store", type="string", dest="log", default=None, help="The prefix of the output files (default = the root part of the specified 'out' file)")
    optional.add_option("-s", "--symm", action="store", type="string", dest="symmetry", default="P1", help="The symmetry operation to apply i.e. P1, P2 etc (default = P1)")
    optional.add_option("-r", "--residues", action="store", type="int", dest="residues", default=0, help="The number of residues in the asymmetric part of the protein, only required for Gasbor")
    optional.add_option("-p", "--pdb", action="store", type="string", dest="pdb", default=None, help="A pdb file to superimpose the damfilt.pdb file onto.")
    optional.add_option("-i", "--image", action="store_false", dest="image", default=True, help="Default is to generate an image of the model, calling this option with no argument will turn this off.")


    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    

    options = eval(str(options))
    job = DummyAtom(options)
    if job.dammin:
        job.RunDammin()
    elif job.dammif:
        job.RunDammif()
    elif job.gasbor:
        job.RunGasbor()
    else:
        pass
    job.RunDamaver()
    job.SuperImpose()
                                          

