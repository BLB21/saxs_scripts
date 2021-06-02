#!/opt/anaconda3/bin/python

'''
Created on Apr 18, 2016

@author: nathan
'''





import glob
import logging
from optparse import OptionParser
from optparse import OptionGroup
import os
import re
from subprocess import check_output, STDOUT
import sys
import matplotlib.pyplot as plt
#import __main__
#__main__.pymol_argv = ['pymol', '-qc']
#sys.path.append('/usr/local/Cellar/pymol/1.7.2.1/lib/python2.7/site-packages/')
#import pymol

class SaxsCalc():
    """Run crysol or Foxs and plot the results
    
    This function of this class is to be able to run crysol or FoXs on
    a pdb file at and collate the results into a report
    """
    
    '''
    Constructor
    '''
    __version__ = '1.05b'
    def __init__(self):
        ###start a log file
        self.logger = logging.getLogger('SaxsCalc')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new SaxsCalc job')        

        ####SET SOME PARAMETERS
        self.datdata = {'Q': [], 'DATA': [], 'FIT': []}
        self.datfile = None
        self.datfile_name = None
        self.pdbfile = None
        self.clean_files = True
        self.number_of_points = 51
        self.subtract_constant = False
        self.chi_score = None
        self.created_files = []
        
    def AddPdbFile(self, pdbfile):
        if os.path.isfile(pdbfile) and pdbfile[-4:] == '.pdb':
            self.logger.info('Added pdb file '+str(pdbfile))
            self.pdbfile = pdbfile
            return True
        else:
            self.logger.error('You did not enter a valid pdb file')
            return False

    def ReturnPdbFile(self):
        return self.pdbfile

    def ReturnDatFileName(self):
        return self.datfile_name
    
    def AddDatFile(self, datfile):
        if os.path.isfile(datfile) and datfile[-4:] == '.dat':
            self.logger.info('Added dat file '+str(datfile))
            self.datfile = datfile
            self.datfile_name = datfile
            return True
        else:
            self.logger.error('You did not enter a valid dat file')
            return False

    def AutoLimitDatFile(self):
        self.logger.info('cutting out dodgy points from low Q')
        first_good_point = 1
        #test autorg is on the path
        command = 'which autorg'
        output = check_output(command, shell=True).rstrip()
        if not os.path.isfile(output):
            sys.exit('autorg is not on your path, use the -a option to turn off points truncation')
        #parse the dat data
        dat_data = {'Q': [], 'I': [], 'E': []}
        if self.datfile:
            with open(self.datfile) as f:
                filedata = f.readlines()
            for line in filedata:
                line = line.rsplit()
                try:
                    dat_data['Q'].append(float(line[0]))
                    dat_data['I'].append(float(line[1]))
                    dat_data['E'].append(float(line[2]))
                except:
                    pass
            #find out where its good from
            command = 'autorg '+self.datfile
            output = check_output(command, shell=True).split('\n'.encode('utf-8'))
            for line in output:
                if line[:6] == 'Points':
                    try:
                        first_good_point = int(line.split()[1]) - 1
                        self.logger.info(f'First {first_good_point-1} will be removed')
                    except:
                        self.logger.error('could not determine good points from autorg')
            #write out the good data
            if first_good_point > 1:
                self.logger.info('will not include points below: '+str(first_good_point))
                outstring_array = []
                for index, q in enumerate(dat_data['Q']):
                    outstring_array.append('{0: <16.9f}{1: <16.9f}{2: <16.9f}'.format(
                                q,
                                dat_data['I'][index],
                                dat_data['E'][index]))
                outfile_name = 'saxscalc.dat'
                outfile = open(outfile_name, 'w')
                outfile.write('\n'.join(outstring_array[first_good_point:]))
                outfile.close()
                self.datfile = outfile_name
                self.files_for_deletion.append(outfile_name)
            
    def SetNumberOfPoints(self, number_of_points=51):
        if type(number_of_points) == type(1):
            self.logger.info('Number of points set to: '+str(number_of_points))
            self.number_of_points = number_of_points
            return True
        else:
            self.logger.error('SetNumberOfPoints needs an integer as an argument')
            return False

    def SetCleanFiles(self, clean_files=True):
        if type(clean_files) == type(True):
            if clean_files:
                self.logger.info('Will remove files after running')
            else:
                self.logger.info('Will not remove files after running')
            self.clean_files = clean_files
            return True
        else:
            self.logger.error('SetCleanFiles needs a boolean (True or False)')
            return False

    def SetSubtractConstant(self, subtract_constant=True):
        if type(subtract_constant) == type(True):
            if subtract_constant:
                self.logger.info('Will subtract a constant')
            else:
                self.logger.info('Will not subtract a constant')
            self.subtract_constant = subtract_constant
            return True
        else:
            self.logger.error('SetSubtractConstant needs a boolean (True or False)')
            return False

    def OutputChiScore(self):
        self.logger.info('Chi score is: '+str(self.chi_score))
        return self.chi_score

    def DeleteFiles(self):
        self.logger.info('Deleting output files')
        for file in self.created_files:
            os.remove(file)

    def PlotTheFit(self, outfile=False):
        self.logger.info('Plotting the data')
        fig = plt.figure(figsize=(10,6), dpi= 100, facecolor='w', edgecolor='k')
        ax1 = fig.add_subplot(111, autoscale_on=True)
        ax1.plot(self.datdata['Q'], self.datdata['DATA'], marker='.', color='0.55', markersize = 4.0, linestyle='None')
        ax1.plot(self.datdata['Q'], self.datdata['FIT'], marker='None', color='0', linewidth=2.0, linestyle='-')
        ax1.set_xlabel('$Q (\AA^{-1}$)', fontsize=20, color='black')
        ax1.set_ylabel('$I(0)$', fontsize=20, color='black')
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.autoscale(enable=True, axis='both', tight=None)
        if outfile:
            plt.savefig(outfile)
        plt.show()
        
    def RunFoxs(self):
        self.logger.info('Running FoXS')
        #test foxs is on the path
        command = 'which foxs'
        output = check_output(command, shell=True).rstrip()
        if not os.path.isfile(output):
            sys.exit('foxs is not on your path, try the -j crysol option')

        command = 'foxs -s '+str(self.number_of_points)+' '
        
        if self.subtract_constant:
            command = command + '-o '

        if not self.pdbfile:
            self.logger.error('Cannot execute run foxs without a pdb file')
            return False
        else:
            command = command + str(self.pdbfile)
            os.chdir(os.path.split(os.path.abspath(self.pdbfile))[0])

        if self.datfile:
            command = command +' '+str(self.datfile)

        filelist_before = os.listdir(os.getcwd())
        output = check_output(command, shell=True, stderr=STDOUT)

        filelist_after = os.listdir(os.getcwd())
        self.created_files = self.created_files + list(set(filelist_after) - set(filelist_before))

        output = output.split('\n'.encode('utf-8'))
        self.output = output
        pattern = re.compile('.* Chi\^2 = .*'.encode('utf-8'))
        for line in output:
            if re.match(pattern, line):
                line = line.split()
                try:
                    self.chi_score = "%.3f" % (float(line[line.index('='.encode('utf-8'))+1]))
                except:
                    self.logger.error('Could not find chi in the foxs output')
                    self.chi_score = "****"
        if self.datfile:
            #PARSE THE FIT FILE
            pattern = re.compile('.*fit')
            fitfile = False
            for f in self.created_files:
                if re.match(pattern, f):
                    fitfile = f
            if not fitfile:
                fs = sorted(os.listdir(os.getcwd()), key=os.path.getctime)
                fs.reverse()
                for f in fs:
                    if re.match(pattern, f):
                        fitfile = f
                        break
            if fitfile:
                with open(fitfile) as f:
                    filedata = f.readlines()
                for line in filedata:
                    line = line.rsplit()
                    try:
                        self.datdata['Q'].append(float(line[0]))
                        self.datdata['DATA'].append(float(line[1]))
                        self.datdata['FIT'].append(float(line[3]))
                    except:
                        pass
                
        #TEST THE OUTPUT
        if len(self.datdata['Q']) > 0:
            return True
        else:
            return False


        
    def RunCrysol(self):
        self.logger.info('Running Crysol')
        #test crysol is on the path
        command = 'which crysol'
        output = check_output(command, shell=True).rstrip()
        if not os.path.isfile(output):
            sys.exit('crysol is not on your path, try the -j foxs option')

        command = 'crysol -ns '+str(self.number_of_points)+' '
        
        if self.subtract_constant:
            command = command + '-cst '

        if not self.pdbfile:
            self.logger.error('Cannot execute runCrysol without a pdb file')
            return False
        else:
            command = command + str(self.pdbfile)

        if self.datfile:
            command = command +' '+str(self.datfile)
        
        filelist_before = os.listdir(os.getcwd())
        output = check_output(command, shell=True)
        filelist_after = os.listdir(os.getcwd())
        self.created_files = self.created_files + list(set(filelist_after) - set(filelist_before))

        if self.datfile:
            #PARSE THE CRYSOL OUTPUT
            output = output.split('\n'.encode('utf-8'))
            self.output = output            
            try:
                index = [i for i, item in enumerate(output) if re.search('.*Chi.*'.encode('utf-8'), item)][-1]
                self.chi_score = "%.3f" % (float(output[index].split(':'.encode('utf-8'))[-1]))
            except:
                self.logger.error('Could not find chi in the crysol output')
                self.chi_score = "****"
    
            #PARSE THE FIT FILE
            fitfile = [ m for m in self.created_files if m[-4:] == '.fit'][0]
            with open(fitfile) as f:
                filedata = f.readlines()
            for line in filedata:
                line = line.rsplit()
                try:
                    self.datdata['Q'].append(float(line[0]))
                    self.datdata['DATA'].append(float(line[1]))
                    self.datdata['FIT'].append(float(line[2]))
                except:
                    pass
    
            #REMOVE EXTRAPOLATED POINTS AT START OF FIT
            while self.datdata['DATA'][0] == self.datdata['DATA'][1]:
                del self.datdata['Q'][0]
                del self.datdata['DATA'][0]
                del self.datdata['FIT'][0]
            del self.datdata['Q'][0]
            del self.datdata['DATA'][0]
            del self.datdata['FIT'][0]
                
            #TEST THE OUTPUT
            if len(self.datdata['Q']) > 0:
                return True
            else:
                return False


if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.argv.append('-h')
        
    '''
    parse command line options
    '''
    
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-p", "--pdbfile", action="store", type="string", dest="pdbfile", help="The pdb file or files you want to calculate SAXS from. Can glob the list but put your pattern in single quotes i.e. ...-p '*.pdb'")

    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-d", "--datfile", action="store", type="string", dest="datfile", default=None, help="The dat file you want to compare your models with. Takes only one dat file.")
    optional.add_option("-o", "--outfile", action="store", type="string", dest="outfile", default=None, help="The image file of the fit you want to save. Can be .png or .pdf, default is None.")
    optional.add_option("-k", "--keep", action="store_true", dest="keep", default=False, help="Keep the output files. Default is false i.e. delete them.")
    optional.add_option("-s", "--subtract", action="store_true", dest="subtract", default=False, help="Subtract constant, if flag is present will subtract, default is not to subtract.")
    optional.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, help="Do not plot the fit. Default is to plot.")
    optional.add_option("-a", "--all_points", action="store_true", dest="all_points", default=False, help="Use all the points, default is to cut low Q points based on autorg.")
    optional.add_option("-w", "--water", action="store_true", dest="all_points", default=False, help="Do not include hydration layer, default is hydration.")    
    optional.add_option("-j", "--job", action="store", type="string", dest="job", default='crysol', help="The calculator to use, options can be crysol or foxs, default is crysol")
    optional.add_option("-n", "--no_points", action="store", type="int", dest="no_points", default=51, help="The number of points in the output fit (default: 51)")

    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    if not options.pdbfile:
        sys.argv = [sys.argv[0], '-h']
        (options, args) = parser.parse_args()

    for pdbfile in glob.glob(options.pdbfile):
        job = SaxsCalc()
        job.AddPdbFile(pdbfile)
        if options.datfile:
            job.AddDatFile(options.datfile)
        job.SetSubtractConstant(options.subtract)
        job.SetNumberOfPoints(options.no_points)
        if not options.all_points:
            job.AutoLimitDatFile()
        if options.job == 'foxs':
            job.RunFoxs()
        else:
            job.RunCrysol()
        if not options.keep:
            job.DeleteFiles()
        if options.datfile:
            if not options.quiet:
                job.PlotTheFit(options.outfile)
            else:
                job.logger.info('No plot requested')
            print(f'{job.ReturnPdbFile()} vs {job.ReturnDatFileName()}: {job.OutputChiScore()}')
