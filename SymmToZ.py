#!/usr/local/bin/python
'''
Created on Sept 25, 2015

@author: nathan
'''


import logging
import math
import numpy
from os.path import isfile as isfile
from optparse import OptionParser
from optparse import OptionGroup
import sys
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d


sys.path.append('/Users/nathan/Documents/PYTHON')
from readwrite import PDB

if len(sys.argv) < 2:
    sys.argv.append('-h')
    
class SymmToZ(object):
    """Orients symmetry axis along the Z-fold axis
    
    Takes a pdb file with symmetry relating chains and orients chains to Z axis

    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('SymmToZ')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new SymmToZ job')        
        try:
            self._options = dict(options)
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()

        self.chains = []
        self.residues = []
        self.is_dna = False

    def ParsePDB(self):
        if isfile(self._options['pdbfile']) and self._options['pdbfile'][-4:] == '.pdb':

            self.logger.info('Reading and parsing '+self._options['pdbfile'])
            self.pdb = PDB(self._options['pdbfile'])
            self.pdb.parse_file()

    def GetChainsAndResn(self):
            #GET ALL CHAINS FROM PDB AND WHAT RESIDUE NUMBERS ARE IN THEM
            found_chains = {}
            for key in self.pdb.hashdata.keys():
                if self.pdb.hashdata[key]['record_type'] == 'ATOM':
                    if self.pdb.hashdata[key]['atom_name'] == 'CA' or self.pdb.hashdata[key]['atom_name'] == 'P':
                        if self.pdb.hashdata[key]['chain'] in found_chains.keys():
                            found_chains[self.pdb.hashdata[key]['chain']].append(self.pdb.hashdata[key]['residue_no'])

                        else:
                            found_chains[self.pdb.hashdata[key]['chain']] = [self.pdb.hashdata[key]['residue_no']]
                    else:pass
                else:pass

            #IF THE USER HAS SPECIFIED WHICH CHAINS TO USE THEN REMOVE CHAINS NOT IN THE LIST
            if self._options['chains'] != '*':
                try:
                    chains = self._options['chains'].split(',')
                    if len(chains) < 2:
                        sys.exit('You need to specify chains as a comma delimited list i.e. A,B and there needs to be at least 2')
                    else:pass
                except:
                    sys.exit('The chains option must be a string in the form A,B,C')
                for chain in found_chains.keys():
                    if not chain in chains:
                        found_chains.pop(chain)
            else:pass

            #GET A LIST OF COMMON RESIDUES FOUND IN ALL CHAINS
            self.common_residues = set.intersection(*map(set, list(found_chains.viewvalues())))
            self.chains = found_chains.keys()
            self.logger.info('Will use chains: '+','.join(sorted(self.chains)))
            self.logger.info('There are: '+str(len(self.common_residues))+' residues common to all chains')

            #RESTRICT TO JUST THE RESIDUES IN THE SPECIFIED RANGE IF A RANGE HAS BEEN SPECIFIED
            if self._options['numbers']:
                try:
                    from_residue = int(self._options['numbers'].split('-')[0])
                    to_residue = int(self._options['numbers'].split('-')[1])
                    restricted_list = []
                    for residue in self.common_residues:
                        if from_residue < residue < to_residue:
                            restricted_list.append(residue)
                    self.common_residues = restricted_list
                    self.logger.info('Restricted residue range from '+str(from_residue)+' to '+str(to_residue))
                except:
                    self.logger.error('Cannot parse residue range you tried to specify')

    def IsDNA(self):
        for key in self.pdb.hashdata.keys():
            if self.pdb.hashdata[key]['record_type'] == 'ATOM':
                if self.pdb.hashdata[key]['chain'] in self.chains:
                    if self.pdb.hashdata[key]['atom_name'] == 'P':
                        return True
        return False

    def FitLine(self):
        chain_sets = {}
        for chain in self.chains:
            chain_sets[chain] = []
        for key in self.pdb.hashdata.keys():
            if self.pdb.hashdata[key]['record_type'] == 'ATOM':
                if self.pdb.hashdata[key]['chain'] in self.chains:
                    if self.pdb.hashdata[key]['atom_name'] == 'CA' or self.pdb.hashdata[key]['atom_name'] == 'P':
                        if self.pdb.hashdata[key]['residue_no'] in self.common_residues:
                            chain_sets[self.pdb.hashdata[key]['chain']].append( (self.pdb.hashdata[key]['x'], self.pdb.hashdata[key]['y'], self.pdb.hashdata[key]['z']) )
        if self._options['reverse']:
            for index, chain in enumerate(sorted(chain_sets.keys())):
                if (index % 2 == 0):
                    self.logger.info('Reversing order of chain: '+chain)
                    chain_sets[chain] = list(reversed(chain_sets[chain]))
        resn_sets = []
        for index, coord in enumerate(chain_sets[chain_sets.keys()[0]]):
            coord_set = []
            for chain in chain_sets.keys():
                coord_set.append(chain_sets[chain][index])
            resn_sets.append(coord_set)

        line_points = []
        for item in resn_sets:
            x = numpy.array([i[0] for i in item]).mean()
            y = numpy.array([i[1] for i in item]).mean()
            z = numpy.array([i[2] for i in item]).mean()
            line_points.append([x, y, z])
        
        data = numpy.array(line_points)
        datamean = data.mean(axis=0)
        uu, dd, vv = numpy.linalg.svd(data - datamean)
        if not self._options['quiet']:
            ax = m3d.Axes3D(plt.figure())
            ax.scatter3D(*data.T)
            plt.show()
        else:
            self.logger.info('Scatter plot is off')

        #CENTRE ON ORIGIN
        self.pdb.Translate(str(datamean[0]/-1)+','+str(datamean[1]/-1)+','+str(datamean[2]/-1))
        
        #MAKE THE ROTATION MATRIX
        R = numpy.matrix('0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0')
        vector_from = vv[0]
        vector_from = vector_from / numpy.linalg.norm(vector_from)
        vector_to = numpy.array([0.0, 0.0, 1.0])#THE Z AXIS
        vector_to = vector_to / numpy.linalg.norm(vector_to)
        axis = numpy.cross(vector_from,vector_to)
        axis_len = numpy.linalg.norm(axis)
        axis = axis / axis_len
        x = axis[0]
        y = axis[1]
        z = axis[2]
        angle = math.acos(numpy.dot(vector_from, vector_to))
        ca = math.cos(angle)
        sa = math.sin(angle)

        R[0,0] = 1.0 + (1.0 - ca)*(x**2 - 1.0)
        R[0,1] = -z*sa + (1.0 - ca)*x*y
        R[0,2] = y*sa + (1.0 - ca)*x*z
        R[1,0] = z*sa+(1.0 - ca)*x*y
        R[1,1] = 1.0 + (1.0 - ca)*(y**2 - 1.0)
        R[1,2] = -x*sa+(1.0 - ca)*y*z
        R[2,0] = -y*sa+(1.0 - ca)*x*z
        R[2,1] = x*sa+(1.0 - ca)*y*z
        R[2,2] = 1.0 + (1.0 - ca)*(z**2 - 1.0)

        #APPLY THE ROTATION MATRIX
        for key in self.pdb.hashdata.keys():
            if self.pdb.hashdata[key]['record_type'] == 'ATOM':
                A = numpy.array([self.pdb.hashdata[key]['x'],self.pdb.hashdata[key]['y'],self.pdb.hashdata[key]['z']])
                A = numpy.dot(A, R.T)#Apply transpose of the rotation matrix
                self.pdb.hashdata[key]['x'] = A.item(0)
                self.pdb.hashdata[key]['y'] = A.item(1)
                self.pdb.hashdata[key]['z'] = A.item(2)

    def WriteFile(self):
        outfile = open(self._options['outfile'], 'w')
        outfile.write(self.pdb.return_file(justatoms=True))
        outfile.close()
        self.logger.info('written transformed pdb file to: '+self._options['outfile'])
        
if __name__ == '__main__':
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-p", "--pdbfile", action="store", type="string", dest="pdbfile", help="The pdb file you want to use") #A STRING
    required.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="The pdb file you want to write out") #A STRING

    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-c", "--chains", action="store", type="string", dest="chains", default="*", help="A comma delimited list of the chains you want to use i.e. A,B,C. Default is to use all chains") #A STRING
    optional.add_option("-n", "--numbers", action="store", type="string", dest="numbers", default=None, help="The residue range you want to use in the form 1-213, default is to use all.") #A STRING
    optional.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, help="don't show the scatter plot") #A BOOLEAN
    optional.add_option("-r", "--reverse", action="store_true", dest="reverse", default=False, help="reverse the order of atoms in every second chain (first res of chain A matches last res of chain B etc)") #A BOOLEAN    

    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()



    options = eval(str(options))
    job = SymmToZ(options)
    job.ParsePDB()
    job.GetChainsAndResn()
    job.FitLine()
    job.WriteFile()

