#!/usr/local/bin/python
'''
Created on Jul 31, 2014

@author: nathan
'''

from operator import itemgetter
import sys
import logging
import os
import re
import copy
from optparse import OptionParser
from optparse import OptionGroup
from subprocess import check_output

import numpy

from mendeleev import element

class PDB():
    """Read and write pdb files
    
    The pdb class contains some functions for reading, writing and 
    performing calculations on pdb files.
    """
    
    '''
    Constructor
    '''
    def __init__(self, pdbfile):
        ###start a log file
        self.logger = logging.getLogger('PDB')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        if len(self.logger.handlers) < 1:
            self.logger.addHandler(streamhandler)
        while len(self.logger.handlers) > 1:
            self.logger.removeHandler(logger.handlers[0])

        self.logger.info('Starting a new PDB job')        

        if os.path.isfile(pdbfile) and pdbfile[-4:] == '.pdb':
            self.pdbfile = pdbfile
        else:
            sys.exit(str(pdbfile)+' either does not exist or is not of type ".pdb"')

        self.pdb_dict = {}

    def Read(self):
        self.logger.info('Reading and parsing pdb file: '+str(self.pdbfile))

        ######################################
        #PDB DEFINITION FROM wwPDB GUIDLINES #
        # VERSION 3.30 31/07/14              #
        ######################################
        self.pdb_definition = {
            'record_type': ( 0 , 6 ),
            'serial_no': ( 6 , 11 ),
            'atom_name': ( 12 , 16 ),
            'alternate': ( 16 , 17 ),
            'residue': ( 17 , 20 ),
            'chain': ( 21 , 22 ),
            'residue_no': ( 22 , 26 ),
            'icode': ( 26 , 27 ),
            'x': ( 30 , 38 ),
            'y': ( 38 , 46 ),
            'z': ( 46 , 54 ),
            'occupancy': ( 54 , 60 ),
            'bfactor': ( 60 , 66 ),
            'element': ( 76 , 78 ),
            'charge': ( 78 , 80 ),
            'string': ( 6 , '' )}
        self.integer_records = ['serial_no', 'residue_no']
        self.float_records = {'bfactor': 2, 'occupancy': 2, 'x': 3, 'y': 3, 'z': 3}

        file = open(self.pdbfile, 'r')
        pdblines = file.readlines()

        for line in pdblines:
            index = pdblines.index(line)
            self.pdb_dict[index] = {}
            record_type = line[self.pdb_definition['record_type'][0]:self.pdb_definition['record_type'][1]].rstrip()
            if record_type == 'ATOM' or record_type == 'HETATM':
                for record in self.pdb_definition.keys():
                    try:
                        if record in self.integer_records:
                            self.pdb_dict[index][record] = int(line[self.pdb_definition[record][0]:self.pdb_definition[record][1]].strip())
                        elif record in self.float_records.keys():
                            self.pdb_dict[index][record] = float(line[self.pdb_definition[record][0]:self.pdb_definition[record][1]].strip())
                        else:
                            self.pdb_dict[index][record] = line[self.pdb_definition[record][0]:self.pdb_definition[record][1]].strip()
                    except:
                        pass
            else:
                self.pdb_dict[index]['record_type'] = line[self.pdb_definition['record_type'][0]:self.pdb_definition['record_type'][1]].rstrip()
                self.pdb_dict[index]['string'] = line[self.pdb_definition['string'][0]:-1]+line[-1]

    def Write(self, justatoms=False):
        if justatoms:
            justatoms = True
            self.logger.info('Will only return the atom lines')

        output_lines = []
        #ljust_records = [ 'record_type', 'atom_name', 'string' ]
        ljust_records = [ 'record_type', 'string' ]
        
        for index in sorted(self.pdb_dict.keys()):
            try:
                if len(self.pdb_dict[index]['atom_name']) < 4:
                    self.pdb_dict[index]['atom_name'] = ' '+self.pdb_dict[option][index]['atom_name']
            except:
                pass
            line = list(' '*80)
            for record in self.pdb_dict[index]:
                try:
                    length = self.pdb_definition[record][1] - self.pdb_definition[record][0]
                except:
                    length = 80 - self.pdb_definition[record][0]
                if record in self.float_records:
                    if record in ljust_records:
                        content = list((('%.'+str(self.float_records[record])+'f') % self.pdb_dict[index][record]).ljust(length))
                    else:
                        content = list((('%.'+str(self.float_records[record])+'f') % self.pdb_dict[index][record]).rjust(length))
                        
                else:
                    if record in ljust_records:
                        content = list(str(self.pdb_dict[index][record]).ljust(length))
                    else:
                        if record == 'atom_name':
                            temp_string = str(self.pdb_dict[index][record]).ljust(3)
                        else:
                            temp_string = str(self.pdb_dict[index][record])
                        content = list(temp_string.rjust(length))
            
                try:
                    line[self.pdb_definition[record][0]:self.pdb_definition[record][1]] = content[0:]
                except:
                    line[self.pdb_definition[record][0]:80] = content[0:]
            if justatoms:
                if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                    output_lines.append(''.join(line))
            else:
                output_lines.append(''.join(line))


        return '\n'.join(output_lines)+'\n'


    def MolecularWeight(self):
        self.logger.info('Calculating molecular weight')
        total_mass = 0
        for index in self.pdb_dict.keys():
            if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                if self.pdb_dict[index]['element'] != '':
                    atom_name = self.pdb_dict[index]['element']
                else:
                    atom_name = self.pdb_dict[index]['atom_name'][0]

                try:
                    mass = element(atom_name).mass
                    total_mass += mass

                except:
                    self.logger.error('Failed to find mass for element '+str(atom_name))
        return total_mass
        
    def CentreOnOrigin(self):
        self.logger.info('Will move centre of mass to 0,0,0')

        #calculate centre of mass
        total_mass = 0
        total_x = 0
        total_y = 0
        total_z = 0
        for index in self.pdb_dict.keys():
            if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                if self.pdb_dict[index]['element'] != '':
                    atom_name = self.pdb_dict[index]['element']
                else:
                    atom_name = self.pdb_dict[index]['atom_name'][0]

                try:
                    mass = element(atom_name).mass
                    total_mass += mass
                    total_x += (mass * self.pdb_dict[index]['x'])
                    total_y += (mass * self.pdb_dict[index]['y'])
                    total_z += (mass * self.pdb_dict[index]['z'])

                except:
                    self.logger.error('Failed to find mass for element '+str(atom_name))

        self.logger.info('Total mass is {:.2f} KDa'.format(total_mass / 1000))
        average_x = total_x / total_mass
        average_y = total_y / total_mass
        average_z = total_z / total_mass
        self.logger.info('Centre of mass is currently: {0:.2f},{1:.2f},{2:.2f}'.format(average_x,average_y,average_z))
        for index in self.pdb_dict.keys():
            if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                self.pdb_dict[index]['x'] -= average_x
                self.pdb_dict[index]['y'] -= average_y
                self.pdb_dict[index]['z'] -= average_z
        self.logger.info('Moved molecule to 0,0,0')

    def Scale(self, factor):
        try:
            factor = float(factor)
        except:
            sys.exit('Scale factor must be an integer or float')

        for index in self.pdb_dict.keys():
           if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
               for axis in ['x','y','z']:
                   self.pdb_dict[index][axis] *= factor

    def DistanceBetween(self, first, second):
        try:
            first = int(first)
            second = int(second)
        except:
            self.logger.error('The DistanceBetween function requires two integers for the atom numbers')
            sys.exit()
        first_index = []
        second_index = []
        for index in self.pdb_dict.keys():
            try:
                if self.pdb_dict[index]['serial_no'] == first:
                    first_index.append(first)
                if self.pdb_dict[index]['serial_no'] == second:
                    second_index.append(second)
            except:
                pass

        if len(first_index) == 0 or len(second_index) == 0:
            self.logger.error('One of the atom numbers that you gave does not exist')
            sys.exit()
        if len(first_index) > 1:
            self.logger.error('There was more than one atom with the number '+str(first)+', will use first instance')
        if len(second_index) > 1:
            self.logger.error('There was more than one atom with the number '+str(second)+', will use first instance')

        distance = numpy.sqrt( numpy.power(self.pdb_dict[first_index[0]]['x'] - self.pdb_dict[second_index[0]]['x'], 2) +\
                               numpy.power(self.pdb_dict[first_index[0]]['y'] - self.pdb_dict[second_index[0]]['y'], 2) +\
                               numpy.power(self.pdb_dict[first_index[0]]['z'] - self.pdb_dict[second_index[0]]['z'], 2) )

        return distance

    def Rotate(self, option):
        if re.search('[xyzXYZ]', option[0]):
            try:
                axis = option[0].upper()
                theta = float(option[1:])
            except:
                self.logger.error("Rotation axis and angle must be in form 'x10' or 'Z2.5'")
                sys.exit()
        else:
            self.logger.error("Rotation axis and angle must be in form 'x10' or 'Z2.5'")
            sys.exit()
            
        def rotation_matrix(axis,theta):
            theta = numpy.radians(theta)
            if axis == 'X':
                return numpy.array([[1,0,0],[0,numpy.cos(theta),-numpy.sin(theta)],[0,numpy.sin(theta),numpy.cos(theta)]])
            elif axis == 'Y':
                return numpy.array([[numpy.cos(theta),0,numpy.sin(theta)],[0,1,0],[-numpy.sin(theta),0,numpy.cos(theta)]])
            elif axis == 'Z':
                return numpy.array([[numpy.cos(theta),-numpy.sin(theta),0],[numpy.sin(theta),numpy.cos(theta),0],[0,0,1]])
            else:
                sys.exit('error matrix function, required axis is neither x,y or z!')

        self.logger.info('Rotating around '+str(axis)+' by '+str(theta)+' degrees')
        for index in self.pdb_dict.keys():
            if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                v = numpy.array([self.pdb_dict[index]['x'],self.pdb_dict[index]['y'],self.pdb_dict[index]['z']])
                newcoords = (numpy.dot(rotation_matrix(axis, theta),v))
                self.pdb_dict[index]['x'] = float(newcoords[0])
                self.pdb_dict[index]['y'] = float(newcoords[1])
                self.pdb_dict[index]['z'] = float(newcoords[2])

    def Translate(self, option):
        translation = []
        for item in re.split('[,x ]', option):
            try:
                translation.append(float(item))
            except:
                pass
        if len(translation) != 3:
            self.logger.error("The translation matrix should be in the form i.e. '10.2,9.6,0.0'")
            sys.exit()
        self.logger.info('Translating by '+str(translation[0])+'x'+str(translation[1])+'x'+str(translation[2]))
        for index in self.pdb_dict.keys():
            if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                self.pdb_dict[index]['x'] = self.pdb_dict[index]['x'] + translation[0]
                self.pdb_dict[index]['y'] = self.pdb_dict[index]['y'] + translation[1]
                self.pdb_dict[index]['z'] = self.pdb_dict[index]['z'] + translation[2]
                              
    def ChainName(self, new):
        new = str(new).upper()
        if re.match('[A-Z]', new):
            self.logger.info('Changing chain name to "'+new+'"')
        else:
            self.logger.error('Chain name should be a single letter')
            sys.exit()

        for index in self.pdb_dict.keys():
            if self.pdb_dict[index]['record_type'] == 'ATOM' or self.pdb_dict[index]['record_type'] == 'HETATM':
                self.pdb_dict[index]['chain'] = new
                

if __name__ == '__main__':
    job = PDB(sys.argv[1])
    job.Read()
    #job.Scale(0.1)
    print job.Write()
