#!/opt/anaconda3/bin/python
'''
Created on Jul 31, 2014

@author: nathan
'''
from copy import deepcopy
import sys
import logging
import os
import re
import numpy
import periodictable

### A common interface contract for all readwrite classes
class Interface(object):    
    def parse_file(self): raise RuntimeError('Not implemented')
    def return_file(self): raise RuntimeError('Not implemented')
    def input_dict(self): raise RuntimeError('Not implemented')
    def return_dict(self): raise RuntimeError('Not implemented')


class PDB(Interface):
    """Read and write pdb files
    
    The pdb class contains functions for reading pdb files to
    a dictionary object where various calculations can be
    performed on it and then writing it back out again.
    """
    def __init__(self, pdbfile=None):
        ###start a log file
        self.logger = logging.getLogger('readwrite.PDB')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(streamhandler)
            self.logger.info('Starting a new '+str('readwrite.PDB')+' job')        
        ###DEFINE SOME PARAMETERS
        self.type = 'pdb'
        self.centre = (0,0,0)
        self.integer_records = ['serial_no', 'residue_no']
        self.float_records = {'bfactor': 2, 'occupancy': 2, 'x': 3, 'y': 3, 'z': 3}
        self.converter = {
            'GLY': ('Glycine','G', 57.05),
            'PRO': ('Proline','P', 97.12),
            'ALA': ('Alanine','A', 71.09),
            'VAL': ('Valine','V', 99.14),
            'LEU': ('Leucine','L', 113.16),
            'ILE': ('Isoleucine','I', 113.16),
            'MET': ('Methionine','M', 131.19),
            'CYS': ('Cysteine','C', 103.15),
            'PHE': ('Phenylalanine', 'F', 147.18),
            'TYR': ('Tyrosine','Y', 163.18),
            'TRP': ('Tryptophan','W', 186.12),
            'HIS': ('Histidine','H', 137.14),
            'LYS': ('Lysine','K', 128.17),
            'ARG': ('Arginine','R', 156.19),
            'GLN': ('Glutamine','Q', 128.14),
            'ASN': ('Asparagine','N', 115.09),
            'GLU': ('Glutamic Acid','E', 129.12),
            'ASP': ('Aspartic Acid','D', 114.11),
            'SER': ('Serine','S', 87.08),
            'THR': ('Threonine','T', 101.11)
        }
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
        
        self.new_atom_values = {
            'record_type': 'ATOM',
            'serial_no': 1,
            'atom_name': 'CA',
            'alternate': '',
            'residue': 'ALA',
            'chain': 'A',
            'residue_no': 1,
            'icode': '',
            'x': 0.0,
            'y': 0.0,
            'z': 0.0,
            'occupancy': 1.0,
            'bfactor': 20.0,
            'element': 'C',
            'charge': ''
        }


        self.hashdata = {}

        ###Check file exists and is of right type
        if pdbfile == None:
            self.pdbfile = None
            self.logger.info('Will create an empty pdb dictionary')

        elif os.path.isfile(pdbfile) and pdbfile[-4:] == '.pdb':
            self.pdbfile = pdbfile
        else:
            sys.exit(str(pdbfile)+' either does not exist or is not of type ".pdb"')
            
    def input_dict(self, input_dict):
        self.logger.info('Reading in PDB data as a dictionary')
        if type(input_dict) == type({}):
            self.hashdata = input_dict

    def return_dict(self):
        self.logger.info('Returning PDB data as a dictionary')
        return self.hashdata
    
    def addNewAtom(self, coordinates=None, bfactor=None):
        if type(coordinates) in [type(()),type([])]:
            if len(coordinates) == 3:
                self.setNewAtomProperty('x', coordinates[0])
                self.setNewAtomProperty('y', coordinates[1])
                self.setNewAtomProperty('z', coordinates[2])
        if bfactor:
            self.setNewAtomProperty('bfactor', bfactor)
        if len(self.hashdata.keys()) == 0:
               index = 1
        else:
               index = sorted(self.hashdata.keys())[-1]+1
        self.hashdata[index] = deepcopy(self.new_atom_values)
        
    def setNewAtomProperty(self, property=None, value=None):
        if property in self.new_atom_values.keys():
            try:
                if property in self.integer_records:
                    self.new_atom_values[property] = int(value)
                elif property in self.float_records:
                    self.new_atom_values[property] = float(value)
                else:
                    self.new_atom_values[property] = str(value)
                self.logger.info('Set: '+str(property)+' to: '+str(value))
            except:
                self.logger.error('Could not set: '+str(property)+' to: '+str(value))
        else:
            self.logger.error('New atoms have no property: '+str(property))
            
    def CACObyNumber(self, residue_number=None):
        return_value = {}
        for index in self.hashdata.keys():
            if self.hashdata[index]['record_type'] == 'ATOM':
                if self.hashdata[index]['residue_no'] == residue_number and self.hashdata[index]['atom_name'] in ['CA','C','O']:
                    return_value[self.hashdata[index]['atom_name']] = (self.hashdata[index]['x'],self.hashdata[index]['y'],self.hashdata[index]['z'])
        if len(return_value.keys()) == 3:
            return return_value
        else:
            return False
        
    def ReturnSeq(self, code=1):
        seq = []

        for index in self.hashdata.keys():
            if self.hashdata[index]['record_type'] == 'ATOM' and self.hashdata[index]['atom_name'] == 'CA':
                if self.hashdata[index]['residue'] in self.converter.keys():
                    if code == 1:
                        seq.append(self.converter[self.hashdata[index]['residue']][1])
                    else:
                        seq.append(self.hashdata[index]['residue'])
                else:
                    if code == 1:
                        seq.append('X')
                    else:
                        seq.append('XXX')
        if 'X' in seq:
            self.logger.info('There are unknown residue types in the structure')
        return seq

    def ReturnPdbFileName(self):
        if self.pdbfile[-4:] == '.pdb':
            return self.pdbfile[:-4]
        else:
            return 'None'
        
    def ReturnMolecularWeight(self, unit='Kd'):
        self.logger.info('Calculating molecular weight')
        mw = 0
        for aa in self.ReturnSeq(3):
            if aa in self.converter.keys():
                mw += self.converter[aa][2]
            else:
                mw += 110
        if unit == 'Kd':
            return round(mw / 1000, 1)
        elif unit == 'Da':
            return mw
        else:
            self.logger.error("ReturnMolecularWeight accepts units 'Kd' or 'Da'")
            return 0
        
    def ReturnExtinctionCoefficient(self, unit='absorbance'):
        self.logger.info('Calculating extinction coefficient')
        myseq = self.ReturnSeq(1)
        ec = myseq.count('W') * 5500.0 + myseq.count('Y') * 1490.0 + myseq.count('C') * 125.0
        if unit == 'absorbance':
            return ec / self.ReturnMolecularWeight('Da')
        elif unit == 'extinction':
            return ec
        else:
            self.logger.error("ReturnExtinctionCoefficient accepts units 'absorbance' or 'extinction'")
            return 0


    def Invert(self):
        self.logger.info('Inverting the structure')
        for index in self.hashdata.keys():
            if self.hashdata[index]['record_type'] == 'ATOM' or self.hashdata[index]['record_type'] == 'HETATM':
                self.hashdata[index]['z'] = self.hashdata[index]['z'] * -1
        
    def Rotate(self, option, origin=None):
        mytype = None
        translation = None
        if type(option) == type('string'):
            try:
                if not re.search('[xyzXYZ]', option[0]):
                    raise IOError('Wrong format input string for Rotate')
                axis = option[0].upper()
                theta = float(option[1:])
                mytype = 'string'
            except:
                self.logger.error("Rotation axis and angle must be in form 'x10' or 'Z2.5'")
                sys.exit()
        elif type(option) == type([]):
            try:
                for i1 in range(0,3):
                    for i2 in range(0,3):
                        float(option[i1][i2])
                option = numpy.array(option)
                mytype = 'matrix'
            except:
                self.logger.error("Rotation matrix must be a 3x3 nested array of floats")
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
        
        if mytype == 'string':
            self.logger.info('Rotating around '+str(axis)+' by '+str(theta)+' degrees')
        if mytype == 'matrix':
            self.logger.info('Rotating by a user input matrix')
        if origin:
            try:
                translation = [
                    -1 * float(origin[0]),
                    -1 * float(origin[1]),
                    -1 * float(origin[2])
                    ]
                self.Translate(translation)
                self.logger.info('Rotating around origin: '+','.join([str(x) for x in origin]))
            except:
                origin=False
                self.logger.error('Translate origin arg should be x,y,z coords as a tuple or list')
        else:
            self.logger.info('Rotating around origin: 0,0,0')
                                 
            
        for index in self.hashdata.keys():
            if self.hashdata[index]['record_type'] == 'ATOM' or self.hashdata[index]['record_type'] == 'HETATM':
                v = numpy.array([self.hashdata[index]['x'],self.hashdata[index]['y'],self.hashdata[index]['z']])
                if mytype == 'string':
                    newcoords = (numpy.dot(rotation_matrix(axis, theta),v))
                else:
                    newcoords = (numpy.dot(option,v))
                self.hashdata[index]['x'] = float(newcoords[0])
                self.hashdata[index]['y'] = float(newcoords[1])
                self.hashdata[index]['z'] = float(newcoords[2])
        if origin:
            self.Translate([-1*x for x in translation])
            
    def RenameChain(self, old, new):
        if not len(str(new)) == 1:
            self.logger.error('Chain names should be a single character')
        else:
            number_renamed = 0
            for index in self.hashdata.keys():
                if self.hashdata[index]['record_type'] == 'ATOM' or self.hashdata[index]['record_type'] == 'HETATM':
                    if self.hashdata[index]['chain'] == str(old):
                        self.hashdata[index]['chain'] = str(new)
                        number_renamed += 1
            self.logger.info('Renamed chain on '+str(number_renamed)+' residues')
                    
    def Translate(self, option):
        translation = None
        try:
            if type(option) in [type(()),type([])]:
                if len(option) == 3:
                    translation = [float(x) for x in option]
            elif type(option) == type(''):
                    translation = [float(x) for x in re.split('[,x ]', option)]
            else:
                raise IOError('wrong input type for Translate function')
        except:
            self.logger.error('Translation should be 3 floats in a tuple or array or a string with three floats delimited by , or x.')
            return False

        self.logger.info('Translating by '+str(translation[0])+'x'+str(translation[1])+'x'+str(translation[2]))
        for index in self.hashdata.keys():
            if self.hashdata[index]['record_type'] == 'ATOM' or self.hashdata[index]['record_type'] == 'HETATM':
                self.hashdata[index]['x'] = self.hashdata[index]['x'] + translation[0]
                self.hashdata[index]['y'] = self.hashdata[index]['y'] + translation[1]
                self.hashdata[index]['z'] = self.hashdata[index]['z'] + translation[2]
                
    def setCentre(self, centre=(0,0,0)):
        try:
            self.centre = ( float(centre[0]), float(centre[1]), float(centre[2]) )
        except:
            self.logger.error('Failed to set centre of mass, it should be X,Y,Z coords as a tuple')
            
    def returnCentre(self):
        return self.centre
    
    def centreOfMass(self):
        total_mass = 0
        total_x = 0
        total_y = 0
        total_z = 0
        for index in self.hashdata.keys():
            if self.hashdata[index]['record_type'] == 'ATOM' or self.hashdata[index]['record_type'] == 'HETATM':
                if self.hashdata[index]['element'] != '':
                    atom_name = self.hashdata[index]['element']
                else:
                    atom_name = self.hashdata[index]['atom_name'][0]
                try:
                    mass = periodictable.elements.symbol(atom_name).mass
                    total_mass += mass
                    total_x += (mass * self.hashdata[index]['x'])
                    total_y += (mass * self.hashdata[index]['y'])
                    total_z += (mass * self.hashdata[index]['z'])

                except:
                    self.logger.error('Failed to find mass for element '+str(atom_name))

        self.logger.info('Total mass is {:.2f} KDa'.format(total_mass / 1000))
        average_x = total_x / total_mass
        average_y = total_y / total_mass
        average_z = total_z / total_mass
        self.logger.info('Centre of mass is currently: {0:.3f},{1:.3f},{2:.3f}'.format(average_x,average_y,average_z))
        return (average_x,average_y,average_z)
    
    
    def parse_file(self):
        self.logger.info('Reading and parsing pdb file: '+str(self.pdbfile))

        ######################################
        #PDB DEFINITION FROM wwPDB GUIDLINES #
        # VERSION 3.30 31/07/14              #
        ######################################
        if self.pdbfile == None:
            pdblines = []
        else:
            file = open(self.pdbfile, 'r')
            pdblines = file.readlines()

        for line in pdblines:
            index = pdblines.index(line)
            self.hashdata[index] = {}
            record_type = line[self.pdb_definition['record_type'][0]:self.pdb_definition['record_type'][1]].rstrip()
            if record_type == 'ATOM' or record_type == 'HETATM':
                for record in self.pdb_definition.keys():
                    try:
                        if record in self.integer_records:
                            self.hashdata[index][record] = int(line[self.pdb_definition[record][0]:self.pdb_definition[record][1]].strip())
                        elif record in self.float_records.keys():
                            self.hashdata[index][record] = float(line[self.pdb_definition[record][0]:self.pdb_definition[record][1]].strip())
                        else:
                            self.hashdata[index][record] = line[self.pdb_definition[record][0]:self.pdb_definition[record][1]].strip()
                    except:
                        pass
            else:
                self.hashdata[index]['record_type'] = line[self.pdb_definition['record_type'][0]:self.pdb_definition['record_type'][1]].rstrip()
                self.hashdata[index]['string'] = line[self.pdb_definition['string'][0]:-1]+line[-1]


    def return_file(self, justatoms=False):
        self.logger.info('Writing out a formatted PDB file')
        if justatoms:
            justatoms = True
            self.logger.info('Will only return the atom lines')

        output_lines = []
        #ljust_records = [ 'record_type', 'atom_name', 'string' ]
        ljust_records = [ 'record_type', 'string' ]
        
        for index in sorted(self.hashdata.keys()):
            try:
                if len(self.hashdata[index]['atom_name']) < 4:
                    self.hashdata[index]['atom_name'] = ' '+self.hashdata[option][index]['atom_name']
            except:
                pass
            line = list(' '*80)
            for record in self.hashdata[index]:
                try:
                    length = self.pdb_definition[record][1] - self.pdb_definition[record][0]
                except:
                    length = 80 - self.pdb_definition[record][0]
                if record in self.float_records:
                    if record in ljust_records:
                        content = list((('%.'+str(self.float_records[record])+'f') % self.hashdata[index][record]).ljust(length))
                    else:
                        content = list((('%.'+str(self.float_records[record])+'f') % self.hashdata[index][record]).rjust(length))
                        
                else:
                    if record in ljust_records:
                        content = list(str(self.hashdata[index][record]).ljust(length))
                    else:
                        if record == 'atom_name':
                            temp_string = str(self.hashdata[index][record]).ljust(3)
                        else:
                            temp_string = str(self.hashdata[index][record])
                        content = list(temp_string.rjust(length))
            
                try:
                    line[self.pdb_definition[record][0]:self.pdb_definition[record][1]] = content[0:]
                except:
                    line[self.pdb_definition[record][0]:80] = content[0:]
            if justatoms:
                if self.hashdata[index]['record_type'] == 'ATOM' or self.hashdata[index]['record_type'] == 'HETATM':
                    output_lines.append(''.join(line))
            else:
                output_lines.append(''.join(line))


        return '\n'.join(output_lines)+'\n'




class DAT(Interface):
    """Read and write dat files
    
    The DAT class contains functions for reading dat files to
    a dictionary object where various calculations can be
    performed on it and then writing it back out again.
    """
    
    '''
    Constructor
    '''
    def __init__(self, datfile):
        ###start a log file
        self.logger = logging.getLogger('readwrite.DAT')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(streamhandler)
            self.logger.info('Starting a new readwrite.DAT job')        
        self.type = 'dat'
        if os.path.isfile(datfile) and datfile[-4:] == '.dat':
            self.datfile = datfile
        else:
            sys.exit(str(datfile)+' either does not exist or is not of type ".dat"')

        self.hashdata = {}
        self.hi_q_window_range = (0.68,0.79)
        self.lo_q_window_range = (0.02,0.05)
        self.lo_q_window = 0
        self.hi_q_window = 0
        
    def parse_file(self):
        self.logger.info('Reading and parsing dat file: '+str(self.datfile))
        qdata = []
        idata = []
        edata = []
                
        datafile = open(self.datfile, "r")
        saxsdata = datafile.readlines()
        for row in saxsdata:
            try:
                q = float(row.split()[0])
                i = float(row.split()[1])
                e = float(row.split()[2])
                if q > 0:
                    qdata.append(q)
                    idata.append(i)
                    edata.append(e)
            except:pass
        self.hashdata = {'Q': qdata, 'I': idata, 'E': edata}
        x1 = int(round(len(self.hashdata['Q'])*self.lo_q_window_range[0]))
        x2 = int(round(len(self.hashdata['Q'])*self.lo_q_window_range[1]))
        self.lo_q_window = sum(self.hashdata['I'][x1:x2])
        x1 = int(round(len(self.hashdata['Q'])*self.hi_q_window_range[0]))
        x2 = int(round(len(self.hashdata['Q'])*self.hi_q_window_range[1]))
        self.hi_q_window = sum(self.hashdata['I'][x1:x2])
        self.logger.info('Parsed with '+str(len(qdata))+' data points.')

    def ReturnDataColumn(self, column='I'):
        if column in ['Q', 'I', 'E']:
            return self.hashdata[column]
        else:
            self.logger.error('ReturnDataColumn function requires you specify either Q, I or E')
            
    def returnIEData(self, q=None):
        try:
            index = self.hashdata['Q'].index(q)
            return (self.hashdata['I'][index],self.hashdata['E'][index])
        except:
            return False

    def ReturnDataColumn(self, column='I'):
        if column in ['Q', 'I', 'E']:
            return self.hashdata[column]
        else:
            self.logger.error('ReturnDataColumn function requires you specify either Q, I or E')            
    def return_file(self):
        self.logger.info('Writing out a formatted DAT file')
        string_list = []
        string_list.append("%-15s %-18s %-15s" % ("Q(A-1)","I(au)","Error"))
        for q in self.hashdata['Q']:
            if q > 0:
                index = self.hashdata['Q'].index(q)
                q = self.hashdata['Q'][index]
                i = self.hashdata['I'][index]
                e = self.hashdata['E'][index]
                string_list.append("%-15s %-18s %-15s" % (q,i,e))
        return '\n'.join(string_list)


    def input_dict(self, input_dict):
        self.logger.info('Reading in DAT data as a dictionary')
        if type(input_dict) == type({}):
            self.hashdata = input_dict

    def return_dict(self):
        self.logger.info('Returning DAT data as a dictionary')
        return self.hashdata



class FIT(Interface):
    """Read and write dat files
    
    The FIT class contains functions for reading fit files
    (i.e. from Crysol) to a dictionary object where various
    calculations can be performed on it and then writing it
    back out again.
    """
    
    '''
    Constructor
    '''
    def __init__(self, fitfile):
        ###start a log file
        self.logger = logging.getLogger('readwrite.FIT')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(streamhandler)
            self.logger.info('Starting a new readwrite.FIT job')        

        if os.path.isfile(fitfile) and fitfile[-4:] == '.fit':
            self.fitfile = fitfile
        else:
            sys.exit(str(fitfile)+' either does not exist or is not of type ".fit"')

        self.hashdata = {}

    def parse_file(self):
        self.logger.info('Reading and parsing fit file: '+str(self.fitfile))
        qdata = []
        exp_data = []
        model_data = []
                
        datafile = open(self.fitfile, "r")
        saxsdata = datafile.readlines()
        for row in saxsdata:
            try:
                q = float(row.split()[0])
                exp = float(row.split()[1])
                mod = float(row.split()[2])
                if q > 0 and exp > 0:
                    qdata.append(q)
                    exp_data.append(exp)
                    model_data.append(mod)
            except:pass
        self.hashdata = {'Q': qdata, 'obs': exp_data, 'mod': model_data}
        self.logger.info('Parsed '+str(self.fitfile)+' with '+str(len(qdata))+' data points.')

    def return_file(self):
        self.logger.info('Writing out a formatted FIT file')
        string_list = []
        string_list.append("%-15s %-18s %-15s" % ("Q(A-1)","Observed","Model"))
        for q in self.hashdata['Q']:
            if q > 0:
                index = self.hashdata['Q'].index(q)
                q = self.hashdata['Q'][index]
                obs = self.hashdata['obs'][index]
                mod = self.hashdata['mod'][index]
                string_list.append("%-15s %-18s %-15s" % (q,obs,mod))
        return '\n'.join(string_list)

    def input_dict(self, input_dict):
        self.logger.info('Reading in DAT data as a dictionary')
        if type(input_dict) == type({}):
            self.hashdata = input_dict

    def return_dict(self):
        self.logger.info('Returning DAT data as a dictionary')
        return self.hashdata



class OUT(Interface):
    """Read and write out files
    
    The OUT class contains functions for reading out files
    (i.e. from GNOM) to a dictionary object where various
    calculations can be performed on it and then writing it
    back out again.
    """
    
    '''
    Constructor
    '''
    def __init__(self, outfile):
        ###start a log file
        self.logger = logging.getLogger('readwrite.OUT')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(streamhandler)
            self.logger.info('Starting a new readwrite.OUT job')        

        if os.path.isfile(outfile) and outfile[-4:] == '.out':
            self.outfile = outfile
        else:
            sys.exit(str(outfile)+' either does not exist or is not of type ".out"')

        self.hashdata = {}

    def Dmax(self):
        try:
            dmax = self.hashdata['R'][-1]
        except:
            dmax = 0
            self.logger.error('Failed to get Dmax from P(r)')
        return dmax

    def numberOfBins(self):
        try:
            bins = len(self.hashdata['R'])
        except:
            bins = 0
            self.logger.error('Failed to get number of bins from P(r)')
        return bins
        
    def parse_file(self):
        self.logger.info('Reading and parsing out file: '+str(self.outfile))
        rdata = []
        prdata = []
        edata = []

        datafile = open(self.outfile, "r")
        saxsdata = datafile.readlines()
        try:
            index = [i for i, item in enumerate(saxsdata) if re.search(' *R          P\(R\)      ERROR.*', item)][-1]
        except:
            self.logger.error('Outfile '+str(self.outfile)+'.out does not seem to have the expected format for an outfile')
        for row in saxsdata[index:]:
            try:
                r = float(row.split()[0])
                pr = float(row.split()[1])
                e = float(row.split()[2])
                rdata.append(r)
                prdata.append(pr)
                edata.append(e)
            except:pass
        self.hashdata = {'R': rdata, 'PR': prdata, 'E': edata}
        self.logger.info('Parsed '+str(self.outfile)+' with '+str(len(rdata))+' data points.')



    def return_file(self):
        self.logger.info('Writing out a formatted OUT file')
        string_list = []
        string_list.append("%-15s %-18s %-15s" % ("R(A)","P(R)","ERROR"))
        for r in self.hashdata['R']:
            index = self.hashdata['R'].index(r)
            r = self.hashdata['R'][index]
            pr = self.hashdata['PR'][index]
            e = self.hashdata['E'][index]
            string_list.append("%-15s %-18s %-15s" % (r,pr,e))
        return '\n'.join(string_list)

    def input_dict(self, input_dict):
        self.logger.info('Reading in OUT data as a dictionary')
        if type(input_dict) == type({}):
            self.hashdata = input_dict

    def return_dict(self):
        self.logger.info('Returning OUT data as a dictionary')
        return self.hashdata
