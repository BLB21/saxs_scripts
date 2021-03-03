#!/anaconda/bin/python
'''
Created on Apr 03, 2014

@author: nathan
'''

import sys, logging, os, re, math, numpy, random
from optparse import OptionParser
from optparse import OptionGroup
from subprocess import check_output


if len(sys.argv) < 2:
    sys.argv.append('-h')
    
'''
parse command line options
'''

parser = OptionParser()
required = OptionGroup(parser, "Required Arguments")
required.add_option("-f", "--firstfile", action="store", type="string", dest="firstfile", help="The first dat file you want to sum together")
required.add_option("-s", "--secondfile", action="store", type="string", dest="secondfile", help="The second dat file you want to sum together")
optional = OptionGroup(parser, "Optional Arguments")
optional.add_option("-1", "--firstmultiplier", action="store", type="float", dest="firstmultiplier", default=0.5, help="A multiplier to apply to the first dat file before combining, default is 0.5")
optional.add_option("-2", "--secondmultiplier", action="store", type="float", dest="secondmultiplier", default=0, help="A multiplier to apply to the second dat file before combining, default is ( 1.0 - first multiplier )")
optional.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="The name of an output file, the default name is the indeces of your two files with the multipliers.")
parser.add_option_group(required)
parser.add_option_group(optional)
(options, args) = parser.parse_args()

'''
fail if you didn't choose a valid pdb file
'''

if options.secondmultiplier == 0:
    options.secondmultiplier = 1.0 - options.firstmultiplier

if os.path.isfile(options.firstfile) and os.path.isfile(options.secondfile):
    output_dir = os.path.split(os.path.abspath(options.firstfile))[0]
    options.outputdir = output_dir
    first_index = re.split('[\._]', os.path.split(os.path.abspath(options.firstfile))[1])[-2]
    second_index = re.split('[\._]', os.path.split(os.path.abspath(options.secondfile))[1])[-2]
    default_output = str(options.firstmultiplier).replace('.','p')+'x'+first_index+'+'+str(options.secondmultiplier).replace('.','p')+'x'+second_index+'.dat'
    if not options.outfile:
        options.outfile = default_output
    else:
        output_dir = os.path.split(os.path.abspath(options.outfile))[0]
        options.outputdir = output_dir
        options.outfile = os.path.split(os.path.abspath(options.outfile))[1]
else:
    sys.exit('One or both of the dat files you specified do not exist')




class Combine():
    """Combine two dat files with option to apply multipliers to each file

    The main purpose of this class is for the case where two buffers have
    been used to make a serial dilution and you want to interpolate between
    the two to match intervening sample files.

    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('Combine')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new Combine job')        
        try:
            self._options = dict(options)
            print self._options
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()


    def ParseFirst(self):
        qdata = []
        idata = []
        edata = []
                
        datafile = open(self._options.get('firstfile'), "r")
        saxsdata = datafile.readlines()
        for row in saxsdata:
            try:
                q = float(row.split()[0])
                i = float(row.split()[1]) * self._options.get('firstmultiplier')
                e = float(row.split()[2]) * self._options.get('firstmultiplier')
                if q > 0:
                    qdata.append(q)
                    idata.append(i)
                    edata.append(e)
            except:pass
        self.firstdict = {'Q': qdata, 'I': idata, 'E': edata}
        self.logger.info('Parsed '+str(self._options.get('firstfile'))+' with '+str(len(qdata))+' data points.')

    def ParseSecond(self):
        qdata = []
        idata = []
        edata = []
                
        datafile = open(self._options.get('secondfile'), "r")
        saxsdata = datafile.readlines()
        for row in saxsdata:
            try:
                q = float(row.split()[0])
                i = float(row.split()[1]) * self._options.get('secondmultiplier')
                e = float(row.split()[2]) * self._options.get('secondmultiplier')
                if q > 0:
                    qdata.append(q)
                    idata.append(i)
                    edata.append(e)
            except:pass
        self.seconddict = {'Q': qdata, 'I': idata, 'E': edata}
        self.logger.info('Parsed '+str(self._options.get('secondfile'))+' with '+str(len(qdata))+' data points.')

    def Combine(self):
        qdata = []
        idata = []
        edata = []

        for q in self.firstdict['Q']:
            index = self.firstdict['Q'].index(q)
            if q in self.seconddict['Q']:
                index2 = self.seconddict['Q'].index(q)
                i = self.firstdict['I'][index] + self.seconddict['I'][index2]
                e = math.sqrt(math.pow(self.firstdict['E'][index], 2) + math.pow(self.seconddict['E'][index2], 2))
                qdata.append(q)
                idata.append(i)
                edata.append(e)
        self.outputdict = {'Q': qdata, 'I': idata, 'E': edata}
        self.logger.info('Summed the two files with '+str(len(self.outputdict['Q']))+' common points')                

    def OutputDat(self):
        string_list = []
        string_list.append("%-15s %-18s %-15s" % ("Q(A-1)","I(au)","Error"))
        for q in self.outputdict['Q']:
            index = self.outputdict['Q'].index(q)
            i = self.outputdict['I'][index]
            e = self.outputdict['E'][index]
            string_list.append("%-15s %-18s %-15s" % (q,i,e))
 
        outfile_name = self._options.get('outputdir')+'/'+self._options.get('outfile')
        outfile = open(outfile_name, 'w')
        outfile.write('\n'.join(string_list))
        outfile.close()
        self.logger.info('wrote averaged file to: '+outfile_name)

if __name__ == '__main__':
    options = eval(str(options))
    job = Combine(options)
    job.ParseFirst()
    job.ParseSecond()
    job.Combine()
    job.OutputDat()
