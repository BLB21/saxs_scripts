#!/anaconda/bin/python
'''
Created on Jul 17, 2014

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
required.add_option("-s", "--sample", action="store", type="string", dest="sample", help="The dat file of the sample")
required.add_option("-b", "--buffer", action="store", type="string", dest="buffer", help="The dat file of the buffer")
optional = OptionGroup(parser, "Optional Arguments")
optional.add_option("-m", "--multiplier", action="store", type="float", dest="multiplier", default=1.0, help="the multiplier will be applied to the buffer dat before subtraction")
optional.add_option("-a", "--auto", action="store_true", dest="auto", default=False, help="if this flag is present the script will attempt to automatically choose a multiplier")
required.add_option("-o", "--output", action="store", type="string", dest="output", help="The name of an output dat file. Default is 'samplename_sub.dat'")
parser.add_option_group(required)
parser.add_option_group(optional)
(options, args) = parser.parse_args()

'''
fail if you didn't choose valid dat files
'''

if os.path.isfile(options.sample) and os.path.isfile(options.buffer):
    rootname = os.path.splitext(os.path.basename(options.sample))[0]
    dir = os.path.split(os.path.realpath(options.sample))[0]
    if not options.output:
        options.output = dir+'/'+rootname+'_sub.dat'
else:
    sys.exit('One or both of the dat files you specified do not exist')

class Subtract():
    """Subtract a buffer dat file from a sample dat file with optional correction of blanking errors
    
    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('Subtract')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new Subtraction job')        
        try:
            self._options = dict(options)
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()

    def ParseDatFiles(self):
        self.logger.info('Reading the dat files')
        self.imagedict = {}
        for file in ['buffer', 'sample']:
            qdata = []
            idata = []
            edata = []
                    
            datafile = open(self._options[file], "r")
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
            self.imagedict[file] = {'Q': qdata, 'I': idata, 'E': edata}
            self.logger.info('Parsed '+str(file)+' with '+str(len(qdata))+' data points.')
        if len(self.imagedict['buffer']['Q']) == len(self.imagedict['sample']['Q']):
            self.logger.info('buffer and sample have the same number of points, this is a bonus!')
        else:
            self.logger.error('the sample and buffer do not have the same number of points. This is disappointing!')

    def DetermineMultiplier(self):
        if self._options['auto']:
            self.logger.info('Determining a multiplier to correct blanking errors')
            fraction = 0.01
            no_datapoints = int(round((len(self.imagedict['sample']['I']) * fraction)))
            self._options['multiplier'] = numpy.average(self.imagedict['sample']['I'][-no_datapoints:]) / numpy.average(self.imagedict['buffer']['I'][-no_datapoints:])
            self.logger.info('will multiply buffer by '+str(round(self._options['multiplier'],4))+' before subtracting from sample')
        elif self._options['multiplier'] != 1.0:
            self.logger.info('will multiply buffer by '+str(self._options['multiplier'])+' before subtracting from sample')
        else:
            self.logger.info('will not correct sample before buffer subtraction')

    def Subtraction(self):
        self.string_list = []
        self.string_list.append("%-15s %-18s %-15s" % ("Q(A-1)","I(au)","Error"))
        for q in self.imagedict['sample']['Q']:
            index = self.imagedict['sample']['Q'].index(q)
            q = self.imagedict['sample']['Q'][index]
            i = self.imagedict['sample']['I'][index] - ( self._options['multiplier'] * self.imagedict['buffer']['I'][index] )
            e = numpy.sqrt(numpy.power(self.imagedict['sample']['E'][index],2) + numpy.power(self._options['multiplier'] * self.imagedict['buffer']['E'][index],2))
            self.string_list.append("%-15s %-18s %-15s" % (q,i,e))

    def OutputFile(self):
        self.logger.info('wrote output to '+str(self._options['output']))
        outfile = open(self._options['output'], 'w')
        outfile.write('\n'.join(self.string_list))
        outfile.close()

        
                                                       


if __name__ == '__main__':
    options = eval(str(options))
    job = Subtract(options)
    job.ParseDatFiles()
    job.DetermineMultiplier()
    job.Subtraction()
    job.OutputFile()
