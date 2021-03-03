#!/anaconda/bin/python
'''
Created on Apr 03, 2014

@author: nathan
'''

import sys, logging, os, re, math, numpy, random
from subprocess import check_output

class Simulate():
    """Simulate SAXS data from a PDB file
    
    The Simulate class contains various functions for calculating SAXS
    data from a pdb file. The script attempts to generate the kind of
    errors and noise that you would find in real, measured SAXS data
    as an aid to experimental design.
    """
    
    '''
    Constructor
    '''
    def __init__(self, options):
        ###start a log file
        self.logger = logging.getLogger('Simulate')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new Simulate job')        
        try:
            self._options = dict(options)
        except:
            self.logger.error('cound not read in the command line options')
            sys.exit()
#        try:
#            self.qmin = float(min(self._options['qrange'].split('-')))
#            self.qmax = float(max(self._options['qrange'].split('-')))
#        except:
#            self.logger.error('Q min and Q max were not defined correctly. Must# be in the format i.e. 0.01-0.55')
#            sys.exit()

    def RunCrysol(self):
        self.logger.info('Running Crysol')
        filelist_before = os.listdir(os.getcwd())
        command = 'crysol -ns '+str(self._options['no_points'])+' -sm '+str(self.highq)+' '+str(self._options['file'])
        self.output = check_output(command, shell=True)
        filelist_after = os.listdir(os.getcwd())
        self.crysol_files = list(set(filelist_after) - set(filelist_before))

    def ParseIntFile(self):
        for file in self.crysol_files:
            if file[-4:] == '.int':
                intfile=open(file)
                intdata=intfile.readlines()
                self.qdata = []
                self.idata = []
                if self._options['vacuum']:
                    self.logger.info('Outputting the in vacuum data')
                for row in intdata:
                    try:
                        q = float(row.split()[0])
                        if self._options['vacuum']:
                            i = float(row.split()[2]) + ( self._options['background'] * 1E6 )
                        else:
                            i = float(row.split()[1]) + ( self._options['background'] * 1E6 )
                        if q > self.lowq and q < self.highq:
                            self.qdata.append(q)
                            self.idata.append(i)
                        if q == 0:
                            self.izero = i
                    except:
                        pass
                    
                self.logger.info('From crysol,     min Q: '+str(min(self.qdata))+', max Q: '+str(max(self.qdata)))

                try:
                    float(self.izero)
                except:
                    self.logger.error('could not determine I(0) from int file')
                    sys.exit()

        ###GET MW
        try:
            self.output = self.output.split('\n')
            pattern = re.compile('.*Molecular weight:.*')
            for line in self.output:
                if re.match(pattern, line):
                    self.molwt = float(line.split()[2]) / 1000
        except:
            self.logger.error('could not get MW from crysol output')
            sys.exit()


    def SFtoPhotons(self):
        pixel_size = 0.000172 # in meters
        camera_length = self._options['camera'] / 1000 # in meters
        solid_angle = (pixel_size / camera_length)**2 #steradians
        attenuation = self._options['attenuation']
        full_flux = 2.00E+13
        intensity = full_flux * attenuation
        exposure_time = self._options['time']
        electron_radius = 2.82E-15
        avogadros = 6.02E+23
        path_length = self._options['pathlength']/1000 # in meters
        concentration = self._options['concentration']
        mol_weight = self.molwt *1000
        water_sf = 6.6
        water_conc = 1000
        water_mw = 18
        self.proteinphotons = []
        self.i_photons = []
        for i in self.idata:
            protein = solid_angle * intensity * exposure_time * (electron_radius**2)*avogadros*1000*concentration*path_length/mol_weight*i
            water = solid_angle * intensity * exposure_time * (electron_radius**2)*avogadros*1000*water_conc*path_length/water_mw*water_sf

            self.i_photons.append(protein + water)
            self.proteinphotons.append(protein)

    def ScaleIntensities(self):
        ###CALCULATE EXPECTED I(0) BASED ON MW AND CONC
        avogadros = 6.022E+23
        io_calc = ( self.molwt * (self._options['sld'] * self._options['vbar'])**2 * self._options['concentration'] ) / avogadros
        self.logger.info('calculated molecular weight at '+str(self.molwt))
        try:
            scale_factor = self.izero / io_calc
            scaled_idata = []
            for i in self.idata:
                i = (i / scale_factor)
                scaled_idata.append(i)
            self.idata = scaled_idata
        except:
            scaled_idata = []
            for i in self.idata:
                scaled_idata.append(0)
            self.idata = scaled_idata
                

    def SimulateImage(self):
        self.logger.info('Simulating pilatus image')
        beamx = 559.94
        beamy = 467.5
        sizex = 981
        sizey = 1043
        pixelsize = 0.172
        camera_length = self._options['camera']
        energy = self._options['energy']
        wavelength = 12.398521 / energy
        beamstop_radius = 15 #(pixels)
        beamstop_arm_gradient = -0.4673
        beamstop_arm_width = 9
        horizontal_module_boundaries = [ ( 195, 212 ), ( 407, 424 ), ( 619, 636 ), ( 831, 848 ) ]
        vertical_module_boundaries = [ ( 487, 494 ) ]

        ###A LIST OF ALL PIXEL COORDINATES
        pixel_array = []
        for xcoord in range(1, sizex+1):
            for ycoord in range(1, sizey+1):
                pixel_array.append( ( xcoord, ycoord ) )

        ###REMOVE PIXELS THAT ARE MASKED BY THE HORIZONTAL MODULE BOUNDARIES
        delete_list = []
        for mask in horizontal_module_boundaries:
            for pixel in pixel_array:
                if pixel[1] > mask[0] and pixel[1] < mask[1]:
                    delete_list.append(pixel)
        s = set(delete_list)
        new_array = [x for x in pixel_array if x not in s]
        pixel_array = new_array

        ###REMOVE PIXELS THAT ARE MASKED BY THE VERTICAL MODULE BOUNDARIES
        delete_list = []
        for mask in vertical_module_boundaries:
            for pixel in pixel_array:
                if pixel[0] > mask[0] and pixel[0] < mask[1]:
                    delete_list.append(pixel)
        s = set(delete_list)
        new_array = [x for x in pixel_array if x not in s]
        pixel_array = new_array

        ###DELETE ITEMS BEHIND BACKSTOP ARM
        delete_list = []
        for pixel in pixel_array:
            if pixel[0] > beamx:
                lowerlimit = beamstop_arm_gradient * pixel[0] + ((beamy+262.3) - (beamstop_arm_width/2))
                upperlimit = beamstop_arm_gradient * pixel[0] + ((beamy+262.3) + (beamstop_arm_width/2))
                if pixel[1] > lowerlimit and pixel[1] < upperlimit:
                    delete_list.append(pixel)
        s = set(delete_list)
        new_array = [x for x in pixel_array if x not in s]
        pixel_array = new_array

        ###DELETE ITEMS BEHIND THE BACKTOP SHADOW
        delete_list = []
        for pixel in pixel_array:
            distance = math.sqrt(( pixel[0] - beamx )**2 + ( pixel[1] - beamy )**2 )
            if distance < beamstop_radius:
                delete_list.append(pixel)
        s = set(delete_list)
        new_array = [x for x in pixel_array if x not in s]
        pixel_array = new_array
        self.q_array = []
        for pixel in pixel_array:
            ###CALCULATE PIXEL DISTANCE FROM DIRECT BEAM
            distance_pixels = math.sqrt(( pixel[0] - beamx )**2 + ( pixel[1] - beamy )**2 )
            distance_mm = distance_pixels * pixelsize
            angle_theta = math.asin( distance_mm / camera_length )/2#in rad
            angle_q = 4*math.pi*math.sin(angle_theta)/wavelength
            self.q_array.append(angle_q)

        self.lowq = min(self.q_array)
        self.highq = max(self.q_array)
        self.logger.info('From simulation, min Q: '+str(round(self.lowq, 6))+', max Q: '+str(round(self.highq, 6)))

    def GenerateErrors(self):
        self.logger.info('Generating errors and adding noise to data')
        bins = numpy.linspace(self.lowq, self.highq, len(self.qdata))
        binned = numpy.histogram(self.q_array, bins=len(self.qdata))
        self.binsizes = binned[0]

        self.noisy_i = []
        self.edata = []
        for q in self.qdata:
            index = self.qdata.index(q)
            i_inv_cm = self.idata[index]
            abs_cal_scale = 0.00072*2  #0.0142
            i_photons = (i_inv_cm / abs_cal_scale)*self._options['time']
            
            error = (2*(math.sqrt(i_photons)/math.sqrt(self.binsizes[index])) * abs_cal_scale) / self._options['time']



            new_i = random.gauss(self.idata[index], error)
            self.noisy_i.append(new_i)
            self.edata.append(error)

    def GenerateErrors2(self):
        self.logger.info('Generating errors and adding noise to data')
        bins = numpy.linspace(self.lowq, self.highq, len(self.qdata))
        binned = numpy.histogram(self.q_array, bins=len(self.qdata))
        self.binsizes = binned[0]

        self.noisy_i = []
        self.edata = []
        for q in self.qdata:
            index = self.qdata.index(q)
            if q == 0.300068:
                print str(q)+','+str(self.binsizes[index])

            i_inv_cm = self.idata[index]
            background_error = 10E-7
            photons_per_bin = self.i_photons[index] * self.binsizes[index]
            proteinphotons_per_bin = self.proteinphotons[index] * self.binsizes[index]
            total_photons = 2* photons_per_bin - proteinphotons_per_bin
            error_in_photons = math.sqrt(total_photons + 0.0030**2*total_photons**2)
            try:
                photon_to_invcm_scale = (proteinphotons_per_bin / self.idata[index])
                error = (error_in_photons / photon_to_invcm_scale) + background_error
            except:
                error = background_error
#            std_error = error / math.sqrt(self.binsize[index]-1)
            #error = (2*(math.sqrt(self.i_photons[index])/math.sqrt(self.binsizes[index])))
            #scale = self.i_photons[index]/self.idata[index]
            #error = 64* error/scale + background_error #double it because of buffer subtraction
            new_i = random.gauss(self.idata[index], error/2)
            self.noisy_i.append(new_i)
            self.edata.append(error)
        
    def OutputDatFile(self):
        string_list = []
        string_list.append("Simulated data from "+self._options['file'])
        string_list.append("%-15s %-18s %-15s" % ("Q(A-1)","I(au)","Error"))
        for q in self.qdata:
            index = self.qdata.index(q)
            q = self.qdata[index]
            i = self.noisy_i[index]
            e = self.edata[index]
            string_list.append("%-15s %-18s %-15s" % (q,i,e))
        outfile = open(self._options['outfile'], 'w')
        outfile.write('\n'.join(string_list))
        outfile.close()
        self.logger.info('Output file '+str(self._options['outfile']))

                           

    def CleanUpFiles(self):
        for file in self.crysol_files:
            os.remove(file)
        self.logger.info('Deleted all files made by Crysol for this run')

if __name__ == '__main__':
    from optparse import OptionParser
    from optparse import OptionGroup
    
    if len(sys.argv) < 2:
        sys.argv.append('-h')
        
    '''
    parse command line options
    '''
    
    parser = OptionParser()
    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-f", "--file", action="store", type="string", dest="file", help="The pdb file you want to use to simulate data")
    optional = OptionGroup(parser, "Optional Arguments")
    optional.add_option("-c", "--concentration", action="store", type="float", dest="concentration", default=1.0, help="The concentration in mg/ml of the protein (default 1 mg/ml)")
    optional.add_option("-t", "--time", action="store", type="float", dest="time", default=10.0, help="The exposure time for the sample, (default 10 secs)")
    
    optional.add_option("-s", "--sld", action="store", type="float", dest="sld", default=2.67E+10, help="The scattering length density of the protein (default is 2.67E+10)")
    optional.add_option("-v", "--vbar", action="store", type="float", dest="vbar", default=0.73, help="The partial specific volume of the protein (default is 0.73)")
    
    optional.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="The name of an output file, the default name is the rootname of your pdb file with .dat at the end.")
    
    optional.add_option("-n", "--number", action="store", type="int", dest="no_points", default=455, help="The number of points in the output file, (default 455)")
    
    optional.add_option("-l", "--length", action="store", type="float", dest="camera", default=1492, help="The length of the SAXS camera in mm, (default 1492)")
    
    optional.add_option("-p", "--pathlength", action="store", type="float", dest="pathlength", default=1.5, help="The pathlength, ie width of capillary in mm, (default 1.5)")
    
    optional.add_option("-e", "--energy", action="store", type="float", dest="energy", default=11, help="The wavelength in KeV, (default 11)")
    optional.add_option("-a", "--attenuation", action="store", type="float", dest="attenuation", default=0.3, help="The fractional percentage that the beam was attenuated (default 0.3)")
    
    optional.add_option("-q", "--quiet", action="store_false", dest="quiet", default=True, help="don't print to st. out or display graphics")
    optional.add_option("-d", "--dehydrated", action="store_true", dest="vacuum", default=False, help="Use the in vacuum column instead of the hydrated one. Default is hydrated")    
    optional.add_option("-b", "--background", action="store", type="float", dest="background", default=0, help="simulate a % error in background, (default is 0). i.e. 1.01 would oversubtract a background and 0.99 would undersubtract by 1%")
    parser.add_option_group(required)
    parser.add_option_group(optional)
    (options, args) = parser.parse_args()
    
    '''
    fail if you didn't choose a valid pdb file
    '''
    
    if os.path.isfile(options.file):
        rootname = os.path.splitext(os.path.basename(options.file))[0]
        cwd = os.path.split(os.path.realpath(options.file))[0]
        os.chdir(cwd)
        if not options.outfile:
            options.outfile = rootname+'.dat'
        pass
    else:
        sys.exit('The pdb file you specified does not exist')

    options = eval(str(options))
    job = Simulate(options)
    job.SimulateImage()
    job.RunCrysol()
    job.ParseIntFile()
    job.SFtoPhotons()
    job.ScaleIntensities()
    job.GenerateErrors2()
    job.OutputDatFile()
    job.CleanUpFiles()
    
