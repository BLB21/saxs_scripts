#!/anaconda/bin/python
import logging, re, os, sys, threading, time, warnings, subprocess, numpy#, matplotlib
import matplotlib.pyplot as plt
from scipy import stats as stats
from matplotlib.backends.backend_pdf import PdfPages


class SaxsPlot():
    """Plots saxs data in various ways
    
    The SaxsPlot libraries allow publication quality plotting of
    various SAXS data.

    
    """
    
    '''
    Constructor
    '''
    def __init__(self):
        self.output_images = {}
        ###start a log file
        self.logger = logging.getLogger('SaxsPlot')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new SaxsPlot job')


    def GetImages(self):
        ############################################
        #USER SUPPLIED ARGS: USE THESE AS FILE LIST#
        ############################################
        self.imagelist = []
        if len(sys.argv) > 1:
            for arg in sys.argv[1:]:
                if os.path.isfile(arg):
                    self.imagelist.append(os.path.abspath(arg))

        self.logger.info('User has chosen '+str(len(self.imagelist))+' files.')

        if len(self.imagelist) > 0:
            self.output_directory=os.path.dirname(self.imagelist[0])
            self.jobtype = self.imagelist[0].split('.')[-1]
            for path in self.imagelist:
                if path.split('.')[-1] != self.jobtype:
                    self.logger.error('If you enter multiple files they must all be of the same type')
                    sys.exit()
            if self.jobtype == 'fir':
                self.jobtype = 'fit'
            if self.jobtype == 'dat' and len(self.imagelist) > 20:
                self.jobtype = 'col'
            if self.jobtype == 'dat' or self.jobtype == 'out' or self.jobtype == 'fit':
                self.logger.info('will run code for displaying files of type '+str(self.jobtype))
            elif self.jobtype == 'col':
                self.logger.info('will run code for processing a size exclusion series')
            else:
                self.logger.error("can only parse files of type 'dat', 'out' or 'fit/fir'")
                sys.exit()
        else:
            self.logger.error('User did not choose any images, nothing to do!')
            sys.exit()

    def ParseDat(self):
        self.imagedict = {}
        for file in self.imagelist:
            fileroot = os.path.splitext(os.path.basename(file))[0]
            qdata = []
            idata = []
            edata = []
                    
            datafile = open(file, "r")
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
            self.imagedict[fileroot] = {'Q': qdata, 'I': idata, 'E': edata}
            self.logger.info('Parsed '+str(fileroot)+' with '+str(len(qdata))+' data points.')

    def ParseFit(self):
        self.imagedict = {}
        for file in self.imagelist:
            fileroot = os.path.splitext(os.path.basename(file))[0]
            qdata = []
            exp_data = []
            model_data = []
                    
            datafile = open(file, "r")
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
            self.imagedict[fileroot] = {'Q': qdata, 'I': exp_data, 'E': model_data}
            self.logger.info('Parsed '+str(fileroot)+' with '+str(len(qdata))+' data points.')



    def ParseOut(self):
        self.imagedict = {}
        for file in self.imagelist:
            fileroot = os.path.splitext(os.path.basename(file))[0]
            qdata = []
            idata = []
            edata = []

            datafile = open(file, "r")
            saxsdata = datafile.readlines()
            try:
                index = [i for i, item in enumerate(saxsdata) if re.search(' *R          P\(R\)      ERROR.*', item)][-1]+2
                index2 = [i for i, item in enumerate(saxsdata[index:]) if re.search(' *Reciprocal space: Rg =.*', item)][0] + index - 1
                saxsdata = saxsdata[index:index2]
                for row in saxsdata:
                    try:
                        q = float(row.split()[0])
                        i = float(row.split()[1])
                        e = float(row.split()[2])
                        qdata.append(q)
                        idata.append(i)
                        edata.append(e)
                    except:pass
                self.imagedict[fileroot] = {'Q': qdata, 'I': idata, 'E': edata}
                self.logger.info('Parsed '+str(fileroot)+' with '+str(len(qdata))+' data points.')

            except:
                self.logger.error('Outfile '+str(fileroot)+'.out does not seem to have the expected format for an outfile')

    def ScaleIntensities(self):
	self.logger.info('Scaling intensity data')
        time.sleep(5)
        if len(self.imagedict.keys()) > 0:
            highest_max = 0
            for file in self.imagedict.keys():
                self.imagedict[file]['max'] = max(self.imagedict[file]['I'])
	        if self.imagedict[file]['max'] > highest_max:
	            highest_max = self.imagedict[file]['max']
	    for file in self.imagedict.keys():
                self.logger.info('scaling '+str(file)+' intensities by '+str(highest_max/self.imagedict[file]['max']))
		self.imagedict[file]['I'][:] = [x*(highest_max/self.imagedict[file]['max']) for x in self.imagedict[file]['I']]

    def Plot(self):
        linetypes = ['-','--','-.',':']
        linecolours = ['k','b','g','r','m','c']
        if len(self.imagedict.keys()) < 5:
            linetype_array = ['k-', 'k--', 'k-.', 'k:']
        else:
            linetype_array = []
            for type in linetypes:
                for colour in linecolours:
                    linetype_array.append(str(colour)+str(type))
        fig = plt.figure()

        if self.jobtype == 'col':
#            plt.interactive(True)
#            fig = plt.figure()
            ax1 = fig.add_subplot(111, autoscale_on=True)
            ax1.plot(self.filenumbers, self.ios, 'k-')
            ax1.set_xlabel('file number')
            # Make the y-axis label and tick labels match the line color.
            ax1.set_ylabel('I(0)', color='k')
            for tl in ax1.get_yticklabels():
                tl.set_color('k')

            ax2 = ax1.twinx()
            ax2.plot(self.filenumbers, self.rgs, color='1', ls='None', marker='o', markersize=5)
            ax2.set_ylabel('Rg (A)', color='k')
            for tl in ax2.get_yticklabels():
                tl.set_color('k')
#            plt.title('Click to close!')
            plt.show()
#            plt.waitforbuttonpress()
#            plt.clf()
            self.logger.info('made a plot of Rg and I(0) vs file number')
#            print 'before'
#            time.sleep(5)
#            print 'after'
        else:
            ax = fig.add_subplot(111, autoscale_on=True)
            count = 0
            if len(self.imagedict.keys()) > 0:
                for file in self.imagedict.keys():
                    self.logger.info('Adding '+str(file)+' to plot')
                    from_point=self.imagedict[file]['I'].index(max(self.imagedict[file]['I']))

                    if self.jobtype == 'fit':
                        ax.plot(self.imagedict[file]['Q'][from_point:],self.imagedict[file]['I'][from_point:],color='0.5', ls='.', marker='.', markersize=2)
                        ax.plot(self.imagedict[file]['Q'][from_point:],self.imagedict[file]['E'][from_point:],'k-')
                    else:
                        ax.plot(self.imagedict[file]['Q'],self.imagedict[file]['I'],linetype_array[count])
                    count += 1

                if self.jobtype == 'dat' or self.jobtype == 'fit':
                    ax.set_yscale('log')
                    ax.set_xscale('log')
                    ax.set_xlabel('Q ($\AA^{-1}$)')
                    ax.set_ylabel('Intensity')
                elif self.jobtype == 'out':
                    ax.set_ylim(ymin=0)
                    ax.set_xlabel('Q ($\AA$)', fontsize=16)
                    ax.set_ylabel('P(R)', fontsize=16)
                
            plt.show()
            
    def saveImage(self):
        ######################
        #save the output file#
        ######################
        outfilename = self.output_directory+'/'+self.imagedict.keys()[0]+'.pdf'
        pp = PdfPages(outfilename)
        pp.savefig(fig)
        pp.close()
        self.logger.info('saved file '+outfilename)

    def PlotLowQ(self):
        self.needs_blanked = True
	self.logger.info('the HPLC dat files need to be blanked')

        filenumber_array = []
        intensity_array = []
        for file in sorted(self.imagedict.keys()):
            try:
                #filenumber = int(file.split('_')[-1])
                filenumber = int(re.split('[-._]',file)[-1])
                intensity = sum(self.imagedict[file]['I'][10:20])
                filenumber_array.append(filenumber)
                intensity_array.append(intensity)
            except:
                self.logger.error('could not parse low Q data for '+str(file))
        lowi = min(intensity_array)
        highi = max(intensity_array)

        fig = plt.figure()
        ax = fig.add_subplot(111, autoscale_on=True)
        ax.plot(filenumber_array,intensity_array,color='black',lw=2)
        ax.set_xlabel('file number')
        ax.set_ylabel('Average Intensity at low Q')
        plt.interactive(True)
        plt.show()
        self.logger.info('made a plot of the low Q data')

        #CHOOSE THE LOWER AVERAGING WINDOW
        plt.title('select the lower edge of the first averaging window',fontsize=16)
        plt.show()
        bottom_lower = round(plt.ginput(1)[0][0])
        plt.vlines(bottom_lower,lowi,highi,colors='green',linestyle='solid')
        plt.title('select the upper edge of the first averaging window',fontsize=16)
        top_lower = round(plt.ginput(1)[0][0])
        plt.vlines(top_lower,lowi,highi,colors='green',linestyle='solid')
        lower_average = round((top_lower + bottom_lower)/2)
        self.logger.info('User selected a window below the peak from '+str(bottom_lower)+'-'+str(top_lower))

        #CHOOSE THE UPPER AVERAGING WINDOW
        plt.title('select the lower edge of the second averaging window',fontsize=16)
        bottom_upper = round(plt.ginput(1)[0][0])
        plt.vlines(bottom_upper,lowi,highi,colors='red',linestyle='solid')
        plt.title('select the upper edge of the second averaging window',fontsize=16)
        top_upper = round(plt.ginput(1)[0][0])
        plt.vlines(top_upper,lowi,highi,colors='red',linestyle='solid')
        upper_average = round((top_upper + bottom_upper)/2)
        self.logger.info('User selected a window above the peak from '+str(bottom_upper)+'-'+str(top_upper))

        #CHOOSE THE REGION TO BE BLANKED
        plt.title('select the lower edge of the region to be blanked',fontsize=16)
        bottom_region = round(plt.ginput(1)[0][0])
        plt.vlines(bottom_region,lowi,highi,colors='blue',linestyle='solid')
        plt.title('select the upper edge of the region to be blanked',fontsize=16)
        top_region = round(plt.ginput(1)[0][0])
        plt.vlines(top_region,lowi,highi,colors='blue',linestyle='solid')
        self.logger.info('User selected a window of files to be blanked from '+str(bottom_region)+'-'+str(top_region))
        plt.title('Click again when you are ready to do the blanking',fontsize=16)
        plt.waitforbuttonpress()

        #AVERAGE THE FIRST WINDOW
        window_dict = {}

        for file in sorted(self.imagedict.keys()):
            #filenumber = int(file.split('_')[-1])
            filenumber = int(re.split('[-._]',file)[-1])
            try:
                if int(filenumber) >= bottom_lower and int(filenumber) <= top_lower:
                    for q in self.imagedict[file]['Q']:
                        index = self.imagedict[file]['Q'].index(q)
                        if q in window_dict.keys():
                            window_dict[q]['I'].append(self.imagedict[file]['I'][index])
                            window_dict[q]['E'].append(self.imagedict[file]['E'][index])
                        else:
                            window_dict[q] = {'I': [self.imagedict[file]['I'][index]], 'E': [self.imagedict[file]['E'][index]]}
            except:
                self.logger.error('could not parse file '+str(file)+' for the lower averaging window')
        qdata = []
        idata = []
        edata = []
        
        for q in sorted(window_dict.keys()):
            qdata.append(q)
            idata.append(sum(window_dict[q]['I'])/len(window_dict[q]['I']))
            edata.append(sum(window_dict[q]['E'])/len(window_dict[q]['E']))
            
        self.output_images['lower_blank'] = {'Q': qdata, 'I': idata, 'E': edata}


        #AVERAGE THE SECOND WINDOW
        window_dict = {}
        for file in sorted(self.imagedict.keys()):
            #filenumber = int(file.split('_')[-1])
            filenumber = int(re.split('[-._]',file)[-1])
            try:
                if int(filenumber) >= bottom_upper and int(filenumber) <= top_upper:
                    for q in self.imagedict[file]['Q']:
                        index = self.imagedict[file]['Q'].index(q)
                        if q in window_dict.keys():
                            window_dict[q]['I'].append(self.imagedict[file]['I'][index])
                            window_dict[q]['E'].append(self.imagedict[file]['E'][index])
                        else:
                            window_dict[q] = {'I': [self.imagedict[file]['I'][index]], 'E': [self.imagedict[file]['E'][index]]}
            except:
                self.logger.error('could not parse file '+str(file)+' for the upper averaging window')


        qdata = []
        idata = []
        edata = []
        
        for q in sorted(window_dict.keys()):
            qdata.append(q)
            idata.append(sum(window_dict[q]['I'])/len(window_dict[q]['I']))
            edata.append(sum(window_dict[q]['E'])/len(window_dict[q]['E']))
            
        self.output_images['upper_blank'] = {'Q': qdata, 'I': idata, 'E': edata}

        #SUBTRACT A LINEAR INTERPOLATION OF THE TWO AVERAGING WINDOWS FROM EACH FILE IN THE PEAK WINDOW
        #LINEAR REGRESSION ON TWO WINDOWS
        warnings.simplefilter('ignore', RuntimeWarning)
        upperslope, upperintercept, r_value, p_value, std_err = stats.linregress([lower_average,upper_average],[0,1])
        lowerslope, lowerintercept, r_value, p_value, std_err = stats.linregress([lower_average,upper_average],[1,0])
        for file in sorted(self.imagedict.keys()):
            #filenumber = int(file.split('_')[-1])
            filenumber = int(re.split('[-._]',file)[-1])
            if int(filenumber) >= bottom_region and int(filenumber) <= top_region:
                fraction_upper = (filenumber * upperslope) + upperintercept
                fraction_lower = (filenumber * lowerslope) + lowerintercept
                qdata = []
                idata = []
                edata = []
                for q in self.imagedict[file]['Q']:
                    index = self.imagedict[file]['Q'].index(q)
                    if q in self.output_images['lower_blank']['Q'] and q in self.output_images['upper_blank']['Q']:
                        i = self.imagedict[file]['I'][index] - (( self.output_images['lower_blank']['I'][index] * fraction_lower ) + ( self.output_images['upper_blank']['I'][index] * fraction_upper ))
                        e = self.imagedict[file]['E'][index] + (( self.output_images['lower_blank']['E'][index] * fraction_lower ) + ( self.output_images['upper_blank']['E'][index] * fraction_upper ))
                        qdata.append(q)
                        idata.append(i)
                        edata.append(e)
                    else:
                        self.logger.error('Q value '+str(q)+' was found in file '+str(file)+' but was missing in at least one of the averaged windows')
                        self.logger.error('it is possible you might have different Q ranges within your files')
                self.output_images[str(file)+'_blanked'] = {'Q': qdata, 'I': idata, 'E': edata}
                self.logger.info('Blanked file '+str(file))
                                    

    def OutputDat(self):
        try:
            if self.needs_blanked:
	        self.logger.info('Writing out the blanked files')
	    else:
                return
        except:
            pass

        for file in self.output_images.keys():
            string_list = []
            string_list.append("%-15s %-18s %-15s" % ("Q(A-1)","I(au)","Error"))
            for q in self.output_images[file]['Q']:
                if q > 0:
                    index = self.output_images[file]['Q'].index(q)
                    q = self.output_images[file]['Q'][index]
                    i = self.output_images[file]['I'][index]
                    e = self.output_images[file]['E'][index]
                    string_list.append("%-15s %-18s %-15s" % (q,i,e))

            outfile = open(self.output_directory+'/'+str(file)+'.dat', 'w')
            outfile.write('\n'.join(string_list))
            outfile.close()
            self.logger.info('wrote averaged file to: '+self.output_directory+'/'+str(file)+'.dat')
    def GetRgIo(self):
        if self.needs_blanked:
	    filelist = self.output_images
	else:
	    filelist = self.imagelist

        self.filenumbers = []
        self.rgs = []
        self.ios = []
        for file in sorted(filelist):
            try:
                filenumber = int(re.split('[_\-./]', file)[-2])
                self.logger.info('running autorg on file '+str(file))
                if self.needs_blanked:
                    command = 'autorg '+str(self.output_directory)+'/'+str(file)+'.dat'
		else:
		    command = 'autorg '+str(file)
                p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
                data = p.stdout.readlines()
                rg = 0
                io = 0

                for line in data:
                    if line.split()[0] == 'Rg':
                        rg = line.split()[2]
                        rg_error = line.split()[4]
                    if line.split()[0] == 'I(0)':
                        io = line.split()[2]
                        io_error = line.split()[4]
                    if line.split()[0] == 'Quality:':
                        quality = line.split()[1].replace('%', '')
                self.filenumbers.append(filenumber)
                self.rgs.append(rg)
                self.ios.append(io)
            except:
                if file == 'lower_blank' or file == 'upper_blank':
                    pass
                else:
                    self.logger.error('Could not run autoRg on the file '+str(file))

    def Average(self):
        filelist = self.imagedict.keys()
        index_list = []
        for file in filelist:
            index_list.append(int(file.split('_')[-1]))
        output_name = '_'.join(filelist[0].split('_')[0:-1])+'_'+str(min(index_list)).zfill(len(filelist[0].split('_')[-1]))+'-'+str(max(index_list)).zfill(len(filelist[0].split('_')[-1]))



        if len(filelist) < 2:
            self.logger.error('Averaging requires more than one input file')
            return
        qdata = self.imagedict[self.imagedict.keys()[0]]['Q']
        idata = []
        edata = []

        for q in qdata:
            i = []
            e = []
            for file in filelist:
                if q in self.imagedict[file]['Q']:
                    index = self.imagedict[file]['Q'].index(q)
                    i.append(self.imagedict[file]['I'][index])
                    e.append(self.imagedict[file]['E'][index])
            idata.append(sum(i)/len(i))
            edata.append(sum(e)/len(e))
        
        self.output_images[output_name] = {'Q': qdata, 'I': idata, 'E': edata}

    def LogBin(self):

        for file in self.imagedict.keys():
            self.logger.info('Log binning '+str(file))
            noBins = 100
            n = [0, ]
            while 0 in n:
                logBins = numpy.logspace(numpy.log10(min(self.imagedict[file]['Q'])),numpy.log10(max(self.imagedict[file]['Q'])),num=noBins)            
                n, _ = numpy.histogram(self.imagedict[file]['Q'],bins=logBins)
                noBins -= 5
            binnedInts, _ = numpy.histogram(self.imagedict[file]['Q'],bins=logBins, weights=self.imagedict[file]['I'])
            binnedErrs, _ = numpy.histogram(self.imagedict[file]['Q'],bins=logBins, weights=self.imagedict[file]['E'])

            meanInts = binnedInts / n
            meanErrs = binnedErrs / n / n
            logBins = [((x+logBins[i+1])/2.0) for i, x in enumerate(logBins[:-1])]
            self.output_images[str(file)+'_log'] = {'Q': logBins, 'I': meanInts.tolist(), 'E': meanErrs.tolist()}

    def BoxCar(self, window=None):
        if window == None:
            window = 5
            self.logger.info('No window size specified for BoxCar averaging, will use 5')
        elif type(window) != type(1):
            self.logger.error('window size for boxcar averaging is not an integer, will use 5 instead! Help me, help you...')
            window = 5
        elif window%2==0:
            window = window + 1
            self.logger.error('window size needs to be an odd number, added 1 to it')
        else:
            window = window
        self.logger.info('A window of '+str(window)+' will be used for BoxCar averaging')

        
        for file in sorted(self.imagedict.keys()):
            index = sorted(self.imagedict.keys()).index(file)
            half_window = ( window - 1 ) / 2
            while min(range(index - (half_window), index + (half_window + 1))) < 0 or max(range(index - (half_window), index + (half_window + 1))) > len(self.imagedict.keys())-1:
                half_window -= 1
            averaging_range = range(index - (half_window), index + (half_window + 1))
            filenames = []
            for item in averaging_range:
                filenames.append(sorted(self.imagedict.keys())[item])
            output_name = '_'.join(re.split('[._]',filenames[0])[0:-1])+'_'+re.split('[._]',filenames[0])[-1]+'-'+re.split('[._]',filenames[-1])[-1]

            qdata = self.imagedict[filenames[0]]['Q']
            idata = []
            edata = []

            for q in qdata:
                i = []
                e = []
                for filename in filenames:
                    if q in self.imagedict[filename]['Q']:
                        index = self.imagedict[filename]['Q'].index(q)
                        i.append(self.imagedict[filename]['I'][index])
                        e.append(self.imagedict[filename]['E'][index])
                idata.append(sum(i)/len(i))
                edata.append(sum(e)/(len(e)*numpy.sqrt(len(e))))#average divided by sqrt of number of samples
        
            self.output_images[output_name] = {'Q': qdata, 'I': idata, 'E': edata}

                    
        
if __name__ == '__main__':
    plot = SaxsPlot()
    plot.GetImages()
    if plot.jobtype == 'col':
       plot.ParseDat()        
       plot.PlotLowQ()
       plot.OutputDat()
       plot.GetRgIo()
       app.Yield()
       plot.Plot()
    elif plot.jobtype == 'dat':
        plot.ParseDat()
	plot.ScaleIntensities()
        plot.Plot()
    elif plot.jobtype == 'out':
        plot.ParseOut()
	plot.ScaleIntensities()
        plot.Plot()
    elif plot.jobtype == 'fit':
        plot.ParseFit()
        plot.Plot()
    else:
        plot.logger.error('Undefined job type, exitting')
        
