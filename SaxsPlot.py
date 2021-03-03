#!/usr/local/bin/python3

import wx, wx.lib.dialogs, logging, re, os, sys, time
import wx.lib.agw.genericmessagedialog as GMD
import matplotlib.pyplot as plt
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
        if len(sys.argv) > 1:
            self.imagelist = []
            for arg in sys.argv[1:]:
                if os.path.isfile(arg):
                    self.imagelist.append(os.path.abspath(arg))
        else:
        ########################################################
        #USER DID NOT SUPPLY ARGS: PROMPT WITH FILE OPEN DIALOG#
        ########################################################
            wildcard = "Dat files (*.dat)|*.dat|" \
                   "P(r) files (*.out)|*.out|" \
                   "fit files (*.fit)|*.fit|" \
                   "fir files (*.fir)|*.fir|" \
                       "All files (*.*)|*.*"
            directory='/Users/nathan/Documents/SAXS/'        
            self.imagelist = []
            count = 0

            while len(self.imagelist) == 0 and count < 3:
                response = wx.lib.dialogs.openFileDialog(title='Choose SAXS files (.out, .fit/.fir or .dat)', directory=directory, wildcard=wildcard, style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR )
                if not response.paths:
                    count += 1
                elif len(response.paths) == 0:
                    self.logger.info('You did not select any images, try again.')
                    count += 1
                else:
                    self.imagelist = response.paths

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
                #index2 = [i for i, item in enumerate(saxsdata[index:]) if re.search(' *Reciprocal space: Rg =.*', item)][0] + index - 1
                #saxsdata = saxsdata[index:index2]
                saxsdata = saxsdata[index:]
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
                self.logger.error(f'Outfile {fileroot}.out does not seem to have the expected format for an outfile')
                

    def ScaleIntensities(self):
        title='Scale data?'
        message='Would you like to scale the intensities so that the various plots overlay?'
        dlg = GMD.GenericMessageDialog(None, message, title, agwStyle=wx.ICON_INFORMATION | wx.YES_NO)
        response = dlg.ShowModal()
        dlg.Destroy()
        if response == wx.ID_YES:
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
        else:
            self.logger.info('Will proceed with unscaled data')



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

        ######################
        #save the output file#
        ######################
        title='Save PDF?'
        message=f'Would you like to save the plot as {self.output_directory}/{self.imagedict.keys()[0]}.pdf'
        dlg = GMD.GenericMessageDialog(None, message, title, agwStyle=wx.ICON_INFORMATION | wx.YES_NO)
        response = dlg.ShowModal()
        dlg.Destroy()
        if response == wx.ID_YES:
            outfilename = self.output_directory+'/'+self.imagedict.keys()[0]+'.pdf'
            pp = PdfPages(outfilename)
            pp.savefig(fig)
            pp.close()
            self.logger.info('saved file '+outfilename)
        else:
            self.logger.info('Did not output an image file')

        
if __name__ == '__main__':
    app = wx.App()
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
        
