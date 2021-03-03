#!/usr/local/bin/python3
import logging
import os.path
import seaborn
import matplotlib.pyplot as plt

class PlotFit():
    """Plot SAXS fit
    
    The PlotFit class facilitates plotting of the fit between
    experimental data and a model.

    
    """
    
    '''
    Constructor
    '''
    def __init__(self):
        ###start a log file
        self.logger = logging.getLogger('PlotFit')
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(module)s: %(message)s',"[%Y-%m-%d %H:%M:%S]")
        streamhandler = logging.StreamHandler()
        streamhandler.setFormatter(formatter)
        self.logger.addHandler(streamhandler)
        self.logger.info('Starting a new PlotFit job')

        ###SOME PARAMETERS
        self.fit_files = []
        self.column_q = 0
        self.column_obs = 1
        self.column_cal = 2
        self.column_err = 3
        self.units = 'angstrom'
        self.log_y = True
        self.log_x = False
        
        
    def addFitFile(self, input_file=None):
        data_array = {'Q': [], 'OBS': [], 'CAL': [], 'ERR': [], 'file': str(input_file), 'n_points': 0}
        if input_file:
            if os.path.isfile(input_file):
                with open(input_file, 'r') as f:
                    for line in f.readlines():
                        line = line.split()
                        try:
                            q = float(line[self.column_q])
                            if self.column_obs == None:
                                obs = 0.0
                            else:
                                obs = float(line[self.column_obs])
                            if self.column_cal == None:
                                exp = 0.0
                            else:
                                exp = float(line[self.column_cal])
                            if self.column_err == None:
                                err = 0.0
                            else:
                                err = float(line[self.column_err])
                            data_array['Q'].append(q)
                            data_array['OBS'].append(obs)
                            data_array['CAL'].append(exp)
                            data_array['ERR'].append(err)
                            data_array['n_points']+=1 
                        except:
                            pass
                
                if len(data_array['Q']) > 0:
                    self.logger.info(f'{str(input_file)} parsed with {str(data_array["n_points"])} data points')
                    self.fit_files.append(data_array)
                    return True
                else:
                    self.logger.error(f'str({input_file} parse yielded no data points')
            else:
                self.logger.error(f'{str(input_file)} does not exist')
                return False
        else:
            self.logger.error('addFitFile function needs a filename as input')
            return False


    def plot(self):
        seaborn.set(style='whitegrid')
        seaborn.residplot(self.fit_files[0]['Q'], self.fit_files[0]['OBS'], color='g')
        seaborn.lineplot(self.fit_files[0]['Q'], self.fit_files[0]['CAL'], color='k')
        plt.show()
        
if __name__ == '__main__':
    job = PlotFit()
    job.addFitFile('/Users/nathan/Documents/SAXS/HEENA_2019/NMA_test/pool/bsa_nma_m22_bsa_2_average_average.dat')
    job.plot()
