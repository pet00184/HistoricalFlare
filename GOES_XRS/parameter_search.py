import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy import stats as st
import math


class ParameterSearch:
    
    flare_fits = 'GOES_XRS_historical.fits'
    
    def __init__(self, parameter_array, save_string, directory):
        '''Saves .fits file data to Astropy Table structure (works similarly to regular .fits, but also lets you
        parse the data by rows.)
        '''
        fitsfile = fits.open(self.flare_fits)
        self.data = Table(fitsfile[1].data)[:]
        self.header = fitsfile[1].header
        self.param_grid = np.array(parameter_array)
        self.savestring = save_string
        self.directory = directory
        self.calculated_flarelist = [] #has format of [flare #, flare ID, max foxsi, mean foxsi, max hic, mean hic] for each tuple
        self.launches_df = pd.DataFrame(columns=('Flare_Number', 'Flare_ID', 'Trigger_Time', 'Max_FOXSI', 'Mean_FOXSI', 'Max_HiC', 'Mean_HiC', 'Max_FOXSI_C5', 
                        'Mean_FOXSI_C5', 'Max_HiC_C5', 'Mean_HiC_C5', 'Flare_C5', 'Flare_Class', 'Flare_Max_Flux', 'Peak_Time', 'Start_to_Peak', 'Background_Flux'))
        
    def loop_through_parameters(self, arrays_to_check):
        ''' Loops through each parameter, and performes launch analysis on each flare. This is the function you will
        call for each runthrough!!
        '''
        for parameter in self.param_grid:
            print(f'starting parameter search for {parameter}')
            self.loop_through_flares(arrays_to_check, parameter)
            print('now trying to perform post-loop functions')
            self.perform_postloop_functions()
            self.save_DataFrame(parameter)
            print('dataframe saved!')
            self.calculated_flarelist = []
            self.launches_df = self.launches_df.iloc[0:0]
            

################ Flare Loop Functions ############################################################################   
    def loop_through_flares(self, arrays_to_check, parameter):
       ''' Loops through each flare in the calculated array that is being checked. For simplest example, array is just
       self.data['xrsb'], but I'll add more functions to calculate derivatives, temps, etc. 
   
       Input: 
       arrays_to_check: array of flares to loop through, when checking of parameter was met. (For the simple xrsb example
           this is self.data['xrsb])
       parameter: the parameter currently being used (in simple example, this is xrsb flux level)
       '''
       for i, flare in enumerate(arrays_to_check):
           self.flareloop_check_if_value_surpassed(flare, parameter, i)
           if self.triggered_bool: 
               self.calculate_observed_xrsb(i)         

    def flareloop_check_if_value_surpassed(self, arrays, parameters, i):
        ''' Process used to loop through flares when there is only a value being checked, and whether the curve
        surpasses that value. This will start as a blank slate to be used for xrsb and xrsa levels, and could also
        be used once a new array of the derivative is calculated, temp, etc. at each point.
    
        Input: 
        array = list of arrays (flare) to be checked (xrsa, xrsb or a computed temp/derivative etc.)
        parametesr = list of values that if surpassed triggers a launch.
        ***these are now lists for multiple parameters that will be zipped!!
    
        Returns: 
        triggered_bool = True if this flare triggers a launch, otherwise is False.
        indeces of the trigger, foxsi obs start/end and hic obs start/end to be used for computing observed flux
        '''
        self.triggered_bool = False
        df = pd.DataFrame()
        if isinstance(parameters, np.float64):
            triggered_check = np.where(arrays > parameters)[0]
        elif isinstance(parameters, np.int64):
            triggered_check = np.where(arrays > parameters)[0]
        else:
            for arr, p in zip(arrays, parameters):
                df[f'param {p}'] = np.array(arr) >= p
            truth_df = df.all(1)
            triggered_check = np.where(truth_df == True)[0]
        #triggered_check = np.where(array >= parameter)[0]
        if not len(triggered_check)==0:
            self.triggered_bool = True
            self.trigger_index = triggered_check[0]
            self.foxsi_obs_start = self.trigger_index + 3 + 3 + 2 #latency + launch prep + launch time
            self.foxsi_obs_end = self.foxsi_obs_start + 6
            self.hic_obs_start = self.foxsi_obs_start + 2
            self.hic_obs_end = self.hic_obs_start + 6

    def calculate_observed_xrsb(self, i):
        ''' Slices out the FOXSI and HiC observation windows of the current flare # (i), and calculated the max and mean 
        observed fluxes. Only do this if triggered bool is true! (will add in when doing the loop)
    
        Appends [i, flare ID, foxsi max, foxsi mean, hic max, hic mean] to the flarelist so that the tuples can be zipped and 
        moved to a pandas DF after all the flares are looped through.
        '''
        foxsi_obs_xrsb = self.data['xrsb'][i][self.foxsi_obs_start:self.foxsi_obs_end]
        if len(foxsi_obs_xrsb)==0: #dealing with the observation range being outside the flare (probably the next flare)
            foxsi_max_observed = math.nan
            foxsi_mean_observed = math.nan
        else:
            foxsi_max_observed = np.max(foxsi_obs_xrsb)
            foxsi_mean_observed = np.mean(foxsi_obs_xrsb)
        hic_obs_xrsb = self.data['xrsb'][i][self.hic_obs_start:self.hic_obs_end]
        if len(hic_obs_xrsb)==0:
            hic_max_observed = math.nan
            hic_mean_observed = math.nan
        else:    
            hic_max_observed = np.max(hic_obs_xrsb)
            hic_mean_observed = np.mean(hic_obs_xrsb)
        flare_ID = self.data['flare ID'][i]
        trigger_time = self.data['time'][i][self.trigger_index]
        self.calculated_flarelist.append([i, flare_ID, trigger_time, foxsi_max_observed, foxsi_mean_observed, hic_max_observed, hic_mean_observed])

################ Post-Loop Functions ############################################################################      
    def perform_postloop_functions(self):
        ''' Once completed, a finished DataFrame should have info saved for all launches.
        '''
        self.save_flarelist_to_df()
        self.save_fitsinfo_to_df()
        self.calculate_c5_bool()
        self.drop_na()
    
    def save_flarelist_to_df(self):
        ''' This is done outside of the loop (for all iterations). Saves calculated values to DataFrame for each flare.
        '''
        self.launches_df[['Flare_Number', 'Flare_ID', 'Trigger_Time', 'Max_FOXSI', 'Mean_FOXSI', 'Max_HiC', 'Mean_HiC']] = self.calculated_flarelist  
    
    def save_fitsinfo_to_df(self):
        ''' Saves the flare class, peak flux, start to peak time, and if flare is above C5 bool info from the FITS file
        using the flare number. 
        '''
        for f, flare_id in enumerate(self.launches_df['Flare_ID']):
            launched_flare = np.where(flare_id == self.data['flare ID'])[0][0]
            self.launches_df.loc[f, 'Flare_Class'] = self.data['class'][launched_flare]
            self.launches_df.loc[f, 'Flare_Max_Flux'] = self.data['peak flux'][launched_flare]
            self.launches_df.loc[f, 'Start_to_Peak'] = self.data['start to peak time'][launched_flare]
            self.launches_df.loc[f, 'Flare_C5'] = self.data['above C5'][launched_flare]
            self.launches_df.loc[f, 'Flare_C5_10min'] = self.data['above C5 10min'][launched_flare]
            self.launches_df.loc[f, 'Background_Flux'] = self.data['background flux'][launched_flare]
            self.launches_df.loc[f, 'Peak_Time'] = self.data['UTC peak time'][launched_flare]
    
    def calculate_c5_bool(self):
        ''' Saves True/False boolean results, for if the Flare and the observed flux is above C5. (Done for the flare 
        itself, max/mean FOXSI and max/mean HiC)
        '''    
        self.launches_df['Max_FOXSI_C5'] = self.launches_df['Max_FOXSI'] > 5e-6
        self.launches_df['Mean_FOXSI_C5'] = self.launches_df['Mean_FOXSI'] > 5e-6
        self.launches_df['Max_HiC_C5'] = self.launches_df['Max_HiC'] > 5e-6
        self.launches_df['Mean_HiC_C5'] = self.launches_df['Mean_HiC'] > 5e-6
        
    def drop_na(self):
        print('before drop NA')
        print(len(self.launches_df['Flare_ID']))
        self.launches_df = self.launches_df.dropna()
        print('after drop NA')
        print(len(self.launches_df['Flare_ID']))
        
    
    def save_DataFrame(self, parameter):
        self.launches_df.to_csv(f'{self.directory}/{parameter}_{self.savestring}_results.csv')
    


        
    
        
        
        