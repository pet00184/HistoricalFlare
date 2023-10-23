import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy import stats as st
import math
import os


class ParameterSearch:
    
    flare_fits = 'GOES_XRS_historical.fits'
    
    def __init__(self, parameter_array, directory):
        '''Saves .fits file data to Astropy Table structure (works similarly to regular .fits, but also lets you
        parse the data by rows.)
        '''
        fitsfile = fits.open(self.flare_fits)
        self.data = Table(fitsfile[1].data)[:]
        self.header = fitsfile[1].header
        self.param_grid = np.array(parameter_array)
        self.directory = directory
        self.calculated_flarelist = [] #has format of [flare #, flare ID, max foxsi, mean foxsi, max hic, mean hic] for each tuple
        self.launches_df = pd.DataFrame(columns=('Flare_Number', 'Flare_ID', 'Trigger_Time', 'Cancelled?', 'Max_FOXSI', 'Mean_FOXSI', 'Max_HiC', 'Mean_HiC', 'Max_FOXSI_and_HiC_C5', 'Max_FOXSI_C5', 
                        'Mean_FOXSI_C5', 'Max_HiC_C5', 'Mean_HiC_C5', 'Flare_C5', 'Flare_Class', 'Flare_Max_Flux', 'Peak_Time', 'Start_to_Peak', 'Background_Flux'))
        
        self.param_combo_results = pd.DataFrame(columns=('Param Combo', 'TN', 'TN_canc', 'TP', 'TP_noc5_obs', 'FN', 'FN_canc', 'FP_c5_noobs', 
                        'FP_noc5_noobs', 'Precision', 'Recall', 'Trigger-to-Launch'))
                        
        if not os.path.exists('Launches'):
            os.makedirs(f'{self.directory}/Launches')
        if not os.path.exists('GoodEggParams'):
            os.makedirs(f'{self.directory}/GoodEggParams')
        if not os.path.exists('BadEggParams'):
            os.makedirs(f'{self.directory}/BadEggParams')
        
    def loop_through_parameters(self, arrays_to_check):
        ''' Loops through each parameter, and performes launch analysis on each flare. This is the function you will
        call for each runthrough!!
        '''
        for j, parameter in enumerate(self.param_grid):
            print(f'starting parameter search for {parameter}')
            self.loop_through_flares(arrays_to_check, parameter)
            print('now trying to perform post-loop functions')
            self.perform_postloop_functions(parameter, j)
            self.save_launch_DataFrame(parameter)
            self.calculated_flarelist = []
            self.launches_df = self.launches_df.iloc[0:0]
        self.save_and_sort_comboDataFrames()
            

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
           # print(flare)
           # print(flare[0])
           self.flareloop_check_if_value_surpassed(flare[0], parameter, i)
           if self.triggered_bool: 
               self.calculate_observed_xrsb_and_cancellation(i)         

    def flareloop_check_if_value_surpassed(self, arrays, parameters, i):
        ''' Process used to loop through flares when there is only a value being checked, and whether the curve
        surpasses that value. This will start as a blank slate to be used for xrsb and xrsa levels, and could also
        be used once a new array of the derivative is calculated, temp, etc. at each point.
        
        ADDING CANCELLATION: I am still saving what we "would have" observed if we didn't cancel, and just doing 
        a bool for cancelled. This way, we can still get some information on what we are cancelling on. We are doing a 
        simple cancellation of only cancelling if the xrsa flux is decreasing during the pre-launch window.
    
        Input: 
        array = list of arrays (flare) to be checked (xrsa, xrsb or a computed temp/derivative etc.)
        parametesr = list of values that if surpassed triggers a launch.
        ***these are now lists for multiple parameters that will be zipped!!
    
        Returns: 
        triggered_bool = True if this flare triggers a launch, otherwise is False.
        indeces of the trigger, foxsi obs start/end and hic obs start/end to be used for computing observed flux
        CANCELLATION bool, so that we know if we would have cancelled the launch or not.
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
        if not len(triggered_check)==0:
            self.triggered_bool = True
            self.trigger_index = triggered_check[0]
            self.foxsi_obs_start = self.trigger_index + 3 + 3 + 2 #latency + launch prep + launch time
            self.foxsi_obs_end = self.foxsi_obs_start + 6
            self.hic_obs_start = self.foxsi_obs_start + 2
            self.hic_obs_end = self.hic_obs_start + 6
            

    def calculate_observed_xrsb_and_cancellation(self, i):
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
        if (self.trigger_index + 3) < len(self.data['xrsa'][i]):
            cancellation_bool = (self.data['xrsa'][i][self.trigger_index + 3] - self.data['xrsa'][i][self.trigger_index]) < 0
        else:
            cancellation_bool = math.nan
        self.calculated_flarelist.append([i, flare_ID, cancellation_bool, trigger_time, foxsi_max_observed, foxsi_mean_observed, hic_max_observed, hic_mean_observed])

################ Post-Loop Functions ############################################################################      
    def perform_postloop_functions(self, parameter, j):
        ''' Once completed, a finished DataFrame should have info saved for all launches.
        '''
        self.save_flarelist_to_df()
        self.save_fitsinfo_to_df()
        self.calculate_c5_bool()
        self.drop_na()
        self.save_launch_and_obs_arrays()
        self.save_cf_input(j)
        self.save_precision_recall_ttlscore(parameter, j)
    
    def save_flarelist_to_df(self):
        ''' This is done outside of the loop (for all iterations). Saves calculated values to DataFrame for each flare.
        '''
        self.launches_df[['Flare_Number', 'Flare_ID', 'Cancelled?', 'Trigger_Time', 'Max_FOXSI', 'Mean_FOXSI', 'Max_HiC', 'Mean_HiC']] = self.calculated_flarelist  
    
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
        self.launches_df['Max_FOXSI_and_HiC_C5'] = (self.launches_df['Max_FOXSI'] > 5e-6) & (self.launches_df['Max_HiC'] > 5e-6)
        
    def drop_na(self):
        ''' Drops rows with Nan for observation times. This helps get rid of double counting, since sometimes the 
        next flare is triggered on the previous flare ID.
        '''
        print('before drop NA')
        print(len(self.launches_df['Flare_ID']))
        self.launches_df = self.launches_df.dropna()
        print('after drop NA')
        print(len(self.launches_df['Flare_ID']))
        
    def save_launch_and_obs_arrays(self):
        '''This saves arrays of true/false for the baseline truth (was the flare above C5), the trigger prediction
        (there was a trigger), the launch prediction (was there a launch) and the launch/obs prediction 
        (launch and obs) and cancellation. From this, we can save the scores needed
        for the confusion matrices. 
        '''
        #defining the array that determines if a flare is above C5
        self.all_abovec5 = self.data['above C5'] #all flares above C5- this is the baseline "truth"
        
        #make an array that determines whether there was a trigger
        all_flareID = self.data['flare ID']
        flareID_triggers = np.array(self.launches_df['Flare_ID'])
        trigger_indx = []
        for flare in flareID_triggers:
            indx = np.where(flare==all_flareID)[0][0]
            trigger_indx.append(indx)
        self.trigger_array = np.array([False]*len(all_flareID))
        np.put(self.trigger_array, trigger_indx, True)
        
        #make an array that determines whether there was a launch
        all_flareID = self.data['flare ID']
        flareID_launches = np.array(self.launches_df['Flare_ID'][self.launches_df['Cancelled?']==False])
        launch_indx = []
        for flare in flareID_launches:
            indx = np.where(flare==all_flareID)[0][0]
            launch_indx.append(indx)
        self.launch_array = np.array([False]*len(all_flareID))
        np.put(self.launch_array, launch_indx, True)
        
        #make an array that defines where flux above C5 was observed by both FOXSI and HIC
        flareID_observations = np.array(self.launches_df['Flare_ID'][(self.launches_df['Cancelled?']==False) & (self.launches_df['Max_FOXSI_and_HiC_C5']==True)])
        obs_indx = []
        for flare in flareID_observations:
            indx = np.where(flare==all_flareID)[0][0]
            obs_indx.append(indx)
        self.obs_array = np.array([False]*len(all_flareID))
        np.put(self.obs_array, obs_indx, True)
        
        #make an array that defines where the launch was cancelled
        flareID_cancellations = np.array(self.launches_df['Flare_ID'][self.launches_df['Cancelled?']==True])
        canc_indx = []
        for flare in flareID_cancellations:
            indx = np.where(flare==all_flareID)[0][0]
            canc_indx.append(indx)
        self.canc_array = np.array([False]*len(all_flareID))
        np.put(self.canc_array, canc_indx, True)
        
    def save_cf_input(self, j):
        ''' This is where we save the values for each box of the 4x2 confusion matrix. I'll now save all of this as
        columns in the launch dataframe. More work on the front end, but plotting will be easier!
        '''
        #self.TN = len(np.where((self.all_abovec5==False) & (self.trigger_array==False))[0])
        self.param_combo_results.loc[j, 'TN'] = len(np.where((self.all_abovec5==False) & (self.trigger_array==False))[0])
        self.param_combo_results.loc[j, 'TN_canc'] = len(np.where((self.all_abovec5==False) & (self.canc_array==True))[0])
        self.param_combo_results.loc[j, 'TP'] = len(np.where((self.all_abovec5==True) & (self.launch_array==True) & (self.obs_array==True))[0])
        self.param_combo_results.loc[j, 'TP_noc5_obs'] = len(np.where((self.all_abovec5==False) & (self.launch_array==True) & (self.obs_array==True))[0])
        self.param_combo_results.loc[j, 'FN'] = len(np.where((self.all_abovec5==True) & (self.trigger_array==False))[0])
        self.param_combo_results.loc[j, 'FN_canc'] = len(np.where((self.all_abovec5==True) & (self.canc_array==True))[0])
        self.param_combo_results.loc[j, 'FP_c5_noobs'] = len(np.where((self.all_abovec5==True) & (self.launch_array==True) & (self.obs_array==False))[0])
        self.param_combo_results.loc[j, 'FP_noc5_noobs'] = len(np.where((self.all_abovec5==False) & (self.launch_array==True) & (self.obs_array==False))[0])
        
    def save_precision_recall_ttlscore(self, parameter, j):
        ''' Saves the precision, recal and "trigger-to-launch" scores for each parameter. These will be used for sorting
        the parameters into the "good egg" and "bad egg" groups for plotting and analysis.
        '''
        precision = (self.param_combo_results.loc[j, 'TP'] + self.param_combo_results.loc[j, 'TP_noc5_obs'])/(self.param_combo_results.loc[j, 'TP'] + self.param_combo_results.loc[j, 'TP_noc5_obs'] +  self.param_combo_results.loc[j, 'FP_c5_noobs'] + self.param_combo_results.loc[j, 'FP_noc5_noobs'])
        recall = (self.param_combo_results.loc[j, 'TP'] + self.param_combo_results.loc[j, 'TP_noc5_obs'])/(self.param_combo_results.loc[j, 'TP'] + self.param_combo_results.loc[j, 'TP_noc5_obs'] + self.param_combo_results.loc[j, 'FN'] + self.param_combo_results.loc[j, 'FN_canc'])
        ttl_score = (self.param_combo_results.loc[j, 'TP'] + self.param_combo_results.loc[j, 'TP_noc5_obs'] + 
                        self.param_combo_results.loc[j, 'FP_c5_noobs'] + self.param_combo_results.loc[j, 'FP_noc5_noobs'])/(
                        self.param_combo_results.loc[j, 'TP'] + self.param_combo_results.loc[j, 'TP_noc5_obs'] +
                        self.param_combo_results.loc[j, 'FP_c5_noobs'] + self.param_combo_results.loc[j, 'FP_noc5_noobs'] +
                        self.param_combo_results.loc[j, 'TN_canc'] + self.param_combo_results.loc[j, 'FN_canc'])
        fbeta = (1.25*precision*recall)/((0.25*precision)+recall) 
        self.param_combo_results.loc[j, 'Precision'] = precision
        self.param_combo_results.loc[j, 'Recall'] = recall
        self.param_combo_results.loc[j, 'Trigger-to-Launch'] = ttl_score
        self.param_combo_results.loc[j, 'Param Combo'] = parameter
        
    def save_launch_DataFrame(self, parameter):
        self.launches_df.to_csv(f'{self.directory}/Launches/{parameter}_results.csv')
        print('launch dataframe saved!')

    def save_and_sort_comboDataFrames(self):
        ''' Sorting the parameter result csv files into a "good egg" and "bad egg" folder, depending on cutoffs for the
        recall and ttl scores. 
        '''   
        recall_min_score = 0.4
        ttl_min_score = 0.4
        good_egg = (self.param_combo_results['Recall']>0.4) & (self.param_combo_results['Trigger-to-Launch'] > 0.4)
        self.good_egg_combos = self.param_combo_results[good_egg]
        self.good_egg_combos.sort_values(by='Precision')
        self.bad_egg_combos = self.param_combo_results[~good_egg]
        self.bad_egg_combos.sort_values(by='Precision')
        
        self.good_egg_combos.to_csv(f'{self.directory}/GoodEggParams/Good_combo_results.csv')
        self.bad_egg_combos.to_csv(f'{self.directory}/BadEggParams/Bad_combo_results.csv')
        print('combo param df saved')
        
        
             
    


        
    


        
    
        
        
        