import parameter_search as ps
import new_plotting as pr
import os
from astropy.io import fits
import numpy as np
import multiprocessing as mp
import functools
import warnings
import sys
import pandas as pd
warnings.filterwarnings("ignore")

calculated_params = 'GOES_computed_parameters.fits' #change depending on if you put your fits file somewhere else!
cparam_hdu = fits.open(calculated_params)
cparam = cparam_hdu[1].data

flare_fits = 'GOES_XRS_historical_finalversion.fits'
fitsfile = fits.open(flare_fits)
flare_data = fitsfile[1].data


################# Dictionary of all Parameters ################################################################        
params = {
    'xrsb': [[0, 1e-6, 2.5e-6, 5e-6, 7.5e-6], flare_data['xrsb']],
    'xrsa': [[0, 2.5e-7, 3e-7, 3.5e-7, 4e-7, 4.5e-7, 5e-7], flare_data['xrsa']],
    'bkgrdinc': [[5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-6, 3e-6, 5e-6], cparam['Increase above Background']],
    'bkgrdincfrac': [[.1, .5, 1,2,3], cparam['Increase above Background Fraction']],
    '1mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 1-min Differences']],
    '2mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 2-min Differences']],
    '3mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 3-min Differences']],
    '4mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 4-min Differences']],
    '5mindiff': [[0, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 5-min Differences']],
    '1mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 1-min Differences']],
    '2mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 2-min Differences']],
    '3mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 3-min Differences']],
    '4mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 4-min Differences']],
    '5mindiffxrsa': [[0, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 5-min Differences']],
    '1mindiffpct': [[5, 10, 20, 30, 40, 50], cparam['XRSB 1-min Differences %']],
    '2mindiffpct': [[5, 10, 20, 30, 40, 50], cparam['XRSB 2-min Differences %']],
    '3mindiffpct': [[5, 10, 20, 30, 40, 50], cparam['XRSB 3-min Differences %']],
    '4mindiffpct': [[5, 10, 20, 30, 40, 50], cparam['XRSB 4-min Differences %']],
    '5mindiffpct': [[5, 10, 20, 30, 40, 50], cparam['XRSB 5-min Differences %']],
    '1mindiffpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 1-min Differences %']],
    '2mindiffpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 2-min Differences %']],
    '3mindiffpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 3-min Differences %']],
    '4mindiffpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 4-min Differences %']],
    '5mindiffpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 5-min Differences %']],
    '1temp_diff': [[.05, .1, .25, .5, 1], cparam['Temp 1-min Differences']],
    '2temp_diff': [[.05, .1, .25, .5, 1], cparam['Temp 2-min Differences']],
    '3temp_diff': [[.05, .1, .25, .5, 1], cparam['Temp 3-min Differences']],
    '4temp_diff': [[.05, .1, .25, .5, 1], cparam['Temp 4-min Differences']],
    '5temp_diff': [[.05, .1, .25, .5, 1], cparam['Temp 5-min Differences']],
    'temp': [[0, .01, .05, .1, .15], cparam['Temperature (xrsa/xrsb)']], 
    'em': [[0, 2e48, 3e48, 4e48, 5e48, 6e48, 7e48], flare_data['em']],
    '1em_diff': [[1e47, 3e47, 5e47], flare_data['1-min em diff']], 
    '2em_diff': [[1e47, 3e47, 5e47, 7e47], flare_data['2-min em diff']], 
    '3em_diff': [[1e47, 3e47, 5e47, 7e47], flare_data['3-min em diff']], 
    '4em_diff': [[1e47, 3e47, 5e47, 7e47, 1e48], flare_data['4-min em diff']], 
    '5em_diff': [[0, 1e47, 5e47, 7e47, 1e48, 3e48], flare_data['5-min em diff']], 
    
    }
################################################################################################################

    
def run_paramsearch(param_directory, param_array_list, param_combo_list):
    #tag, param_combo_list = param_combo_list
    
    run = Running_ParamSearch(param_directory, tag)
    run.run_paramsearch(param_combo_list, param_array_list)
    
def msi_slurm(param_grid, param_array_list, out_dir):
    #getting the number of cores for the slurm job
    try:
        num_cores = int(sys.argv[1])
    except IndexError:
        num_cores = os.cpu_count()
    print('num cores used:', num_cores)
    #making fancy param combinations list
    param_combo_list1 = np.array(param_grid).T.reshape(-1, len(param_grid))
    print('Total Params':, len(param_combo_list1))
    #splitting the array so all available cores are used
    splitup = np.array_split(param_combo_list1, num_cores) 
    #splitup = [[i, s] for (i, s) in enumerate(splitup)] #dont need this rn since I will be doing the other stuff later!
    #doing the multiple run!
    call_me = functools.partial(multi_run, out_dir, param_array_list)
    with mp.Pool(num_cores) as p:
        p.map(call_me, splitup)
        
def make_    
    
    
    
    
    
    

if __name__ == '__main__':
    
    
#trying a sort of big kahuna on MSI!!
param_grid = np.meshgrid(params['xrsb'][0], params['xrsa'][0], params['5mindiffxrsa'][0], params['temp'][0], params['em'][0])
param_array_list1 = list(zip(params['xrsb'][1], params['xrsa'][1], params['5mindiffxrsa'][1], params['temp'][1], params['em'][1]))
#making directories
os.makedirs('larger_test', exist_ok=True)
out_dir = os.path.join('larger_test', 'flux_5minfluxxrsa_temp_em')
msi_slurm(param_grid, param_array_list1, out_dir)
save_best_dfs(out_dir)
    


### graveyard:
# def save_best_dfs(out_dir):
#     good_combo_path = os.path.join(out_dir, 'GoodEggParams')
#     good_files = os.listdir(good_combo_path)
#     print(good_files)
#     good_df = pd.concat([pd.read_csv(os.path.join(good_combo_path, f)) for f in good_files ], ignore_index=True)
#     bad_combo_path = os.path.join(out_dir, 'BadEggParams')
#     bad_files = os.listdir(bad_combo_path)
#     bad_df = pd.concat([pd.read_csv(os.path.join(bad_combo_path, f)) for f in bad_files ], ignore_index=True)
#
#     hip_good = good_df[good_df['Precision'] > 0.85]
#     hip_bad = bad_df[bad_df['Precision'] > 0.85]
#     hip_df = pd.concat([hip_good, hip_bad], ignore_index=True)
#     hip_df.sort_values(by=['Precision'], ascending=False)
#     hip_df.to_csv(os.path.join(out_dir, 'HighPrecisionParams.csv'))
#
#     best_params = hip_good[hip_good['Recall'] > 0.4]
#     best_params.sort_values(by=['Precision'], ascending=False)
#     best_params.to_csv(os.path.join(out_dir, 'BestParams.csv'))
#
#     good_df.sort_values(by=['Precision'], ascending=False)
#     good_df.to_csv(os.path.join(good_combo_path, 'Good_combo_results.csv'))
#
#     bad_df.sort_values(by=['Precision'], ascending=False)
#     bad_df.to_csv(os.path.join(bad_combo_path, 'Bad_combo_results.csv'))


# class Running_ParamSearch:
#
#     def __init__(self, param_directory):
#         self.param_directory = param_directory
#         os.makedirs(self.param_directory, exist_ok=True)
#
#     def run_paramsearch(self, combos, combo_arrays):
#         param_search = ps.ParameterSearch(combos, self.param_directory, self.tag)
#         param_search.loop_through_parameters(combo_arrays)

    #################### trying to figure out a slurm job!!!!! ###################################
    #just emission measure and temp:
    # param_grid = np.meshgrid(params['temp'][0], params['em'][0])
    # param_array_list1 = list(zip(params['temp'][1], params['em'][1]))
    # #making directories
    # os.makedirs('slurmyslurmy', exist_ok=True)
    # out_dir = os.path.join('slurmyslurmy', 'temp_em_tester')
    # msi_slurm(param_grid, param_array_list1, out_dir)
    # save_best_dfs(out_dir)
    
