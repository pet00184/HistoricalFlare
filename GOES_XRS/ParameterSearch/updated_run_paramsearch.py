import updated_paramsearch as ps
import updated_save_scores as ss
import os
from astropy.io import fits
import numpy as np
import pandas as pd
import functools
import warnings
import sys
import multiprocessing as mp
warnings.filterwarnings("ignore")

calculated_params = '../GOES_computed_parameters.fits' #change depending on if you put your fits file somewhere else!
cparam_hdu = fits.open(calculated_params)
cparam = cparam_hdu[1].data

flare_fits = '../GOES_XRS_historical_finalversion.fits'
fitsfile = fits.open(flare_fits)
flare_data = fitsfile[1].data


################# Dictionary of all Parameters ################################################################        
params = {
    'xrsb': [[0, 1e-6, 2.5e-6, 5e-6, 7.5e-6], flare_data['xrsb'], 'W/m^2'],
    'xrsa': [[0, 2.5e-7, 3e-7, 3.5e-7, 4e-7, 4.5e-7, 5e-7], flare_data['xrsa'], 'W/m^2'],
    'bkinc': [[5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-6, 3e-6, 5e-6], cparam['Increase above Background'], 'W/m^2'],
    'bkincfrac': [[.1, .5, 1,2,3], cparam['Increase above Background Fraction'], ''],
    '1minxrsb': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 1-min Differences'], 'W/m^2'],
    '2minxrsb': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 2-min Differences'], 'W/m^2'],
    '3minxrsb': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 3-min Differences'], 'W/m^2'],
    '4minxrsb': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 4-min Differences'], 'W/m^2'],
    '5minxrsb': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 5-min Differences'], 'W/m^2'],
    '1minxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 1-min Differences'], 'W/m^2'],
    '2minxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 2-min Differences'], 'W/m^2'],
    '3minxrsa': [[0, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 3-min Differences'], 'W/m^2'],
    '4minxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 4-min Differences'], 'W/m^2'],
    '5minxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 5-min Differences'], 'W/m^2'],
    '1minpctxrsb': [[5, 10, 20, 30, 40, 50], cparam['XRSB 1-min Differences %'], '%'],
    '2minpctxrsb': [[5, 10, 20, 30, 40, 50], cparam['XRSB 2-min Differences %'], '%'],
    '3minpctxrsb': [[5, 10, 20, 30, 40, 50], cparam['XRSB 3-min Differences %'], '%'],
    '4minpctxrsb': [[5, 10, 20, 30, 40, 50], cparam['XRSB 4-min Differences %'], '%'],
    '5minpctxrsb': [[5, 10, 20, 30, 40, 50], cparam['XRSB 5-min Differences %'], '%'],
    '1minpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 1-min Differences %'], '%'],
    '2minpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 2-min Differences %'], '%'],
    '3minpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 3-min Differences %'], '%'],
    '4minpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 4-min Differences %'], '%'],
    '5minpctxrsa': [[5, 10, 20, 30, 40, 50], cparam['XRSB 5-min Differences %'], '%'],
    '1mintemp': [[.05, .1, .25, .5, 1], cparam['Temp 1-min Differences'], ''],
    '2mintemp': [[.05, .1, .25, .5, 1], cparam['Temp 2-min Differences'], ''],
    '3mintemp': [[.05, .1, .25, .5, 1], cparam['Temp 3-min Differences'], ''],
    '4mintemp': [[.05, .1, .25, .5, 1], cparam['Temp 4-min Differences'], ''],
    '5mintemp': [[.05, .1, .25, .5, 1], cparam['Temp 5-min Differences'], ''],
    'temp': [[0, .01, .05, .1, .15], cparam['Temperature (xrsa/xrsb)'], ''], 
    'em': [[0, 5e47, 1e48, 2e48, 3e48, 4e48, 5e48, 6e48, 7e48], flare_data['em'], 'cm^-3'],
    '1minem': [[1e47, 3e47, 5e47], flare_data['1-min em diff'], 'cm^-3'], 
    '2minem': [[1e47, 3e47, 5e47, 7e47], flare_data['2-min em diff'], 'cm^-3'], 
    '3minem': [[1e47, 3e47, 5e47, 7e47], flare_data['3-min em diff'], 'cm^-3'], 
    '4minem': [[1e47, 3e47, 5e47, 7e47, 1e48], flare_data['4-min em diff'], 'cm^-3'], 
    '5minem': [[0, 1e47, 5e47, 7e47, 1e48, 3e48], flare_data['5-min em diff'], 'cm^-3'], 
    
    }
    
def make_param_info(keys_list):
    param_combinations = np.array(np.meshgrid(*[params[key][0] for key in keys_list])).T.reshape(-1, len(keys_list))
    param_arrays = list(zip(*[params[key][1] for key in keys_list]))
    param_units = [params[key][2] for key in keys_list]
    return keys_list, param_combinations, param_arrays, param_units  

################################################################################################################
    
def run_paramsearch(param_directory, param_names, param_units, param_array_list, param_combo_list):
    #tag, param_combo_list = param_combo_list #keep this saved so I can remember it for the param combos stuff
    param_search = ps.ParameterSearch(param_names, param_units, param_array_list, param_combo_list, param_directory)
    param_search.loop_through_parameters()
    
def run_multiprocessing_paramsearch(keys_list, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    param_names, param_combinations, param_arrays, param_units = make_param_info(keys_list)
    #getting the number of cores for the slurm job
    try:
        num_cores = int(sys.argv[1])
    except IndexError:
        num_cores = os.cpu_count()
    print('num cores used:', num_cores)
    print('Total Params:', len(param_combinations))
    #splitting the array so all available cores are used
    splitup = np.array_split(param_combinations, num_cores) 
    #splitup = [[i, s] for (i, s) in enumerate(splitup)] #dont need this rn since I will be doing the other stuff later!
    #doing the multiple run!
    call_me = functools.partial(run_paramsearch, out_dir, param_names, param_units, param_arrays)
    with mp.Pool(num_cores) as p:
        p.map(call_me, splitup)
        
##################################################################################################################

def run_savescores(out_dir, param_names, param_units, launches_and_tag):
    tag, launch_df_list = launches_and_tag
    save_scores = ss.SaveScores(out_dir, launch_df_list, tag, param_names, param_units)
    save_scores.loop_through_param_combos() 
    
def run_multiprocessing_savescores(keys_list, out_dir):
    param_names, _, _, param_units = make_param_info(keys_list)
    launches_list = np.array(os.listdir(os.path.join(out_dir, 'Launches')))
    print(launches_list)
    try:
        num_cores = int(sys.argv[1])
    except IndexError:
        num_cores = os.cpu_count()
    print('num cores used:', num_cores)
    print('Number of combinations:', len(launches_list))  
    splitup = np.array_split(launches_list, num_cores)
    splitup = [[i, s] for (i, s) in enumerate(splitup)]
    call_me = functools.partial(run_savescores, out_dir, param_names, param_units)
    with mp.Pool(num_cores) as p:
        p.map(call_me, splitup)
    make_large_df(keys_list, out_dir)

def make_large_df(keys_list, out_dir):
    if os.path.exists(os.path.join(out_dir, 'AllParameterScores.csv')):
        os.remove(os.path.join(out_dir, 'AllParameterScores.csv'))
    score_files = [f for f in os.listdir(out_dir) if f.endswith('.csv')]
    total_score_df = pd.concat([pd.read_csv(os.path.join(out_dir, file), index_col=0) for file in score_files], ignore_index=True)
    total_score_df = total_score_df.sort_values(by=keys_list)
    total_score_df = total_score_df.reset_index(drop=True)
    total_score_df.to_csv(os.path.join(out_dir, 'AllParameterScores.csv'))
    print('All parameter scores saved.')
    for score_file in score_files:
        os.remove(os.path.join(out_dir, score_file))
    
##################################################################################################################  

if __name__ == '__main__':
    keys_list = ['xrsa', 'xrsb', '3minxrsa', 'temp', 'em']
    os.makedirs('RESULTS', exist_ok=True)
    savestring = "_".join(keys_list)
    out_dir = os.path.join('RESULTS', savestring)
    run_multiprocessing_paramsearch(keys_list, out_dir)
    run_multiprocessing_savescores(keys_list, out_dir)
