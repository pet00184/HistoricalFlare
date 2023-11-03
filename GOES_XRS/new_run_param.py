import parameter_search as ps
import new_plotting as pr
import os
from astropy.io import fits
import numpy as np
import math
import warnings
warnings.filterwarnings("ignore")

calculated_params = 'GOES_computed_parameters.fits' #change depending on if you put your fits file somewhere else!
cparam_hdu = fits.open(calculated_params)
cparam = cparam_hdu[1].data

flare_fits = 'GOES_XRS_historical_finalversion.fits'
fitsfile = fits.open(flare_fits)
flare_data = fitsfile[1].data


################# Dictionary of all Parameters ################################################################        
params = {
    'xrsb': [[1e-6, 2.5e-6, 5e-6, 7.5e-6], flare_data['xrsb']],
    'xrsa': [[2.5e-7, 3e-7, 3.5e-7, 4e-7, 4.5e-7, 5e-7], flare_data['xrsa']],
    'bkgrdinc': [[5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-6, 3e-6, 5e-6], cparam['Increase above Background']],
    'bkgrdincfrac': [[.1, .5, 1,2,3], cparam['Increase above Background Fraction']],
    '1mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 1-min Differences']],
    '2mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 2-min Differences']],
    '3mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 3-min Differences']],
    '4mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 4-min Differences']],
    '5mindiff': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSB 5-min Differences']],
    '1mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 1-min Differences']],
    '2mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 2-min Differences']],
    '3mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 3-min Differences']],
    '4mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 4-min Differences']],
    '5mindiffxrsa': [[1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6], cparam['XRSA 5-min Differences']],
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
    '1temp_diff': [[.05, .1, .25, .5, 1], cparam[f'Temp 1-min Differences']],
    '2temp_diff': [[.05, .1, .25, .5, 1], cparam[f'Temp 2-min Differences']],
    '3temp_diff': [[.05, .1, .25, .5, 1], cparam[f'Temp 3-min Differences']],
    '4temp_diff': [[.05, .1, .25, .5, 1], cparam[f'Temp 4-min Differences']],
    '5temp_diff': [[.05, .1, .25, .5, 1], cparam[f'Temp 5-min Differences']],
    'temp': [[.01, .05, .1, .15], cparam[f'Temperature (xrsa/xrsb)']], 
    'em': [[2e48, 3e48, 4e48, 5e48, 6e48, 7e48], flare_data['em']],
    '1em_diff': [[1e47, 3e47, 5e47], flare_data['1-min em diff']], 
    '2em_diff': [[1e47, 3e47, 5e47, 7e47], flare_data['2-min em diff']], 
    '3em_diff': [[1e47, 3e47, 5e47, 7e47], flare_data['3-min em diff']], 
    '4em_diff': [[1e47, 3e47, 5e47, 7e47, 1e48], flare_data['4-min em diff']], 
    '5em_diff': [[1e47, 5e47, 7e47, 1e48, 3e48], flare_data['5-min em diff']], 
    
    }
################################################################################################################


class Running_ParamSearch:
    
    def __init__(self, param_directory, multi=False):
        self.multi = multi
        self.param_directory = param_directory
        if not os.path.exists(self.param_directory):
            os.mkdir(self.param_directory)
            
    def run_paramsearch(self, combos, combo_arrays):
        param_search = ps.ParameterSearch(combos, self.param_directory, self.multi)
        param_search.loop_through_parameters(combo_arrays)
    
        plotting = pr.Assessing_Data(f'{self.param_directory}/GoodEggParams/Good_combo_results.csv', 
                        f'{self.param_directory}/BadEggParams/Bad_combo_results.csv', self.param_directory)
        plotting.make_all_plots()
    
def single_run(param_directory, param, array):
    run = Running_ParamSearch(param_directory)
    run.run_paramsearch(param, array)
    
def multi_run(param_directory, param_combo_list, param_array_list):
    combo_params = np.array(param_combo_list).T.reshape(-1, len(param_combo_list))
    combo_arrays = param_array_list
    run = Running_ParamSearch(param_directory, multi=False)
    run.run_paramsearch(combo_params, combo_arrays)
    
    
    
if __name__ == '__main__':
    # #going through each value independently
#     for key, value in params.items():
#         print(f'DOING {key} RUN!!')
#         single_run(os.path.join('Single_msi_test', key), value[0], value[1])
    
    
    #single_run(os.path.join('New_Single', 'em'), params['em'][0], params['em'][1])
    # single_run(os.path.join('New_Single', '1em_diff'), params['1em_diff'][0], params['1em_diff'][1])
    # single_run(os.path.join('New_Single', '2em_diff'), params['2em_diff'][0], params['2em_diff'][1])
    # single_run(os.path.join('New_Single', '3em_diff'), params['3em_diff'][0], params['3em_diff'][1])
    # single_run(os.path.join('New_Single', '4em_diff'), params['4em_diff'][0], params['4em_diff'][1])
    # single_run(os.path.join('New_Single', '5em_diff'), params['5em_diff'][0], params['5em_diff'][1])
        
    #xrsb and others param search:
    # for key, value in params.items():
    #     if key=='xrsb':
    #         continue
    #     print(f'DOING XRSB + {key} RUN!!!')
    #     param_combo_list = np.meshgrid(params['xrsb'][0], value[0],)
    #     param_array_list = list(zip(params['xrsb'][1], value[1]))
    #     multi_run(os.path.join('TwoParam_V4', key), param_combo_list, param_array_list)
        
    # # #triple param test
    # param_combo_list = np.meshgrid(params['xrsb'][0], params['xrsa'][0], params['temp'][0])
    # param_array_list = list(zip(params['xrsb'][1], params['xrsa'][1], params['temp'][1]))
    # multi_run(os.path.join('ThreeParam_V4', 'xrsb_xrsa_temp'), param_combo_list, param_array_list)
    
    # #quad param test:
    # param_combo_list = np.meshgrid(params['xrsb'][0], params['xrsa'][0], params['temp'][0], params['em'][0])
    # param_array_list = list(zip(params['xrsb'][1], params['xrsa'][1], params['temp'][1], params['em'][1]))
    # multi_run(os.path.join('FourPlusParam_V4', 'xrsb_xrsa_temp_em'), param_combo_list, param_array_list)
    
    #just emission measure and temp:
    param_combo_list1 = np.meshgrid(params['temp'][0], params['em'][0])
    param_array_list1 = list(zip(params['temp'][1], params['em'][1]))
    multi_run(os.path.join('Tests_1031', 'temp_em'), param_combo_list1, param_array_list1)
    