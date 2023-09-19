import parameter_search as ps
import plotting_results as pr
import os
from astropy.io import fits
import numpy as np
import math

calculated_params = 'GOES_computed_parameters.fits' #change depending on if you put your fits file somewhere else!
cparam_hdu = fits.open(calculated_params)
cparam = cparam_hdu[1].data


####### XRSB Flux Level Parameter Search ################
def xrsb_value_search():
    ''' Parameter Search for the XRSB flux level. Baseline search used mostly to get the parameter search and plotting up and 
    running.

    Values tried: [C1, C2.5, C5, C7.5]

    Results: 
    As expected, the results are not great! The C1/C2.5 triggers are too sensitive, and we launch on too many flares that aren't
    C5 or above. The C7.5 launch is very conservative, making it so most observations are successful, but we only launch on 
    ~10% of the flares.
    '''
    xrsb_level = [1e-6, 2.5e-6, 5e-6, 7.5e-6]  
    param_directory = 'XRSB_FluxValue'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)
    
    param_search = ps.ParameterSearch(xrsb_level, 'xrsb_fulltry', param_directory)
    param_search.loop_through_parameters(param_search.data['xrsb'])

    launches_df_list = ['1e-06_xrsb_fulltry_results.csv', '2.5e-06_xrsb_fulltry_results.csv','5e-06_xrsb_fulltry_results.csv',
                '7.5e-06_xrsb_fulltry_results.csv']
    plot_directory_list = ['1e-06_results', '2.5e-06_results', '5e-06_results', '7.5e-06_results']
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, xrsb_level)
    
def xrsa_value_search():
    ''' Parameter search for the XRSA flux level. This is also somewhat a baseline search, but hopefully it proves more helpful!
    
    '''
    xrsa_level = [1e-7, 1.5e-7, 2e-7, 2.5e-7, 3e-7, 3.5e-7, 4e-7, 4.5e-7, 5e-7, 5.5e-7]
    param_directory = 'XRSA_FluxValue'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)
        
    param_search = ps.ParameterSearch(xrsa_level, 'xrsa_values', param_directory)
    param_search.loop_through_parameters(param_search.data['xrsa'])

    launches_df_list = ['1e-07_xrsa_values_results.csv', '1.5e-07_xrsa_values_results.csv', 
                '2e-07_xrsa_values_results.csv', '2.5e-07_xrsa_values_results.csv', '3e-07_xrsa_values_results.csv', 
                '3.5e-07_xrsa_values_results.csv', '4e-07_xrsa_values_results.csv', '4.5e-07_xrsa_values_results.csv',
                '5e-07_xrsa_values_results.csv', '5.5e-07_xrsa_values_results.csv']
    plot_directory_list = ['1e-07_results', '1.5e-07_results', '2e-07_results', '2.5e-07_results', '3e-07_results',
                '3.5e-07_results', '4e-07_results', '4.5e-07_results', '5e-07_results', '5.5e-07_results',]
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, xrsa_level)
    
def increase_above_background_search(cparam):
    ''' Parameter search for the XRSB increase above background.
    '''    
    increase_level = [5e-7, 7e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-6, 3e-6, 5e-6]
    param_directory = 'IncreaseFromBackground'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)
        
    param_search = ps.ParameterSearch(increase_level, 'increase_above_background', param_directory)
    param_search.loop_through_parameters(cparam['Increase above Background'])
    
    launches_df_list = ['5e-07_increase_above_background_results.csv',
            '7e-07_increase_above_background_results.csv', '1e-06_increase_above_background_results.csv',
            '1.5e-06_increase_above_background_results.csv', '2e-06_increase_above_background_results.csv',
            '2.5e-06_increase_above_background_results.csv',
            '3e-06_increase_above_background_results.csv', '5e-06_increase_above_background_results.csv']
    plot_directory_list = ['5e-07_results', '7e-07_results', '1e-06_results', '1.5e-06_results', '2e-06_results',
             '2.5e-06_results','3e-06_results', '5e-06_results']
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, increase_level)
    
def increase_above_background_frac_search(cparam):
    ''' Parameter search for the XRSB increase above background, %.
    '''    
    increase_level = [1,2,3,4,5,10,12,15,20]
    param_directory = 'IncreaseFromBackgroundFraction'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)

    param_search = ps.ParameterSearch(increase_level, 'increase_above_background_fraction', param_directory)
    param_search.loop_through_parameters(cparam['Increase above Background Fraction'])

    launches_df_list = ['1_increase_above_background_fraction_results.csv', '2_increase_above_background_fraction_results.csv',
                '3_increase_above_background_fraction_results.csv', '4_increase_above_background_fraction_results.csv',
                '5_increase_above_background_fraction_results.csv',
                '10_increase_above_background_fraction_results.csv', '12_increase_above_background_fraction_results.csv',
                '15_increase_above_background_fraction_results.csv', '10_increase_above_background_fraction_results.csv',]
    plot_directory_list = ['1_results', '2_results', '3_results', '4_results','5_results', '10_results', '12_results', '15_results', '20_results',]
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, increase_level)
    
def min_increase_xrsb_search(cparam, n):
    ''' Parameter search for the 1-minute difference in flux between datapoints
    '''
    rough_diff = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6]
    param_directory = f'XRSBDifferences_{n}Min'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)

    param_search = ps.ParameterSearch(rough_diff, f'{n}_diff_XRSB', param_directory)
    param_search.loop_through_parameters(cparam[f'XRSB {n}-min Differences'])

    launches_df_list = [f'1e-08_{n}_diff_XRSB_results.csv', f'5e-08_{n}_diff_XRSB_results.csv', f'1e-07_{n}_diff_XRSB_results.csv',
                f'5e-07_{n}_diff_XRSB_results.csv', f'1e-06_{n}_diff_XRSB_results.csv', f'5e-06_{n}_diff_XRSB_results.csv',]
    plot_directory_list = ['1e-08_results', '5e-08_results', '1e-07_results', '5e-07_results', '1e-06_results', '5e-06_results',]
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, rough_diff)
    
def min_increase_xrsb_search_pct(cparam, n):
    ''' Parameter search for the 1-minute difference in flux between datapoints in percentage!
    '''
    rough_diff = [20, 30, 40, 50, 60, 70]
    param_directory = f'XRSBDifferencesPct_{n}Min'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)

    param_search = ps.ParameterSearch(rough_diff, f'{n}_diff_pct_XRSB', param_directory)
    param_search.loop_through_parameters(cparam[f'XRSB {n}-min Differences %'])

    launches_df_list = [f'20_{n}_diff_pct_XRSB_results.csv',f'30_{n}_diff_pct_XRSB_results.csv', f'40_{n}_diff_pct_XRSB_results.csv',
            f'50_{n}_diff_pct_XRSB_results.csv', f'60_{n}_diff_pct_XRSB_results.csv', f'70_{n}_diff_pct_XRSB_results.csv']
    plot_directory_list = ['20_results', '30_results', '40_results', '50_results', '60_results', '70_results']
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, rough_diff)
    
def temp_search(cparam):
    ''' Parameter search for the basic temperature xrsa/xrsb
    '''
    temp_levels = [.1, .15, .2, .25, .3]
    param_directory = 'Temperature'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)

    param_search = ps.ParameterSearch(temp_levels, f'temperature', param_directory)
    param_search.loop_through_parameters(cparam[f'Temperature (xrsa/xrsb)'])

    launches_df_list = ['0.1_temperature_results.csv', '0.15_temperature_results.csv', '0.2_temperature_results.csv',
            '0.25_temperature_results.csv', '0.3_temperature_results.csv',]
    plot_directory_list = ['0.1_results', '0.15_results', '0.2_results', '0.25_results', '0.3_results']
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, temp_levels)
    
def temp_difference_search(cparam, n):
    ''' Parameter search for difference temperatures (1-min difference xrsa/1-min difference xrsb), etc.
    '''
    rough_diff = [.5, 1, 5, 10]
    param_directory = f'TempDifferences_{n}Min'
    if not os.path.exists(param_directory):
        os.mkdir(param_directory)

    param_search = ps.ParameterSearch(rough_diff, f'{n}_diff_temp', param_directory)
    param_search.loop_through_parameters(cparam[f'Temp {n}-min Differences'])

    launches_df_list = [f'0.5_{n}_diff_temp_results.csv', f'1.0_{n}_diff_temp_results.csv', f'5.0_{n}_diff_temp_results.csv',
            f'10.0_{n}_diff_temp_results.csv',]
    plot_directory_list = ['0.5_results', '1_results', '5_results', '10_results']
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list, rough_diff)
    


if __name__ == '__main__':
   #temp_difference_search(cparam, 5)
