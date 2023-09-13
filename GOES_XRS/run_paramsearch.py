import parameter_search as ps
import plotting_results as pr
import os

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
    
if __name__ == '__main__':
    xrsa_value_search()