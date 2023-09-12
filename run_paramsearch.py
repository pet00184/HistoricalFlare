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
    pr.plotting_results(param_directory, launches_df_list, plot_directory_list)
    
if __name__ == '__main__':
    xrsb_value_search()