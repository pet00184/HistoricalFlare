import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import pandas as pd
import os
from astropy.io import fits


class PlottingResults:
    ''' New plotting class that will utilize the updated parameter combination saving method for nicer looking plots. 
    Things I want to achieve: 
        6. .txt file with nice summary data that I can just copy paste
        7. histograms of launches, cancellations, etc.... maybe make this its own class bc the data is different!!!
    '''
    
    def __init__(self, keys_list, nice_keys_list, out_dir, score_csv):
        self.score_df = pd.read_csv(os.path.join(out_dir, score_csv))
        self.out_dir = out_dir
        self.keys_list = keys_list
        self.nice_keys_list = nice_keys_list
        
        os.makedirs(os.path.join(self.out_dir, 'Plots'), exist_ok=True)
        
    def make_full_pr_plot(self, main_key):
        ''' Precision-Recall plot of all combinations tested during the run. 
        color_key = key used to color code the data, so that it is more easily readable. (In most cases, this will be xrsb)
        '''
        nice_main_key = self.nice_keys_list[np.where(np.array(self.keys_list) == main_key)[0][0]]
        main_key_values = self.score_df[main_key].unique()
        unit = self.score_df.loc[0, f'{main_key}_units']
        colors = cm.rainbow(np.linspace(0, 1, len(main_key_values)))
        fig, ax = plt.subplots()
        ax.set_xlabel(r'Recall = $\frac{\text{Launches Meeting Nominal Success Criteria}}{\text{All C5 Flares (All Potential Launches)}}$', fontsize=14)
        ax.set_ylabel(r'Precision = $\frac{\text{Launches Meeting Nominal Success Criteria}}{\text{All Launches}}$', fontsize=12)
        ax.set_title(f'Precision-Recall Curve \n Parameters: {", ".join(self.nice_keys_list)}', fontsize=13) 
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        for i, value in enumerate(main_key_values):
            main_key_df = self.score_df[self.score_df[main_key]==value].reset_index()
            ax.scatter(main_key_df['Recall'], main_key_df['Precision'], marker='o', c=[colors[i]]*main_key_df.shape[0], label = f'{nice_main_key}={value:.1e} {unit}')
        plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=7, title=f'Sorted by: {nice_main_key}')
        plt.savefig(os.path.join(self.out_dir, 'Plots', f'fullprplot_{main_key}coded_option4.png'), bbox_inches='tight', dpi=250)
        
    def make_optimal_pr_plot(self):
        optimal_df = pd.DataFrame()
        recall_bins = np.arange(0, 1.025, .025)
        for i in range(len(recall_bins)-1):
            recall_binned_df = self.score_df[(self.score_df['Recall'] > recall_bins[i]) & (self.score_df['Recall'] <= recall_bins[i+1])].reset_index(drop=True)
            recall_binned_df = recall_binned_df.sort_values(by='Precision', ascending=False).reset_index()
            if recall_binned_df.shape[0] > 0:
                optimal_df = optimal_df._append(recall_binned_df.iloc[0], ignore_index=True)
        fig, ax = plt.subplots()
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title(f'Optimal Precision-Recall Curve \n Parameters: {", ".join(self.nice_keys_list)}')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        for i, recall in enumerate(optimal_df['Recall']):
            label_list = []
            for j, key in enumerate(self.keys_list):
                key_units = f'{key}_units'
                label_list.append(f'{self.nice_keys_list[j]}={optimal_df.loc[i, key]:.1e} {optimal_df.loc[i, key_units]}')
            ax.plot(recall, optimal_df.loc[i, 'Precision'], marker='o', markersize=5, label=", ".join(label_list))
        plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=7, title=f'Optimal Param Combos:')
        plt.savefig(os.path.join(self.out_dir, 'Plots', 'fullprplot_optimal.png'), bbox_inches='tight', dpi=250)
        plt.close()
        
    def make_singlevarying_pr_plot(self, combo_list, main_key):
        ''' Makes a PR plot for a definited set of combination values, and all the options for a single value. For example,
        this could be used to see how the EM varies the results with fixed XRSA, XRSB and temp.
        imput: 
        combo_array = list (in same order as keys_list) of specific parameter combinations, aside from the one that 
        will be changing.
        main_key = main parameter that will vary. 
        '''
        nice_main_key = self.nice_keys_list[np.where(np.array(self.keys_list) == main_key)[0][0]]
        frozen_keys = self.keys_list
        nice_frozen_keys = self.nice_keys_list
        nice_frozen_keys.pop(np.where(np.array(frozen_keys) == main_key)[0][0])
        frozen_keys.remove(main_key)
        print(frozen_keys)
        print(nice_frozen_keys)
        frozen_key_units = [self.score_df.loc[0, f'{fkey}_units'] for fkey in frozen_keys]
        frozen_keyval_combos = [list(a) for a in zip(frozen_keys, nice_frozen_keys, combo_list, frozen_key_units)]
        frozen_df = self.score_df
        for frozen_combo in frozen_keyval_combos:
            frozen_df = frozen_df[frozen_df[frozen_combo[0]]==frozen_combo[2]].reset_index(drop=True)
        title_list = [f'{frozen_key[1]}={frozen_key[2]:.1e} {frozen_key[3]}' for frozen_key in frozen_keyval_combos]
        fig, ax = plt.subplots()
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title(f'Precision-Recall Curve Varying {nice_main_key}')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        for i, recall in enumerate(frozen_df['Recall']):
            main_units = frozen_df.loc[i, f'{main_key}_units']
            ax.plot(recall, frozen_df.loc[i, 'Precision'], marker='o', markersize=5, label=f'{nice_main_key}={frozen_df.loc[i, main_key]:.1e} {main_units}')
        plt.legend(loc='lower right', fontsize=8, title=f'{nice_main_key} values:')
        plt.text(1.01, 0.8, ("Fixed Parameters: \n" + "\n".join(title_list)), fontsize=10)
        fixed_str = [frozen[0] + str(frozen[2]) for frozen in frozen_keyval_combos]
        fixed_str = "_".join(fixed_str)
        plt.savefig(os.path.join(self.out_dir, 'Plots', f'prplot_varying{main_key}_{fixed_str}.png'), bbox_inches='tight', dpi=250)
        plt.close()
        
    def plot_specific_cf(self, param_dict, savestring, big_paramset=False):
        ''' Using a dictionary of parameters and keys (must be what was used for the chosen run), a CF matrix is plotted
        for that combination. 
        '''
        cf_df = self.score_df
        for key, val in param_dict.items():
            cf_df = cf_df[cf_df[key]==val].reset_index(drop=True)
        if cf_df.shape[0] > 1: raise ValueError('More than one combination left! Double check you have all parameters.')
        confusion_matrix = np.array([[cf_df.loc[0, 'TN'], cf_df.loc[0, 'TN_canc'], 
                        cf_df.loc[0, 'FP_noc5'], cf_df.loc[0, 'TP_noc5']],
                        [cf_df.loc[0, 'FN'], cf_df.loc[0, 'FN_canc'], 
                        cf_df.loc[0, 'FP_c5'], cf_df.loc[0, 'TP']]])
        fig, ax = plt.subplots()
        ax.matshow(confusion_matrix, cmap=plt.cm.Blues, alpha=0.7)
        for k in range(confusion_matrix.shape[0]):
            for j in range(confusion_matrix.shape[1]):
                ax.text(x=j, y=k,s=confusion_matrix[k, j], va='center', ha='center', size='xx-large')
        ylabels = ['', 'False', 'True']
        ax.set_ylabel('Flare above C5?', fontsize=14)
        ax.set_yticklabels(ylabels, rotation=90, fontsize=12)
        ax.xaxis.set_ticks_position("bottom")
        xlabels = ['', 'No Trigger', 'Cancelled Launch','Launch \n No C5 Obs', 'Launch \n C5 Obs']
        ax.set_xlabel('Launch Results', fontsize=14)
        ax.set_xticklabels(xlabels, fontsize=12)
        plt.title(f"Precision = {cf_df.loc[0, 'Precision']:.2f}, Recall = {cf_df.loc[0, 'Recall']:.2f}, \n Launch/Trigger Ratio = {cf_df.loc[0, 'LaunchTriggerRatio']:.2f}, Gordon Score = {cf_df.loc[0, 'Gordon']:.2f}")
        #lastly, need to make the stuff for showing the values!
        if not big_paramset:
            param_units = [cf_df.loc[0, f'{key}_units'] for key in param_dict.keys()]
            param_list = [f'{nice_key}={val} {unit}' for (key, val), unit, nice_key in zip(param_dict.items(), param_units, self.nice_keys_list)]
            plt.text(1.01, 0.8, ("Parameters \n" + "\n".join(param_list)), fontsize=12, transform=ax.transAxes)
        savestring = [f'{key}{val}' for key, val in param_dict.items()]
        savestring = "_".join(savestring)
        plt.savefig(os.path.join(self.out_dir, 'Plots', savestring, 'cf.png'), bbox_inches='tight', dpi=250)
        plt.close()
    
class LaunchPlotting:
    ''' Class for plotting specific launch/cancellation histograms! This will take in a dictionary, similar to how the cf
    is plotted. One challenge is opening the correct file! maybe doing the savestring naming convention...? It probably wouldn't
    work if the file name gets too long....
    '''

    def __init__(self, combo_dict, nice_keys_list, flare_fits, out_dir, savestring):
        self.combo_dict = combo_dict
        self.nice_keys_list = nice_keys_list
        self.launch_combo_list = os.listdir(os.path.join(out_dir, 'Launches'))
        fitsfile = fits.open(flare_fits)
        self.all_flare_data = fitsfile[1].data
        self.savestring = savestring

    def find_correct_launch_file(self):
        param_combo_string = "_".join([str(val) for val in self.combo_dict.values()])
        launch_csv_str = "_".join([param_combo_string, 'results.csv'])
        self.launch_combo_df = pd.read_csv(os.path.join(out_dir, 'Launches', launch_csv_str))
        
    def save_launch_cancellation_dfs(self):
        self.launch_df = self.launch_combo_df[self.launch_combo_df['Cancelled?']==False].reset_index()
        self.cancelled_df = self.launch_combo_df[self.launch_combo_df['Cancelled?']==True].reset_index()
        
    def plot_flare_histogram(self, cancellation=False):
        ''' plots histogram of all launches. Once that is working, I will try to add in something where a histogram of 
        all the flares is there in a low alpha in the background, so we can also see what we missed!!
        '''
        if cancellation:
            df = self.cancelled_df
        else:
            df = self.launch_df
        print(df)
        param_units = [df.loc[0, f'{key}_units'] for key in self.combo_dict.keys()]
        param_list = [f'{nice_key}={val} {unit}' for (key, val), unit, nice_key in zip(self.combo_dict.items(), param_units, self.nice_keys_list)]
        logbins=np.logspace(np.log10(1e-8),np.log10(1e-2), 50)
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        ax.hist(df['Flare_Max_Flux'], bins=logbins, range=(1e-8, 1e-2))
        ax.axvline(5e-6, c='r', lw=2, label='C5')
        ax.set_xscale('log')
        ax.set_xlabel('GOES Flux W/m$^2$', fontsize=14)
        ax.set_ylabel('Number of Flares', fontsize=14)
        ax.legend()
        plt.text(1.01, 0.8, ("Parameters: \n" + "\n".join(param_list)), fontsize=12, transform=ax.transAxes)
        if cancellation:
            ax.set_title(f'Maximum Flare Flux of Cancellations')
            plt.savefig(os.path.join(out_dir, 'Plots', 'maxflux_cancellations_histogram.png'), bbox_inches='tight', dpi=250)
        else:
            ax.set_title(f'Maximum Flare Flux of Launches', fontsize=16)
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxflux_launches_histogram_option1.png'), bbox_inches='tight', dpi=250)
            
    def plot_flare_histogram_includingallflares(self, cancellation=False):
        ''' Here is a histogram of all launches/cancellations, with all the potential flares plotted in the background.
        '''
        if cancellation:
            df = self.cancelled_df
        else:
            df = self.launch_df
        param_units = [df.loc[0, f'{key}_units'] for key in self.combo_dict.keys()]
        param_list = [f'{nice_key}={val} {unit}' for (key, val), unit, nice_key in zip(self.combo_dict.items(), param_units, self.nice_keys_list)]
        logbins=np.logspace(np.log10(1e-8),np.log10(5e-3), 50)
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        ax.hist(self.all_flare_data['peak flux'], bins=logbins, range=(1e-8, 1e-2), color='r', alpha=0.2, label='All Flares')
        if cancellation:
            ax.hist(df['Flare_Max_Flux'], bins=logbins, range=(1e-8, 1e-2), color='b', label='Cancelled Flares')
        else:
            ax.hist(df['Flare_Max_Flux'], bins=logbins, range=(1e-8, 1e-2), color='b', label='Launched Flares')
        ax.axvline(5e-6, c='k', lw=2, label='C5')
        ax.set_xscale('log')
        ax.set_xlabel('GOES Flux W/m$^2$')
        ax.set_ylabel('# of Flares')
        ax.legend()
        plt.text(1.01, 0.8, ("Parameters: \n" + "\n".join(param_list)), fontsize=12, transform=ax.transAxes)
        ax.set_title(f'Maximum Flare Flux')
        if cancellation:
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxflux_cancellations_allflaresshown_histogram.png'), bbox_inches='tight', dpi=250)
        else:
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxflux_launches_allflaresshown_histogram.png'), bbox_inches='tight', dpi=250)
        
    def plot_flare_histogram_includingallflares_candl(self):
        ''' Here is a histogram of all launches/cancellations, with all the potential flares plotted in the background.
        '''
        param_units = [self.launch_df.loc[0, f'{key}_units'] for key in self.combo_dict.keys()]
        param_list = [f'{nice_key}={val} {unit}' for (key, val), unit, nice_key in zip(self.combo_dict.items(), param_units, self.nice_keys_list)]
        logbins=np.logspace(np.log10(1e-8),np.log10(5e-3), 50)
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        ax.hist(self.all_flare_data['peak flux'], bins=logbins, range=(1e-8, 1e-2), color='r', alpha=0.2, label='All Flares')
        ax.hist(self.cancelled_df['Flare_Max_Flux'], bins=logbins, range=(1e-8, 1e-2), color='b', label='Cancelled Flares', alpha=0.5)
        ax.hist(self.launch_df['Flare_Max_Flux'], bins=logbins, range=(1e-8, 1e-2), color='k', label='Launched Flares', alpha=0.7)
        ax.axvline(5e-6, c='k', lw=1, ls='--', label='C5')
        ax.set_xscale('log')
        ax.set_xlabel('GOES Flux W/m$^2$')
        ax.set_ylabel('# of Flares')
        ax.legend()
        plt.text(1.01, 0.8, ("Parameters: \n" + "\n".join(param_list)), fontsize=12, transform=ax.transAxes)
        ax.set_title(f'Maximum Flare Flux')
        plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxflux_candl_allflaresshown_histogram.png'), bbox_inches='tight', dpi=250)
        
    def plot_observation_histograms(self, hic=False, cancellation=False):
        if cancellation:
            df = self.cancelled_df
        else:
            df = self.launch_df
        param_units = [df.loc[0, f'{key}_units'] for key in self.combo_dict.keys()]
        param_list = [f'{nice_key}={val} {unit}' for (key, val), unit, nice_key in zip(self.combo_dict.items(), param_units, self.nice_keys_list)]
        logbins=np.logspace(np.log10(1e-6),np.log10(5e-3), 40)
        fig, ax = plt.subplots(1,1,figsize=(8,6))
        ax.axvline(5e-6, c='k', lw=2, label='C5')
        ax.legend()
        ax.set_xscale('log')
        ax.set_xlabel('GOES Flux W/m$^2$', fontsize=14)
        ax.set_ylabel('Number of Flares', fontsize=14)
        #ax.set_xlim(1e-6, 2e-3)
        plt.text(1.01, 0.8, ("Parameters: \n" + "\n".join(param_list)), fontsize=12, transform=ax.transAxes)
        if hic:
            ax.hist(df['Max_HiC'], bins=logbins, range=(1e-8, 1e-2))
        else:
            ax.hist(df['Max_FOXSI'], bins=logbins, range=(1e-8, 1e-2))
        if hic and cancellation:
            ax.set_title(f'Maximum Flux Missed by HiC due to Cancellation')
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxhic_cancellation_histogram.png'), bbox_inches='tight', dpi=250)
        if not hic and cancellation:
            ax.set_title(f'Maximum Flux Missed by FOXSI due to Cancellation')
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxfoxsi_cancellation_histogram.png'), bbox_inches='tight', dpi=250)
        if not hic and not cancellation:
            ax.set_title(f'Maximum Observed Flux for FOXSI', fontsize=16)
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxfoxsi_launch_histogram_smallerrange3.png'), bbox_inches='tight', dpi=250)
        if hic and not cancellation:
            ax.set_title(f'Maximum Observed Flux for HiC', fontsize=16)
            plt.savefig(os.path.join(out_dir, 'Plots', self.savestring, 'maxhic_launch_histogram.png'), bbox_inches='tight', dpi=250)
    
        
if __name__ == '__main__':
    keys_list = ['xrsa', 'xrsb', '3minxrsa', 'temp', 'em']
    nice_keys_list = ['XRSA', 'XRSB', 'XRSA 3-min Increase', 'Temperature', 'Emission Measure']
    out_dir = os.path.join('RESULTS', 'xrsa_xrsb_3minxrsa_temp_em')
    score_csv = 'AllParameterScores.csv'
    test = PlottingResults(keys_list, nice_keys_list, out_dir, score_csv)
    # test.make_full_pr_plot('xrsb')
    # test.make_full_pr_plot('em')
    # test.make_optimal_pr_plot()
    # combo_list = [4.5e-7, 2.5e-6, 1e-8, .01]
    # test.make_singlevarying_pr_plot(combo_list, 'em')
    cf_dict = {
        'xrsa': 5e-7, #if we want a 0, we might need to check out how that saves in launch files...
        'xrsb': 0.0,
        '3minxrsa': 1e-8,
        'temp': 0.05,
        'em': 3e48
    } #needs to be in the same order as everything!!!
    savestring = [f'{key}{val}' for key, val in cf_dict.items()]
    savestring = '_'.join(savestring)
    os.makedirs(os.path.join(out_dir, 'Plots', savestring), exist_ok=True)
    test.plot_specific_cf(cf_dict, savestring)

    flare_fits = '../GOES_XRS_historical_finalversion.fits'

    ltest = LaunchPlotting(cf_dict, nice_keys_list, flare_fits, out_dir, savestring)
    ltest.find_correct_launch_file()
    ltest.save_launch_cancellation_dfs()
    # ltest.plot_flare_histogram()
    # ltest.plot_flare_histogram(cancellation=True)
    # ltest.plot_flare_histogram_includingallflares_candl()
    ltest.plot_flare_histogram_includingallflares()
    # ltest.plot_flare_histogram_includingallflares(cancellation=True)
    ltest.plot_observation_histograms()
    # ltest.plot_observation_histograms(hic=True)
    # ltest.plot_observation_histograms(cancellation=True)
    # ltest.plot_observation_histograms(hic=True, cancellation=True)
    
        
        