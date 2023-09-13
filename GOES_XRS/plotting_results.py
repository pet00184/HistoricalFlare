import sklearn.metrics as m
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import os

param_directory = 'XRSB_FluxValue'
launches_df_list = ['1e-06_xrsb_fulltry_results.csv', '2.5e-06_xrsb_fulltry_results.csv','5e-06_xrsb_fulltry_results.csv',
             '7.5e-06_xrsb_fulltry_results.csv']
plot_directory_list = ['1e-06_results', '2.5e-06_results', '5e-06_results', '7.5e-06_results']
    
class VisualizingData:
    
    def __init__(self, launches_df, save_directory):
        self.launches_df = pd.read_csv(launches_df)
        self.flare_fits = fits.open('GOES_XRS_historical.fits')
        self.flare_data = self.flare_fits[1].data
        self.flare_hdr = self.flare_fits[1].header
        self.total_launches = len(self.launches_df['Flare_ID'])
        self.total_flares = len(self.flare_data['flare ID'])
        self.save_directory = save_directory
        
    def plot_histogram_flux(self, flux_to_plot, plot_title):
        ''' Plots a histogram in logspace of the peak flux observed.
        parameters:
        ----------------------
        flux_to_plot (str) = column name of which flux to plot. (Either Max/Mean FOXSI/HiC or Flare Max Flux)
        '''
        flux = self.launches_df[flux_to_plot]
        logbins=np.logspace(np.log10(1e-8),np.log10(1e-2), 50)
        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        ax.hist(flux, bins=logbins, range=(1e-8, 1e-2))
        ax.axvline(5e-6, c='r', lw=2, label='C5 Flux')
        ax.axvline(1e-5, c='k', lw=2, label='M1 Flux')
        ax.set_xscale('log')
        ax.legend()
        ax.set_title(f'{plot_title} \n Total Launches = {self.total_launches}')
        plt.savefig(f'{self.save_directory}/{flux_to_plot}_histogram.png', dpi=200)
        plt.close()
        
    def plot_observed_confusionmatrix(self, bool_to_plot, plot_title):
        flare_values = self.launches_df['Flare_C5']
        observed_values = self.launches_df[bool_to_plot]
        conf_matrix = m.confusion_matrix(flare_values, observed_values)
        tn, fp, fn, tp = conf_matrix.ravel()
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1_score = 2 * (precision * recall) / (precision + recall)
        
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        ax.matshow(conf_matrix, cmap=plt.cm.Blues, alpha=0.7)
        for i in range(conf_matrix.shape[0]):
            for j in range(conf_matrix.shape[1]):
                ax.text(x=j, y=i,s=conf_matrix[i, j], va='center', ha='center', size='xx-large')
        plt.xlabel('Observation Above C5', fontsize=18)
        plt.ylabel('Flare Above C5', fontsize=18)
        plt.suptitle(f'Observation Confusion Matrix for {plot_title}', fontsize=18)
        plt.title(f'accuracy={accuracy:.2f}, precision={precision:.2f}, recall={recall:.2f}, F1={f1_score:.2f}', fontsize=12)
        plt.savefig(f'{self.save_directory}/{bool_to_plot}_obsconfmatrix.png', dpi=200)
        plt.close()
        
    def plot_launches_confusionmatrix(self):
        all_abovec5 = self.flare_data['above C5'] #this is the T/F array as the truth for the confusion matrix!
        all_flareID = self.flare_data['flare ID']
        flareID_launches = np.array(self.launches_df['Flare_ID'])
        launch_indx = []
        for flare in flareID_launches:
            indx = np.where(flare==all_flareID)[0][0]
            launch_indx.append(indx)
        
        pred_array = np.array([False]*len(all_flareID))
        np.put(pred_array, launch_indx, True)
        conf_matrix = m.confusion_matrix(all_abovec5, pred_array)
        tn, fp, fn, tp = conf_matrix.ravel()
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        precision = tp / (tp + fp)
        self.recall = tp / (tp + fn)
        f1_score = 2 * (precision * self.recall) / (precision + self.recall) 
        self.fpr = fp / (fp + tn)

        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        ax.matshow(conf_matrix, cmap=plt.cm.Blues, alpha=0.7)
        for i in range(conf_matrix.shape[0]):
            for j in range(conf_matrix.shape[1]):
                ax.text(x=j, y=i,s=conf_matrix[i, j], va='center', ha='center', size='xx-large')
        plt.xlabel('Launched Y/N', fontsize=18)
        plt.ylabel('Flare Above C5', fontsize=18)
        plt.suptitle(f'Confusion Matrix for Launches', fontsize=18)
        plt.title(f'accuracy={accuracy:.2f}, precision={precision:.2f}, recall={self.recall:.2f}, F1={f1_score:.2f}', fontsize=12)
        plt.savefig(f'{self.save_directory}/ynlaunchconfmatrix.png', dpi=200)
        plt.close()
       
    def make_textsummary_file(self, param_type, param_value):
        self.launches_aboveC5 = len(self.launches_df['Flare_C5'][self.launches_df['Flare_C5']==True])
        self.foxsimax_aboveC5 = len(self.launches_df['Flare_C5'][self.launches_df['Max_FOXSI_C5']==True])
        self.foxsimean_aboveC5 = len(self.launches_df['Flare_C5'][self.launches_df['Mean_FOXSI_C5']==True])
        self.hicmax_aboveC5 = len(self.launches_df['Flare_C5'][self.launches_df['Max_HiC_C5']==True])
        self.hicmean_aboveC5 = len(self.launches_df['Flare_C5'][self.launches_df['Mean_HiC_C5']==True])
        with open(f'{self.save_directory}/summary_file.txt', 'w') as f:
            f.write(f'Parameter Type: {param_type} \n')
            f.write(f'Parameter Value: {param_value} \n \n')
            f.write(f'Total Launches: {self.total_launches} of {self.total_flares}  Flares ({(self.total_launches/self.total_flares)*100:.1f}%)\n')
            f.write(f'Of the {self.total_launches} launches: \n')
            f.write(f'\t Flares C5 or above: {self.launches_aboveC5} ({self.launches_aboveC5/self.total_launches*100:.1f}%) \n')
            f.write(f'\t Max observed FOXSI flux C5 or above: {self.foxsimax_aboveC5} ({self.foxsimax_aboveC5/self.total_launches*100:.1f}%) \n')
            f.write(f'\t Mean observed FOXSI flux C5 or above: {self.foxsimean_aboveC5} ({self.foxsimean_aboveC5/self.total_launches*100:.1f}%) \n')
            f.write(f'\t Max observed HiC flux C5 or above: {self.hicmax_aboveC5} ({self.hicmax_aboveC5/self.total_launches*100:.1f}%) \n')
            f.write(f'\t Mean observed HiC flux C5 or above: {self.hicmean_aboveC5} ({self.hicmean_aboveC5/self.total_launches*100:.1f}%) \n')
      
def plotting_results(param_directory, launches_df_list, plot_directory_list, param_level):
    roc_list = []
    for i, launches_df in enumerate(launches_df_list):
        #making launch and save directories, and initializing class
        launch_directory = os.path.join(param_directory, launches_df)
        save_directory = os.path.join(param_directory, plot_directory_list[i])
        if not os.path.exists(save_directory):
            os.mkdir(save_directory)
        plottest = VisualizingData(launch_directory, save_directory)
        #doing histogram plots for all observations
        hist_loop = [['Flare_Max_Flux', 'Peak Flare Flux'], ['Max_FOXSI', 'Max Flux Observed by FOXSI'], ['Mean_FOXSI', 'Mean Flux Observed by FOXSI'],
                ['Max_HiC', 'Max Flux Observed by HiC'], ['Mean_HiC', 'Mean Flux Observed by HiC']]
        for group in hist_loop:
            plottest.plot_histogram_flux(group[0], group[1])
        #doing confusion matrices for observations
        obsconf_loop = [['Max_FOXSI_C5', 'Max FOXSI Flux'], ['Mean_FOXSI_C5', 'Mean FOXSI Flux'], ['Max_HiC_C5', 'Max HiC Flux'], ['Mean_HiC_C5', 'Mean HiC Flux']]
        for group in obsconf_loop:
            plottest.plot_observed_confusionmatrix(group[0], group[1])
        #doing confusion matrices for y/n launch 
        plottest.plot_launches_confusionmatrix()
        roc_list.append([plottest.recall, plottest.fpr, param_level[i]])
        #saving summary text file
        plottest.make_textsummary_file(param_directory, plot_directory_list[i])
        print(f'all plots and summary file complete for {launches_df}!')
    roc_plot(roc_list, param_directory)
    print('roc plot done!')
    
def roc_plot(roc_list, param_directory):
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(f'ROC Space for {param_directory}')
    ax.plot([0, 1], [0, 1], ls='--', c='r', lw=2, label='Random Classifier')
    for param in roc_list:
        ax.plot(param[1], param[0], marker='o', markersize=5, label=str(param[2]))
    ax.legend()
    plt.savefig(f'{param_directory}/ROC_plot.png', dpi=200)
            
# if __name__ == "__main__":
#
#     for i, launches_df in enumerate(launches_df_list):
#         #making launch and save directories, and initializing class
#         launch_directory = os.path.join(param_directory, launches_df)
#         save_directory = os.path.join(param_directory, plot_directory_list[i])
#         plottest = VisualizingData(launch_directory, save_directory)
#         #doing histogram plots for all observations
#         hist_loop = [['Flare_Max_Flux', 'Peak Flare Flux'], ['Max_FOXSI', 'Max Flux Observed by FOXSI'], ['Mean_FOXSI', 'Mean Flux Observed by FOXSI'],
#                 ['Max_HiC', 'Max Flux Observed by HiC'], ['Mean_HiC', 'Mean Flux Observed by HiC']]
#         for group in hist_loop:
#             plottest.plot_histogram_flux(group[0], group[1])
#         #doing confusion matrices for observations
#         obsconf_loop = [['Max_FOXSI_C5', 'Max FOXSI Flux'], ['Mean_FOXSI_C5', 'Mean FOXSI Flux'], ['Max_HiC_C5', 'Max HiC Flux'], ['Mean_HiC_C5', 'Mean HiC Flux']]
#         for group in obsconf_loop:
#             plottest.plot_observed_confusionmatrix(group[0], group[1])
#         #doing confusion matrices for y/n launch
#         plottest.plot_launches_confusionmatrix()
#         #saving summary text file
#         plottest.make_textsummary_file(param_directory, plot_directory_list[i])
        

            
            
            
    