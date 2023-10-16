import sklearn.metrics as m
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from astropy.io import fits
import os
import itertools

class Assessing_Data:
    ''' The goal of this new plotting class is to take into account the scoring of what we define as success (more than before).
    I am going to make a few edited data products that will be saved for each run. I also hope to automate a way of 
    ranking the different parameter searches.
    
    Success: 
    We are rating success as a launch where both HiC and FOXSI observe C5 flux (both maximums are C5 or above). For missed 
    flares, anything that we didn't launch on but was C5 counts as something we would want to launch on (FN). 
    
    Data products made: 
    1. Traditional y/n launch confusion matrix. This is basically the same as what I had before! It might be helpful
    to save, but it is mostly an initial step for the new confusion matrix.
    2. 3x2 funky confusion matrix. This takes into account both launching and observing. By doing this, I can edit the FP and 
    FPR in order to indclude launches where the flare was C5 but we missed C5 flux. By doing so, this will edit the 
    FPR and TPR which will hopefully give us more accurate ROC plots.
    3. The edited ROC plots
    4. A dictionary that includes scores for each parameter combination. This can then be ranked!
    5. Still do a little text summary bc it is helpful
    '''
    
    def __init__(self, launches_df, param_directory, param_level):
        self.launches_df = pd.read_csv(launches_df)
        self.flare_fits = fits.open('GOES_XRS_historical.fits')
        self.flare_data = self.flare_fits[1].data
        self.flare_hdr = self.flare_fits[1].header
        self.total_launches = len(self.launches_df['Flare_ID'])
        self.total_flares = len(self.flare_data['flare ID'])
        self.param_level = param_level
        self.param_directory = param_directory
        
    def save_launch_and_obs_arrays(self):
        '''This saves arrays of true/false for the baseline truth (was the flare above C5), the launch prediction
        (was there a launch) and the launch/obs prediction (launch and obs). From this, we can sort each box!
        '''
        #defining the array that determines if a flare is above C5
        self.all_abovec5 = self.flare_data['above C5'] #all flares above C5- this is the baseline "truth"
        
        #make an array that determines whether there was a launch
        all_flareID = self.flare_data['flare ID']
        flareID_launches = np.array(self.launches_df['Flare_ID'])
        launch_indx = []
        for flare in flareID_launches:
            indx = np.where(flare==all_flareID)[0][0]
            launch_indx.append(indx)
        self.launch_array = np.array([False]*len(all_flareID))
        np.put(self.launch_array, launch_indx, True)
        
        #make an array that defines where flux above C5 was observed by both FOXSI and HIC
        flareID_observations = np.array(self.launches_df['Flare_ID'][self.launches_df['Max_FOXSI_and_HiC_C5']==True])
        obs_indx = []
        for flare in flareID_observations:
            indx = np.where(flare==all_flareID)[0][0]
            obs_indx.append(indx)
        self.obs_array = np.array([False]*len(all_flareID))
        np.put(self.obs_array, obs_indx, True)
        
    def make_confusion_matrix(self):
        ''' This is going to be a lot more hands on since we are not longer working with nxn classification and can just
        use a pre-built confusion matrix function. You should end up with a 3x2 array that has the right value for each
        one in it, to be used for plotting.
        '''
        self.TN = len(np.where((self.all_abovec5==False) & (self.launch_array==False))[0])
        self.TP = len(np.where((self.all_abovec5==True) & (self.obs_array==True))[0])
        self.FN = len(np.where((self.all_abovec5==True) & (self.launch_array==False))[0])
        self.FP_c5_noobs = len(np.where((self.all_abovec5==True) & (self.launch_array==True) & (self.obs_array==False))[0])
        self.FP_noc5_noobs = len(np.where((self.all_abovec5==False) & (self.launch_array==True) & (self.obs_array==False))[0])
        self.noc5_obs = len(np.where((self.all_abovec5==False) & (self.launch_array==True) & (self.obs_array==True))[0])
        self.confusion_matrix = np.array([[self.TN, self.FP_noc5_noobs, self.noc5_obs], [self.FN, self.FP_c5_noobs, self.TP]])
        
    def save_parameter_scores(self):
        ''' Parameters we care about is mostly the precision, since our hope is to minimize the FP rate (times we launched
        but didn't observe C5 flux!). I'm willing to up the FN rate if that increases the precision.
        
        Basically, we want to maximize precision and minimize FPR
        '''
        self.precision = (self.TP + self.noc5_obs) / (self.FP_c5_noobs + self.FP_noc5_noobs + self.TP + self.noc5_obs)
        self.recall = (self.TP + self.noc5_obs) / (self.TP + self.noc5_obs + self.FN)
        self.fpr = (self.FP_c5_noobs + self.FP_noc5_noobs) / (self.FP_c5_noobs + self.FP_noc5_noobs + self.TN)
        self.fbeta_measure = (1.25*self.precision*self.recall)/((0.25*self.precision)+self.recall) #puts more weight on minimizing false positives!!
        
    def plot_confusion_matrix(self):
        ''' This just plots the new and improved confusion matrix for the specific parameter!
        '''
        fig, ax = plt.subplots()
        ax.matshow(self.confusion_matrix, cmap=plt.cm.Blues, alpha=0.7)
        for i in range(self.confusion_matrix.shape[0]):
            for j in range(self.confusion_matrix.shape[1]):
                ax.text(x=j, y=i,s=self.confusion_matrix[i, j], va='center', ha='center', size='xx-large')
        ylabels = ['', 'False', 'True']
        ax.set_ylabel('Flare above C5?', fontsize=14)
        ax.set_yticklabels(ylabels, rotation=90, fontsize=12)
        ax.xaxis.set_ticks_position("bottom")
        xlabels = ['', 'No Launch', 'Launch \n No C5 Obs', 'Launch \n C5 Obs']
        ax.set_xlabel('Launch Results', fontsize=14)
        ax.set_xticklabels(xlabels, fontsize=12)
        plt.suptitle(f'Parameters: {self.param_level}')
        plt.title(f'Precision = {self.precision:.2f}, FBeta = {self.fbeta_measure:.2f}')
        plt.savefig(f'{self.param_directory}/{self.param_level}confmatrix.png', dpi=200, bbox_inches='tight')
        plt.close()
        
    # def save_df(self, df, i):
    #     df.loc(i, 'Precision') = self.precision
    #     df.loc(i, 'Recall') = self.recall
    #     df.loc(i, 'FPR') = self.fpr
    #     df.loc(i, 'Fbeta') = self.fbeta_measure
        
def plot_pr(result_df, param_directory):
    fig, ax = plt.subplots()
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(f'Precision-Recall Curve for {param_directory}')
    for i, param in enumerate(result_df['Param Level']):
        ax.plot(result_df.loc[i, 'Recall'], result_df.loc[i, 'Precision'], marker='o', markersize=5, label=param)
    #ax.legend(loc='lower left')
    plt.savefig(f'{param_directory}/pr_curve.png', dpi=200, bbox_inches='tight')
    plt.close()
    
def sorting_ranking(result_df, param_directory):
    sorted_df = result_df.sort_values(by=['Precision', 'Fbeta'])
    sorted_df.to_csv(f'{param_directory}/ScoreResults.csv')
    
    
        
def plotting_results(param_directory, launches_df_list, param_level):
    print(param_directory)
    result_df = pd.DataFrame(columns=['Param Level', 'Precision', 'Recall', 'FPR', 'Fbeta'])
    #result_list = []
    for i, launches_df in enumerate(launches_df_list):
        #making launch and save directories, and initializing class
        launch_directory = os.path.join(param_directory, launches_df)
        plottest = Assessing_Data(launch_directory, param_directory, param_level[i])
        #doing confusion matrix for param level
        plottest.save_launch_and_obs_arrays()
        plottest.make_confusion_matrix()
        plottest.save_parameter_scores()
        #plottest.save_df(result_df, i)
        plottest.plot_confusion_matrix()
        result_df.loc[i, 'Param Level'] = param_level[i]
        result_df.loc[i, 'Precision'] = plottest.precision
        result_df.loc[i, 'Recall'] = plottest.recall
        result_df.loc[i, 'FPR'] = plottest.fpr
        result_df.loc[i, 'Fbeta'] = plottest.fbeta_measure
        print(f'confusion matrix {i}/{len(launches_df_list)} done')
        #result_list.append([param_level[i], plottest.precision, plottest.recall, plottest.fpr, plottest.fbeta_measure])
    #result_df[['Param Level', 'Precision', 'Recall', 'FPR', 'Fbeta']] = result_list
    plot_pr(result_df, param_directory)
    sorting_ranking(result_df, param_directory)
    
    
if __name__=='__main__':
    csv = 'TryingAlmostEverything_EvenLess/ScoreResults.csv'
    df = pd.read_csv(csv)
    print(np.max(df['Fbeta']))
    print(np.where(df['Fbeta']==np.max(df['Fbeta'])))
    
        


        
        
        