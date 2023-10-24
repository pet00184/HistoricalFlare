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
    ''' Redoing plotting again now that a lot of the data anlysis and scoring is being done real-time! We now have both
    "good egg" and "bad egg" parameters, for an initial separation between things that are useful and things that have
    too low of a recall or cancellation score. 
    
    Now have a folder of confuions matrices, Launch histograms, Cancellation histograms and the PRTTL plot!! (do one 
    2D and one 3D)
    '''
    
    def __init__(self, combo_results_df, bad_results_df, directory):
        self.combo_params = pd.read_csv(combo_results_df)
        self.bad_eggs = pd.read_csv(bad_results_df)
        self.main_directory = directory
        self.save_directory = os.path.join(directory, 'GoodEggParams')
        if not os.path.exists(os.path.join(self.save_directory, 'ConfMatrix')):
            os.makedirs(os.path.join(self.save_directory, 'ConfMatrix'))
        
    def make_all_plots(self):
        for i, param_combo in enumerate(self.combo_params['Param Combo']):
            self.make_confusion_matrix(i)
            self.plot_confusion_matrix(param_combo, i)
        self.plot_pr_2d()
        self.plot_pr_3d()
        self.plot_pr_2d(include_bad=True)
        
    def make_confusion_matrix(self, i):
        self.confusion_matrix = np.array([[self.combo_params.loc[i, 'TN'], self.combo_params.loc[i, 'TN_canc'], 
                        self.combo_params.loc[i, 'FP_noc5_noobs'], self.combo_params.loc[i, 'TP_noc5_obs']],
                        [self.combo_params.loc[i, 'FN'], self.combo_params.loc[i, 'FN_canc'], 
                        self.combo_params.loc[i, 'FP_c5_noobs'], self.combo_params.loc[i, 'TP']]])
        
    def plot_confusion_matrix(self, param_combo, i):
        ''' This just plots the new and improved confusion matrix for the specific parameter!
        '''
        fig, ax = plt.subplots()
        ax.matshow(self.confusion_matrix, cmap=plt.cm.Blues, alpha=0.7)
        for k in range(self.confusion_matrix.shape[0]):
            for j in range(self.confusion_matrix.shape[1]):
                ax.text(x=j, y=k,s=self.confusion_matrix[k, j], va='center', ha='center', size='xx-large')
        ylabels = ['', 'False', 'True']
        ax.set_ylabel('Flare above C5?', fontsize=14)
        ax.set_yticklabels(ylabels, rotation=90, fontsize=12)
        ax.xaxis.set_ticks_position("bottom")
        xlabels = ['', 'No Trigger', 'Cancelled Launch','Launch \n No C5 Obs', 'Launch \n C5 Obs']
        ax.set_xlabel('Launch Results', fontsize=14)
        ax.set_xticklabels(xlabels, fontsize=12)
        plt.suptitle(f'Parameters: {param_combo}')
        plt.title(f"Precision = {self.combo_params.loc[i, 'Precision']:.2f}, Recall = {self.combo_params.loc[i, 'Recall']:.2f}, LTT Ratio = {self.combo_params.loc[i, 'Trigger-to-Launch']:.2f}")
        plt.savefig(f'{self.save_directory}/ConfMatrix/{param_combo}confmatrix.png', dpi=200, bbox_inches='tight')
        plt.close()
        
    def plot_pr_2d(self, include_bad=False):
        fig, ax = plt.subplots()
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.set_title(f'Precision-Recall Curve for {self.main_directory}')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axvspan(0, 0.4, color='red', alpha=.2)
        getting_cmap = cm.get_cmap('viridis', len(self.combo_params['Param Combo']))
        color_list = getting_cmap(np.linspace(0, 1, len(self.combo_params['Param Combo'])))
        for i, param_combo in enumerate(self.combo_params['Param Combo']):
            ax.plot([self.combo_params.loc[i, 'Recall']], [self.combo_params.loc[i, 'Precision']], marker='o', markersize=5, label=param_combo, color=color_list[i])
        if include_bad: 
            for i, param_combo in enumerate(self.bad_eggs['Param Combo']):
                ax.plot([self.bad_eggs.loc[i, 'Recall']], [self.bad_eggs.loc[i, 'Precision']], marker='o', markersize=5, label='_', color='red')
        if len(self.combo_params['Param Combo']) < 10:        
            ax.legend(loc='lower left', fontsize=10)
        if include_bad:
            plt.savefig(f'{self.main_directory}/pr_all2D.png', dpi=200, bbox_inches='tight')
        if not include_bad:
            plt.savefig(f'{self.save_directory}/pr_2D.png', dpi=200, bbox_inches='tight')
        
    def plot_pr_3d(self):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlim(0.4, 1)
        ax.set_ylim(0.4, 1)
        ax.set_zlim(0, 1)
        ax.set_xlabel('Recall')
        ax.set_ylabel('LTT Ratio')
        ax.set_zlabel('Precision')
        #ax.axvspan(0, 0.4, color='red', alpha=.2)
        for i, param_combo in enumerate(self.combo_params['Param Combo']):
            ax.plot(self.combo_params.loc[i, 'Recall'], self.combo_params.loc[i, 'Trigger-to-Launch'], self.combo_params.loc[i, 'Precision'], marker='o', markersize=5, label=param_combo)
        plt.savefig(f'{self.save_directory}/pr_3D.png', dpi=200, bbox_inches='tight')
    
    
        


        
        
        