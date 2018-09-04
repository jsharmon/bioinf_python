#!/usr/bin/env python
"""
Created on Wed Mar 22 12:56:51 2017

@author: jsharmon
"""

import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import pandas as pd

#All of these use the matplotlib package to create the bar plots; look up
#matplotlib documentation online
def make_intra_class_plots(uniqueClassList, data1, saveDir):
    ind_lipid_data = data1
    lipid_names = []
    lipid_mean_vals_c = []
    lipid_mean_vals_d = []
    
    for x in xrange(0, len(uniqueClassList)):
        for i in xrange(0, len(ind_lipid_data)):
            lipid_class = ind_lipid_data['Lipid_Names'][i].split(' ')[0]
            if(lipid_class == uniqueClassList[x]):
                stored_lipid_class = lipid_class
                lipid_names.append(data1['Lipid_Names'][i])
                lipid_mean_vals_c.append(data1['Control_Mean_Intensity'][i])
                lipid_mean_vals_d.append(data1['Disease_Mean_Intensity'][i])
                
        objects = tuple(lipid_names)
        y_pos = np.arange(len(objects))
        plt.bar(y_pos, lipid_mean_vals_c, align = 'center', alpha = 0.5,\
                color = 'blue', label = 'Control Samples')
        plt.xticks(y_pos, objects, rotation = 65)
        plt.ylabel('Concentration (mean mol %)')
        plt.title('Distribution of' + stored_lipid_class +\
                  'in Control vs Disease Samples')
        plt.tight_layout()
        
        plt.bar(y_pos, lipid_mean_vals_d, align = 'center', alpha = 0.5,\
                color = 'orange', label = 'Disease Samples')
        plt.xticks(y_pos, objects, rotation = 65)
        plt.ylabel('Concentration (mean mol %)')
        plt.title('Distribution of ' + stored_lipid_class +\
                  ' in Control vs Disease Samples')
        plt.legend(loc = 'upper left')
        plt.tight_layout()
        
        plt.savefig(saveDir + stored_lipid_class + 'plot.png')
        plt.close()
        stored_lipid_class = ''
        lipid_names = []
        lipid_mean_vals_c = []
        lipid_mean_vals_d = []

def make_inter_class_plots(uniqueClassList, data1, saveDir):

    classPlotObjects = tuple(uniqueClassList)
    classPlotDat = data1
    classPlotResponse_c = classPlotDat['Control Averages']
    classPlotResponse_d = classPlotDat['Disease Averages']
    y_pos_1 = np.arange(len(classPlotObjects))


    plt.bar(y_pos_1, classPlotResponse_c, align = 'center',\
                      alpha = 0.5, color='blue', label='Control Samples')
    plt.xticks(y_pos_1, classPlotObjects, rotation=65)
    plt.ylabel('Concentration (mean mol %)')
    plt.title('Concentration Values for Lipid Classes, Control vs Disease')
    plt.tight_layout()

    y_pos_1 = np.arange(len(classPlotObjects))
    plt.bar(y_pos_1, classPlotResponse_d, align = 'center', alpha = 0.5,\
            color='orange', label='Disease Samples')
    plt.xticks(y_pos_1, classPlotObjects, rotation=65)
    plt.ylabel('Concentration (mean mol %)')
    plt.title('Concentration Values for Lipid Classes, Control vs Disease')
    plt.tight_layout()
    plt.legend(loc = 'upper left')
    plt.savefig(saveDir + 'Inter_Class_Plot.png', dpi=150)
    plt.close()

def which_less_than(v = [], numOrIndex = '', val = 0):
    numValsInList = 0
    indexes = []
    for i in xrange(0, len(v)):
        if(v[i] < val):
            indexes.append(i)
            numValsInList += 1
    if(numOrIndex == 'num'):
        return numValsInList
    elif(numOrIndex == 'index'):
        return indexes
    else:
        return 'Specify Index or Value'

def make_class_pval_plot(uniqueClassList, saveDirCSV, saveDir):
    objects = tuple(uniqueClassList)
    y_pos = np.arange(len(objects))
    pvalDat = pd.read_csv(saveDirCSV + 'Class_Lipid_Statistics.csv')
    responseVar = pvalDat['Class P-Values']
    
    sig_pvals = which_less_than(responseVar, 'index', 0.05)
    
    barlist = plt.bar(y_pos, responseVar, align='center', alpha=0.5, color='green')
    
    for x in xrange(0, len(sig_pvals)):
        barlist[sig_pvals[x]].set_color('r')
        
    plt.xticks(y_pos, objects, rotation = 65)
    plt.ylabel('P-Value')
    plt.title('P-values by Lipid Class')
    plt.tight_layout()
    plt.savefig(saveDir + 'Class_P-Value_Plot.png', dpi=150)
    plt.close()

def make_intra_class_pval_plots(uniqueClassList, saveDirCSV, saveDir):
    ind_lipid_data = pd.read_csv(saveDirCSV + 'Individual_Lipid_Statistics.csv')
    lipid_names = []
    lipid_pvals = []
    
    for x in xrange(0, len(uniqueClassList)):
        for i in xrange(0, len(ind_lipid_data)):
            lipid_class = ind_lipid_data['Lipid_Names'][i].split(' ')[0]
            if(lipid_class == uniqueClassList[x]):
                stored_lipid_class = lipid_class
                lipid_names.append(ind_lipid_data['Lipid_Names'][i])
                lipid_pvals.append(ind_lipid_data['Individual P-Values'][i])
                
        objects = tuple(lipid_names)
        y_pos = np.arange(len(objects))
        
        
        barlist = plt.bar(y_pos, lipid_pvals, align = 'center', alpha = 0.5,\
                color = 'green')
        sig_pvals = which_less_than(lipid_pvals, 'index', 0.05)
        for x in xrange(0, len(sig_pvals)):
            barlist[sig_pvals[x]].set_color('r')
        
        plt.xticks(y_pos, objects, rotation = 65)
        plt.ylabel('P-Value')
        plt.title('P-Value for ' + stored_lipid_class +\
                  ' between Control and Disease Samples')
        plt.tight_layout()
        
        plt.savefig(saveDir + stored_lipid_class + ' p-value plot.png')
        plt.close()
        
        stored_lipid_class = ''
        lipid_names = []
        lipid_pvals = []
    
    
    
    
    
    
    
    
    
    
    
    