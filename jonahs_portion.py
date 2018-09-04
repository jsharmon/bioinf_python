#!/usr/bin/env python
#The line above lets you run this from the command line
'''
This code will take inputData.csv and 
make another csv containing the average intensity values, average stdev of
intensity values, and number of lipids by class for both control and disease
groups.

For testing, use inputData.csv from the desktop as the sys argument
'''

'''
This is how to import packages. Could just import them without the 'as' part,
but it's easier to import them as fewer letters so they're easier to type later
on. You'll also probably want to import math, and:
    from __future__ import division
sys is imported to read inputs from the command line, and os is imported to 
make the folders to organize the output.
'''
import pandas as pd
import control_vs_disease_ttests as cvdt
import individual_lipid_ttests as ilt
import sys
import making_plots as mp
import os

#Assumes the lipids have class followed by a space, then other data
'''
This is how to write your own function in python. You can call it later in the
file and it will reference this part and run through the code that it contains
'''
def get_classes(dataFrame):
    #Look at the pandas documentation online for more info about data frames
    lipid_names = dataFrame['Lipid_Names']
    #this created an empty list - looks at python documentation for info about
    #lists
    class_names = []
    for val in xrange(0, len(lipid_names)):
        #This splits the input string on spaces, then takes the first piece -
        #the lipid class name
        class_names.append(lipid_names[val].split(' ')[0])
        
    return class_names

#This part takes inputs from the command line
input_file = sys.argv[1]
numSamples = int(sys.argv[2])

#This takes the input file path and creates a new variable called saveDir
#to assist with saving all of the output to the same location as the input file
saveDirList = input_file.split('\\')[0:-1]
saveDir = '\\'.join(saveDirList) + '\\'

#Makes a bunch of folders to organize the output
plotFolderPath = saveDir + 'Plots' 
if not os.path.exists(plotFolderPath):
    os.makedirs(plotFolderPath)
    
saveDirPlots = plotFolderPath + '\\'

plotFolderPath1 = saveDir + 'CSV Output' 
if not os.path.exists(plotFolderPath1):
    os.makedirs(plotFolderPath1)
    
saveDirCSV = plotFolderPath1 + '\\'

plotFolderPath2 = saveDirPlots + 'Intra_Class_Plots' 
if not os.path.exists(plotFolderPath2):
    os.makedirs(plotFolderPath2)
    
saveDirIntraClassPlots = plotFolderPath2 + '\\'

plotFolderPath3 = saveDirPlots + 'P-Value_Plots' 
if not os.path.exists(plotFolderPath3):
    os.makedirs(plotFolderPath3)

saveDirPValuePlots = plotFolderPath3 + '\\'

#Reading in the data as a dataframe using pandas
lipid_data = pd.read_csv(input_file)
#Creating a list of the unique class names to use for finding averages
class_list = get_classes(lipid_data)
lipid_data['class'] = pd.Series(class_list, index = lipid_data.index)

uniqueClassList = cvdt.unique(class_list)
#Creating an empty list as long as the class number to account for any number
#of lipid classes
empty_vals = 0.0 * len(uniqueClassList)

'''
This section creates a data frame using pandas. This data frame stores the
average values. avg_class_vals will later be output to a csv file.
'''
    
avg_class_vals = pd.DataFrame({
    'Control Averages' : 
    pd.Series(empty_vals, index = uniqueClassList),
    'Control StErrors' :
    pd.Series(empty_vals, index = uniqueClassList),
    'Control Number' :
    pd.Series(empty_vals, index = uniqueClassList),
    'Disease Averages' : 
    pd.Series(empty_vals, index = uniqueClassList),
    'Disease StErrors' :
    pd.Series(empty_vals, index = uniqueClassList),
    'Disease Number' :
    pd.Series(empty_vals, index = uniqueClassList)})

def get_mean_sd_num(control_or_disease = 'c'):
    '''
    This sums the intensity values for each class and stores them in the
    corresponding class row in intensity_class_vals. It also stores the number of
    lipids in each class that are seen. The first loop runs through the mean values,
    then stores them in avg_class_vals; the second runs through sterror values.
    This will separate the lipids based on which sample they are from (control
    or disease).
    
    Arguments:
        control_or_disease - specify whether the sample is from a control
        or diseased individual
    Returns:
        Nothing; this modifies the existing data frames created earlier
    '''
    
    '''
    The following sections create data frames using pandas that store
    intensity values (intensity_class_vals), and number of lipids of a 
    particular class (num_class_vals)
    '''
    if(control_or_disease == 'c'):
        modifier = 'Control '
    elif(control_or_disease == 'd'):
        modifier = 'Disease '
    else:
        return
    
    intensity_class_vals = pd.DataFrame({
    'Intensity' : 
    pd.Series(empty_vals, index = uniqueClassList)})
    
    num_class_vals = pd.DataFrame({
    'Number' : 
    pd.Series(empty_vals, index = uniqueClassList)})         
    
    #For the mean values
    for i in xrange(0, len(lipid_data['class'])):
        for x in xrange(0, len(uniqueClassList)):
            if((lipid_data['class'][i] == uniqueClassList[x])\
               and (modifier == 'Control ')):
                intensity_class_vals.loc[uniqueClassList[x]] += \
                        lipid_data['Control_Mean_Intensity'][i]
                num_class_vals.loc[uniqueClassList[x]] += 1
            elif((lipid_data['class'][i] == uniqueClassList[x])\
               and (modifier == 'Disease ')):
                intensity_class_vals.loc[uniqueClassList[x]] += \
                        lipid_data['Disease_Mean_Intensity'][i]
                num_class_vals.loc[uniqueClassList[x]] += 1

    for x in xrange(0, len(intensity_class_vals)):
        avg_class_vals[modifier + 'Averages'][x] = \
                      intensity_class_vals['Intensity'][x]
 
    intensity_class_vals = pd.DataFrame({
    'Intensity' : 
    pd.Series(empty_vals, index = uniqueClassList)})
    
    #For the StError values
    for i in xrange(0, len(lipid_data['class'])):
        for x in xrange(0, len(uniqueClassList)):
            if((lipid_data['class'][i] == uniqueClassList[x])\
               and (modifier == 'Control ')):
                intensity_class_vals.loc[uniqueClassList[x]] += \
                        lipid_data['Control_SE_Intensity'][i]
            elif((lipid_data['class'][i] == uniqueClassList[x])\
               and (modifier == 'Disease ')):
                intensity_class_vals.loc[uniqueClassList[x]] += \
                        lipid_data['Disease_SE_Intensity'][i]
        
    for x in xrange(0, len(intensity_class_vals)):
        avg_class_vals[modifier + 'StErrors'][x] = \
                      intensity_class_vals['Intensity'][x]
    
    for x in xrange(0, len(num_class_vals)): 
        avg_class_vals[modifier + 'Number'][x] = num_class_vals['Number'][x]
    
    
#Get control and disease values for csv
get_mean_sd_num('c')
get_mean_sd_num('d')

#This outputs avg_class_vals to a csv file. It contains the mean and SE of
#intensity values for each class of lipids, and the number of lipids in each 
#class, separated into control and disease categories
avg_class_vals.to_csv(saveDirCSV + 'Lipid_Class_Data.csv')

'''
This uses all of the code in the other files to run through other functions
written in those files
'''
cvdt.make_pvals_csv(avg_class_vals, lipid_data, numSamples, saveDirCSV)

ilt.run_tests(input_file, numSamples, saveDirCSV)

mp.make_inter_class_plots(uniqueClassList, avg_class_vals, saveDirPlots)

mp.make_class_pval_plot(uniqueClassList, saveDirCSV, saveDirPValuePlots)

mp.make_intra_class_plots(uniqueClassList, lipid_data, saveDirIntraClassPlots)

mp.make_intra_class_pval_plots(uniqueClassList, saveDirCSV, saveDirPValuePlots)