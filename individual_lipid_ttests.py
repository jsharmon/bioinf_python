#!/usr/bin/env python
"""
Take the directory from the input, use that, save all of the csv files to that
directory to keep everything in the same file
"""
from __future__ import division
import math
import pandas as pd
import scipy.stats as stats

#Will calculate the bounds of a 95% confidence interval to save to a csv
def conf_int_95percent(m1, sd1, n1, m2, sd2, n2):   
    df_top = (((sd1**2)/n1)+((sd2**2)/n2))**2
    df_bottom = (1/(n1-1))*((sd1**2)/n1)**2 + (1/(n2-1))*((sd2**2)/n2)**2
    df = df_top / df_bottom
    t_star = stats.t.ppf(0.025, df)
    
    mean_diff = m2 - m1
    error = math.sqrt(((sd1**2)/n1)+((sd2**2)/n2))
    lowerBound = mean_diff - (t_star * error)
    upperBound = mean_diff + (t_star * error)
    
    return lowerBound, upperBound
    
#Easier to run through this in the second file
def run_tests(fileName, NS, saveDir):
    #Reading in the data as a dataframe using pandas
    ind_lipid_data = pd.read_csv(fileName)
        
    #Generating an empty list comprised of zeroes that is the length of
    #ind_lipid_data
    ind_empty_vals = 0.0 * len(ind_lipid_data)
    
    lipid_names = ind_lipid_data['Lipid_Names']

    ind_pvals = pd.DataFrame({'Individual P-Values': \
                      pd.Series(ind_empty_vals, index = lipid_names),
                          'Control Sample Size': \
                      pd.Series(ind_empty_vals, index = lipid_names),
                          'Disease Sample Size': \
                      pd.Series(ind_empty_vals, index = lipid_names),
                          'Lower Bound of 95% CI': \
                      pd.Series(ind_empty_vals, index = lipid_names),
                          'Upper Bound of 95% CI': \
                      pd.Series(ind_empty_vals, index = lipid_names)})
    
    #Will likely want to include this as a sys.argv input; 10 samples from the paper
    #we are using
    numSamples = NS

    c_mean = []
    c_sd = []
    d_mean = []
    d_sd = []
    
    #assumes input is in the form of inputData.csv
    for x in xrange(0, len(ind_lipid_data)):
        c_mean.append(float(ind_lipid_data['Control_Mean_Intensity'][x]))
        c_sd.append(float(ind_lipid_data['Control_SE_Intensity'][x]) * \
                    math.sqrt(numSamples))
        
        d_mean.append(float(ind_lipid_data['Disease_Mean_Intensity'][x]))
        d_sd.append(float(ind_lipid_data['Disease_SE_Intensity'][x]) * \
                    math.sqrt(numSamples))
        

    for x in xrange(0, len(lipid_names)):
        n1 = numSamples
        n2 = numSamples
        m1 = c_mean[x]
        m2 = d_mean[x]
        sd1 = c_sd[x]
        sd2 = d_sd[x]
    
        #Uses the scipy.stats package to run t-tests
        t_stat, p_val = stats.ttest_ind_from_stats(m1, sd1, \
                        n1, m2, sd2, n2, equal_var=True)
    
        ind_pvals['Individual P-Values'][x] = p_val
        
        if((n1 < 3) or (n2 < 3)):
            ind_pvals['Individual P-Values'][x] = 1
            
        ind_pvals['Control Sample Size'][x] = n1
        ind_pvals['Disease Sample Size'][x] = n2
    
        if((sd1 == 0) | (sd2 == 0) | (m1 == 0) | (m2 == 0)):
            ind_pvals['Lower Bound of 95% CI'][x] = None
            ind_pvals['Upper Bound of 95% CI'][x] = None
        else:
            l, u = conf_int_95percent(m1, sd1, n1, m2, sd2, n2)
        
            #the lower and upper bound were flipped in the csv, corrected here
            ind_pvals['Lower Bound of 95% CI'][x] = u
            ind_pvals['Upper Bound of 95% CI'][x] = l
    
    #Uses pandas to save output to a csv file
    ind_pvals.to_csv(saveDir + 'Individual_Lipid_Statistics.csv')
        















