#!/usr/bin/env python
'''

'''

import scipy.stats as stats
import pandas as pd
import math
from individual_lipid_ttests import conf_int_95percent

'''
This function takes a sequence as a list and outputs a list with only the
unique values, while preserving the initial order of the input. Credit to 
Markus Jarderot, Kevin Guan on Stack Overflow.

Argument:
    seq - the list from which you want to extract unique values
Returns:
    A list of only the unique values, in the original order
'''
def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def make_pvals_csv(data1, data2, NS, saveDir):
    cvd_data = data1
    class_empty_vals = 0.0 * len(cvd_data)
    #Making an empty data frame to store the p-values
    class_pvals = pd.DataFrame({
                'Class P-Values': \
                pd.Series(class_empty_vals, index = unique(list(data2['class']))),
                'Control Sample Size': \
                pd.Series(class_empty_vals, index = unique(list(data2['class']))),
                'Disease Sample Size': \
                pd.Series(class_empty_vals, index = unique(list(data2['class']))),
                'Lower Bound of 95% CI': \
                pd.Series(class_empty_vals, index = unique(list(data2['class']))),
                'Upper Bound of 95% CI': \
                pd.Series(class_empty_vals, index = unique(list(data2['class'])))})
    #This is the number of samples from the original experiment; need to use this
    #for the t-test to account for actual sample size in each class. In the paper
    #we are using, numSamples = 10 for disease and control. 
    #We will probably need to set this as a sys.argv input
    numSamples = NS
    #Going to run a t-test for each class between control and disease, and store
    #the p-values into pvals
    for x in xrange(0, len(cvd_data)):
        class_pvals['Control Sample Size'][x] = numSamples
        class_pvals['Disease Sample Size'][x] = numSamples
        
        c_mean = float(cvd_data['Control Averages'][x])
        c_sd = float(cvd_data['Control StErrors'][x]) * \
                    math.sqrt(numSamples)
                    
        d_mean = float(cvd_data['Disease Averages'][x])
        d_sd = float(cvd_data['Disease StErrors'][x]) * \
        math.sqrt(numSamples)
                    
        n1 = numSamples
        n2 = numSamples
    
        t_stat, p_val = stats.ttest_ind_from_stats(c_mean, c_sd, \
                            n1, d_mean, d_sd, n2, equal_var=True)
        
        lb, ub = conf_int_95percent(c_mean, c_sd, NS, d_mean, d_sd, NS)
        #The numbers were flipped for some reason
        class_pvals['Lower Bound of 95% CI'][x] = ub
        class_pvals['Upper Bound of 95% CI'][x] = lb
        
        class_pvals['Class P-Values'][x] = p_val
        if((n1 < 3) or (n2 < 3)):
            class_pvals['Class P-Values'][x] = 1
        

    #This outputs pvals to a csv of p-values with classes as the labels
    class_pvals.to_csv(saveDir + 'Class_Lipid_Statistics.csv')