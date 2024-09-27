# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 18:48:09 2022

@author: q63762bi
"""
import os
import sys
import numpy as np
import math
import statistics

def sturges(population):
    population_size   = len(population)
    number_of_bins    = math.ceil(1+math.log(population_size,2))
    
    if (number_of_bins<2):
        sys.exit("Stratification impossible due to unsufficient number of bins")
        
    return number_of_bins

def doane(population):
    pop_size       = len(population)
    mean           = statistics.mean(population)
    g1_numer       = 0
    g1_denom       = 0
        
    for value in population:
        g1_numer+= np.power(value-mean,3)/pop_size
        g1_denom+= np.power(value-mean,2)/pop_size
        
    g1_denom       = np.power(g1_denom,3/2)
    g1             = g1_numer/g1_denom
    sigma_g1_numer = 6*(pop_size-2)
    sigma_g1_denom = (pop_size+1)*(pop_size+3)
    sigma_g1       = np.power(sigma_g1_numer/sigma_g1_denom,0.5)
    factor         = math.log(pop_size*(1+abs(g1)/sigma_g1),2)
    number_of_bins = math.ceil(1+factor)
    
    #print (f"g1={g1}")
    
    if (number_of_bins<2):
        sys.exit("Stratification impossible due to unsufficient number of bins")
        
    return number_of_bins

def scott(population):
    population_size   = len(population)
    std               = statistics.stdev(population)
    bin_width         = 3.49*std/(np.power(population_size,1/3))
    number_of_bins    = math.ceil((max(population)-min(population))/bin_width)
    
    #print (f"STD={std}")
    #print ("Range=",max(population)-min(population))
    if (number_of_bins<2):
        sys.exit("Stratification impossible due to unsufficient number of bins")
        
    return number_of_bins

def fd(population):
    population_size   = len(population)
    q3, q1            = np.percentile(population, [75 ,25])
    iqr               = q3 - q1
    bin_width         = 2*iqr/(np.power(population_size,1/3))
    range_            = max(population)-min(population)
    number_of_bins    = math.ceil(range_/bin_width)
    
    #print (f"IQR={iqr}")
    #print (f"Range={range_}")
    #print (f"Bin width={bin_width}")
    if (number_of_bins<2):
        sys.exit("Stratification impossible due to unsufficient number of bins")
        
    return number_of_bins

def rice(population):
    population_size   = len(population)
    number_of_bins    = math.ceil(2*np.power(population_size,1/3))
    
    if (number_of_bins<2):
        sys.exit("Stratification impossible due to unsufficient number of bins")
        
    return number_of_bins

def EpB(population):
    population_size   = len(population)                
    number_of_bins    = math.ceil(np.power(2*population_size,2/5))
    
    if (number_of_bins<2):
        sys.exit("Stratification impossible due to unsufficient number of bins")
        
    return number_of_bins    

def get_strat_properties(population,strat):
    
    # Return number_of_bins, range, IQR, min, max,
    # Std, median, mean, kurtosis
    number_of_bins = None
    min_dist = None
    max_dist = None
    IQRange = None
    range_dist= None
    Std_dist= None
    median_dist = None
    mean_dist = None
    kurtosis = None
    
    # Get number of bins
    if isinstance(strat,str):
        method = strat
        # Get number of bins
        if (method=="Sturges"):
            number_of_bins = sturges(population)
        elif (method=="Scott"):
            number_of_bins = scott(population)
        elif (method=="Rice"):
            number_of_bins = rice(population)
        elif (method=="Doane"):
            number_of_bins = doane(population)
        elif (method=="Equiprob" or method=="Equiprobable bins"):
            number_of_bins = EpB(population)
        elif (method=="Freedman-Diaconis" or method=="FD"):
            number_of_bins = fd(population)
        else:
            number_of_bins = fd(population)
    elif isinstance(strat,int):
        number_of_bins = strat
    else:
        sys.exit("Invalid stratification method")
        
        
        
    result = [number_of_bins, min_dist, max_dist, range_dist, IQRange, Std_dist,\
              median_dist, mean_dist, kurtosis]
    
    return result
        
