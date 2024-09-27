# -*- coding: utf-8 -*-
"""
@Author: Bienfait K. Isamura
         Quantum Chemical Topology (QCT) research group
         Department of Chemistry, UoM, UK
         -----------------------------------------------

This module performs a passive sampling  (PS).
Every PS object is specified by four mandatory attributes, namely:
    1. filename: name of the input file containing the ref data
    2. output_prop: name of the output property
    3. atom: atom of interest
    4. init_method: initialisation method 

The PS() class encompasses several methods that act on each instance.
The sampling is designed in such a way to prevent any sort of overlap
between the training, validation and test sets. The results are written
out in csv files. 

"""
import os
import numpy as np
import statistics
import random
import math
import sys
from utils.read_module import readfile
from utils.stratifiers import sturges, scott, rice, doane, fd, EpB
import random
from samplers.uncertainty import uncertainty_sampling
import shutil
import csv


    
class PS():   #Passive Sampling
    
    def __init__(self,filename, output_prop, atom,init_method="random", uncertainty_method="pairwise",*args):
        
        self.filename        = filename
        self.prop            = output_prop
        self.init_method     = init_method
        self.uncertainty_method = uncertainty_method
        self.atom_of_interest = atom
        #self.tr_set_size    = tr_size
        self.training_set   = None
        self.validation_set = None
        self.test_set = None
        self.pop_size = None
        self.uncertainty_curve_copy = None
        self.uncertainty_file_written = False
        for ar in args:
            if isinstance(ar,int):
                self.val_set_size = ar
            else:
                print("Unexpected optional argument encountered!")   
    
    def set_training_set(self,training_set_ID):
        self.training_set = training_set_ID    

    def get_training_set(self):
        return  self.training_set     

    def set_validation_set(self,validation_set_ID):
        self.validation_set = validation_set_ID    

    def get_validation_set(self):
        return  self.validation_set     

    def set_test_set(self,test_set_ID):
        self.test_set = test_set_ID    

    def get_test_set(self):
        return  self.test_set     

    def set_uncertainty_curve_copy(self,uncertainty_data):
        self.uncertainty_curve_copy = uncertainty_data    

    def get_uncertainty_curve_copy(self):
        return  self.uncertainty_curve_copy     

    def get_training_point_IDs(self,size_train):
        if (self.training_set!=None and len(self.training_set)==size_train):
            result = self.training_set
        else:
            ref_data_filename = self.filename
            file_content = readfile(ref_data_filename,self.prop)
            pop_size = len(file_content[1])
            regional_indices = [x for x in range(pop_size)] # indices of full pop
            size_regional_sample = size_train  # training is done on the entire population
            job_uncertainty = uncertainty_sampling(ref_data_filename,regional_indices,pop_size,size_regional_sample,self.prop,\
                                                 self.init_method, self.uncertainty_method)
            result = job_uncertainty.perform_uncertainty_sampling(self.atom_of_interest,SM="PS",DTYPE="TRAIN")
            self.uncertainty_curve_copy = job_uncertainty.uncertainty_curve
            self.training_set = result
            self.pop_size = pop_size
            
        
        return result
    
    def get_validation_point_IDs(self,size_valid):
        
        # Get pop size
        if self.pop_size!=None:
            pop_size = self.pop_size
        else:
            ref_data_filename = self.filename
            file_content = readfile(ref_data_filename,self.prop)
            pop_size = len(file_content[1])
        
        # Get train set idx
        if self.training_set!=None:
            training_set_IDs = self.training_set
        else:
            training_set_IDs = []
            
        indices = [x for x in range(pop_size)]
        complementary = []
        for idx in indices:
            if idx not in training_set_IDs:
                complementary.append(idx)
                
        random.shuffle(complementary)
        result = complementary[:size_valid]
        self.validation_set = result
        
        return result
        
    def get_test_point_IDs(self,size_test):
        # Get pop size
        if self.pop_size!=None:
            pop_size = self.pop_size
        else:
            ref_data_filename = self.filename
            file_content = readfile(ref_data_filename,self.prop)
            pop_size = len(file_content[1])
        
        # Get train set idx
        if self.training_set!=None:
            training_set_IDs = self.training_set
        else:
            training_set_IDs = []
            
        # Get valid set idx
        if self.training_set!=None:
            validation_set_IDs = self.validation_set
        else:
            validation_set_IDs = []
            
        indices = [x for x in range(pop_size)]
        complementary = []
        for idx in indices:
            #if idx not in training_set_IDs and idx not in validation_set_IDs:
            if idx not in training_set_IDs:
                complementary.append(idx)
                
        random.shuffle(complementary)
        result = complementary[:size_test]
        self.test_set = result
        
        return result
        
        
                
    def build_data_set(self,N,type_="Train", all_prop=False, current_prop=None):
        dataset             = []
        if (current_prop==None):
            input_file_content  = readfile(self.filename,self.prop)
        else:
            input_file_content  = readfile(self.filename,current_prop)
        input_data          = input_file_content[1]
        feats_indices       = input_file_content[3]
        prop_index          = input_file_content[4]
        all_prop_ind        = input_file_content[5]
        
        if (type_.upper()=="TRAIN"):
            points_IDs      = self.get_training_point_IDs(N)
            #Sometimes, there are slightly more points picked up (<1% for 1000 pts)
            #We truncate the sampled population
            if (len(points_IDs)>N):
                points_IDs = points_IDs[:N]            
        elif (type_.upper()=="VALID"):
            points_IDs      = self.get_validation_point_IDs(N)
            #Sometimes, there are slightly more points picked up (<1% for 1000 pts)
            #We truncate the sampled population
            if (len(points_IDs)>N):
                points_IDs = points_IDs[:N] 
                
        elif (type_.upper()=="TEST"):
            points_IDs      = self.get_test_point_IDs(N)
            #Sometimes, there are slightly more points picked up (<1% for 1000 pts)
            #We truncate the sampled population
            if (len(points_IDs)>N):
                points_IDs = points_IDs[:N] 
        else:
            sys.exit(" Dataset type can only be 'T' (training) or 'V'(Validation)")
        k=0
        for i in points_IDs:
            splitted_input_data=input_data[i].split(',')
            dataset.append([])
            for j in feats_indices: 
                dataset[k].append(splitted_input_data[j])
            
            if (all_prop==True):
                for index in all_prop_ind:
                    dataset[k].append(splitted_input_data[index])
            else:
                dataset[k].append(splitted_input_data[prop_index])
            k+=1
            
        return dataset
    def write_data_set(self, output_filename, sample_size, type_="Train", system_name="",write_all_prop=False,current_prop=None):
        if (current_prop==None):
            input_file        = readfile(self.filename,self.prop)
        else:
            input_file        = readfile(self.filename,current_prop)
        header            = input_file[0]
        feats_indices     = input_file[3] 
        prop_index        = input_file[4]
        all_prop_ind      = input_file[5]
        
        #Setting the header of the output file
        csv_header        = [header.split(',')[x] for x in feats_indices]
        if (write_all_prop):
            for index in all_prop_ind:
                csv_header.append(header.split(',')[index]) 
        else:
             csv_header.append(header.split(',')[prop_index])            
             
        if prop_index == len(header.split(",")) - 1:
            csv_header[-1] = csv_header[-1][:-1]
        csv_header = [str(item.split("_")[0]) for item in csv_header] 
        
               
        
        if (type_.upper()=="TRAIN"):
            data_set = self.build_data_set(sample_size,type_="Train",all_prop=write_all_prop)
        elif (type_.upper()=="VALID"):
            data_set = self.build_data_set(sample_size,type_="Valid",all_prop=write_all_prop)
        else:
            data_set = self.build_data_set(sample_size,type_="Test",all_prop=write_all_prop)

        
        for i in range(len(data_set)):
            data_row = data_set[i]
            if prop_index == len(header.split(","))-1:
                data_row[-1]=data_row[-1][:-1]
                data_set[i]=data_row
            data_set[i] = [eval(x) for x in data_set[i]]
        
        data = []
        data.append(csv_header)
        for row in data_set:
            data.append(row)


        # Check output folder
        outdir = os.path.split(output_filename)[0]
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        with open(output_filename,"w") as myfile:
            csv_writer = csv.writer(myfile)
            csv_writer.writerows(data) 
                        
        if type_.upper()=="TRAIN" and not self.uncertainty_file_written:
            output_dir = os.path.join(os.getcwd(),"Uncertainty-curves-PS-"+system_name)
            if not os.path.isdir(output_dir) and not self.uncertainty_file_written:
                os.mkdir(output_dir)
            
            uncertainty_file = os.path.join(output_dir,system_name+"_"+self.atom_of_interest+"_PS_uncertainty.dat")
            with open(uncertainty_file,"w") as uncertainty_curv:
                uncertainty_curv.write("Point  Uncertainty\n")
                for i in range(len(self.uncertainty_curve_copy)):
                    if i!=len(self.uncertainty_curve_copy)-1:
                        uncertainty_curv.write(f"{i+1}  {self.uncertainty_curve_copy[i]}\n")
                    else:
                        uncertainty_curv.write(f"{i+1}  {self.uncertainty_curve_copy[i]}")
            
            self.uncertainty_file_written = True

        
