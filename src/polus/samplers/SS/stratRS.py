# -*- coding: utf-8 -*-
"""
@Author: Bienfait K. Isamura
         Quantum Chemical Topology (QCT) research group
         Department of Chemistry, UoM, UK
         -----------------------------------------------

This module implements the SRS (Stratified Random Sampling) class. 
Every SRS object is specified by three attributes, namely:
    1. filename: name of the input file containing the ref data
    2. output_prop: name of the output property with respect to
                    which the stratification is performed
    3. strat_method: the method to be used for stratification. 
                     Stratification methods are defined in the 
                     stratifiers module (here imported)

The SRS() class encompasses several methods that act on each instance.
Once the stratification has been performed, one can now sample the 
training, validation and test sets. This sampling returns the IDs 
(row number in the input filename) of the points sampled. The last
step is to write the different samples into csv files.

"""

import os
import time
import numpy as np
import statistics
import random
import math
import sys
import shutil
import csv
import random
from polus.utils.read_module import readfile
from polus.files.outputs import WriteJobDetails
from polus.samplers.SS.stratifiers import sturges, scott, rice, doane, fd, EpB
from polus.utils.logging import RaiseWarning
#from uncertainty import uncertainty_sampling



class SRS():  # SRS: Stratified Random Sampling
    
    def __init__(self,filename=None, valTest=False,inputDir=None,systemName="MOL",outputDir=None,trainSize=1000, valSize=1000, testSize=1000,output_prop="iqa", stratByWFN=True,writeAllProps=False,*args):
        
        self.filename       = filename
        self.systemName     = systemName
        self.inputDir       = inputDir
        self.outputDir      = outputDir
        self.prop           = output_prop
        self.stratByWFN     = stratByWFN
        self.trainSize      = trainSize
        self.validSize      = valSize
        self.testSize       = testSize
        self.writeAllProps  = writeAllProps
        self.valTest        = valTest
        self.inputFiles     = None
        self.strat_method   = None
        self.training_set   = None
        self.validation_set = None
        self.test_set       = None
        self.nbr_of_zones   = None

        # Check input and output directory
        if self.outputDir==None:
            self.outputDir=os.getcwd()
        else:
            if not os.path.isdir(self.outputDir):
                os.mkdir(self.outputDir)

        # Check optional arguments
        for ar in args:
            if isinstance(ar,int):
                self.nbr_of_zones = ar
            elif isinstance(ar,str):
                self.strat_method = ar
            else:
                print("Unexpected optional argument encountered!")
                sys.exit()
                
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

    def get_number_of_bins(self):       
        method          = self.strat_method
        file_content    = readfile(self.filename,self.prop)
        population     = list(file_content[2])        
        if ((self.nbr_of_zones == None) or (isinstance(self.strat_method,str))):        
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
        elif ((self.strat_method == None) and (isinstance(self.nbr_of_zones,int))):
             number_of_bins = self.nbr_of_zones
        else:
            sys.exit("Error: Cannot determine the number of bins! Job killed.")
            
        return number_of_bins
       
    def backup_training_set_ID(self, tr_sample_IDs):
        return tr_sample_IDs
    
    def stratify_population(self):
        if not self.stratByWFN:
            data               = readfile(self.filename,self.prop)[2]
        else:
            data               = readfile(self.filename,"wfn_energy")[2]
        nbr_of_zones       = self.get_number_of_bins()
        stratified_pop_IDs = []
        fractions          = []
        
        
        for i in range(nbr_of_zones):
            stratified_pop_IDs.append([])
            if (i!=0):
                fractions.append((i)*100/nbr_of_zones)
        
        percentiles = np.percentile(data,fractions)

        for j in range(len(data)):
            value=data[j]
            if value <= percentiles[0]:
                stratified_pop_IDs[0].append(j)
            else:
                flag=False
                for k in range(1,len(percentiles)):
                    if percentiles[k-1] < value <= percentiles[k]:
                        stratified_pop_IDs[k].append(j)
                        flag=True
                if flag==False:
                    stratified_pop_IDs[-1].append(j)
                    
        return stratified_pop_IDs
    
    def get_training_point_IDs(self, size_tr):
        nbr_of_zones            = self.get_number_of_bins()
        pop_size                = len(readfile(self.filename,self.prop)[2])
        sampled_tr_points_IDs   = []
        sampling_ratio          = size_tr/pop_size
        stratified_pop_IDs      = self.stratify_population()
        regions_size            = min([len(stratified_pop_IDs[i]) for i in range(nbr_of_zones)])
        local_sampling_attempts = math.ceil(sampling_ratio*regions_size)
        
        if (self.training_set==None): 
            if (local_sampling_attempts<=regions_size):
                for i in range(len(stratified_pop_IDs)):
                    cur_region=stratified_pop_IDs[i]
                    cur_reg_size=len(cur_region)
                    j=0
                    while j<local_sampling_attempts:
                        random.seed()
                        rnd_index = random.randint(0,cur_reg_size-1)
                        sampled_tr_points_IDs.append(cur_region[rnd_index])
                        j+=1
            else:
                sys.exit("Wrong stratification pattern.\n Please edit sample (regions) sizes)")
            remaining_points = size_tr - len(sampled_tr_points_IDs)
            if (remaining_points < 0):
                sampled_tr_points_IDs = sampled_tr_points_IDs[:remaining_points] 
            else:
                all_IDs = []
                for k in range(len(stratified_pop_IDs)):
                    for value in stratified_pop_IDs[k]:
                        all_IDs.append(value)
                for p in range(remaining_points):    
                    random.seed()
                    rnd_index = random.randint(0,len(all_IDs)-1)
                    if  rnd_index not in sampled_tr_points_IDs:
                        sampled_tr_points_IDs.append(rnd_index)
                    else:
                        retries = 1
                        while retries <20:
                            rnd_index = random.randint(0,len(all_IDs)-1)
                            if  rnd_index not in sampled_tr_points_IDs:
                                sampled_tr_points_IDs.append(all_IDs[rnd_index])
                                retries = 20
                            else:
                                retries +=1

            random.shuffle(sampled_tr_points_IDs)
            self.training_set = sampled_tr_points_IDs

        else:
            sampled_tr_points_IDs = self.training_set

        return sampled_tr_points_IDs
        
    def get_validation_point_IDs(self, size_val):
        nbr_of_zones            = self.get_number_of_bins()
        pop_size                = len(readfile(self.filename,self.prop)[2])
        sampled_val_points_IDs  = []
        sampling_ratio          = size_val/pop_size
        stratified_pop_IDs      = self.stratify_population()
        regions_size            = min([len(stratified_pop_IDs[i]) for i in range(nbr_of_zones)])
        local_sampling_attempts = math.ceil(sampling_ratio*regions_size)
        
        if (self.training_set==None):
            tr_points_IDs=[]
        else:
            tr_points_IDs = self.training_set
            
        if (self.validation_set==None):
            validation_points_IDs = []
            if (local_sampling_attempts<=regions_size):
                for i in range(len(stratified_pop_IDs)):
                    cur_region=stratified_pop_IDs[i]
                    cur_reg_size = len(cur_region)
                    j=0
                    while j<local_sampling_attempts:
                        random.seed()
                        rnd_index    = random.randint(0,cur_reg_size-1)
                        rnd_point_ID = cur_region[rnd_index]
                        if (rnd_point_ID not in tr_points_IDs) and \
                            rnd_point_ID not in validation_points_IDs:
                            sampled_val_points_IDs.append(rnd_point_ID)
                            validation_points_IDs.append(rnd_point_ID)
                        else:
                            flag=True
                            nretry=0
                            while flag:
                                rnd_index    = random.randint(0,regions_size-1)
                                rnd_point_ID = cur_region[rnd_index]
                                if rnd_point_ID not in tr_points_IDs and \
                                    rnd_point_ID not in validation_points_IDs:
                                    sampled_val_points_IDs.append(rnd_point_ID)
                                    validation_points_IDs.append(rnd_point_ID)
                                    flag=False
                                else:
                                    nretry+=1
                                    if (nretry==len(cur_region)):
                                        sys.exit("Could not build validation set")
                                    else:
                                        flag=True
                    
                        j+=1

            else:
                sys.exit("Wrong stratification pattern.\n Please edit sample (regions) sizes)")
            random.shuffle(sampled_val_points_IDs)
            self.validation_set = sampled_val_points_IDs
        else:
            sampled_val_points_IDs = self.validation_set
            
        return sampled_val_points_IDs
    
    def get_test_point_IDs(self, size_val):
        nbr_of_zones            = self.get_number_of_bins()
        pop_size                = len(readfile(self.filename,self.prop)[2])
        sampled_test_points_IDs  = []
        sampling_ratio          = size_val/pop_size
        stratified_pop_IDs      = self.stratify_population()
        regions_size            = min([len(stratified_pop_IDs[i]) for i in range(nbr_of_zones)])
        local_sampling_attempts = math.ceil(sampling_ratio*regions_size)
        
        if (self.training_set==None):
            tr_points_IDs=[]
        else:
            tr_points_IDs = self.training_set
            
        if (self.validation_set==None):
            validation_points_IDs=[]
        else:
            validation_points_IDs = self.validation_set
            
        if (self.test_set==None):
            test_points_IDs=[]
            if (local_sampling_attempts<=regions_size):
                for i in range(len(stratified_pop_IDs)):
                    cur_region=stratified_pop_IDs[i]
                    cur_reg_size = len(cur_region)
                    j=0
                    while j<local_sampling_attempts:
                        random.seed()
                        rnd_index    = random.randint(0,cur_reg_size-1)
                        rnd_point_ID = cur_region[rnd_index]
                        if (rnd_point_ID not in tr_points_IDs) and \
                            rnd_point_ID not in validation_points_IDs and \
                            rnd_point_ID not in test_points_IDs:
                            sampled_test_points_IDs.append(rnd_point_ID)
                            test_points_IDs.append(rnd_point_ID)
                        else:
                            flag=True
                            nretry=0
                            while flag:
                                rnd_index    = random.randint(0,regions_size-1)
                                rnd_point_ID = cur_region[rnd_index]
                                if rnd_point_ID not in tr_points_IDs and \
                                    rnd_point_ID not in test_points_IDs:
                                    sampled_test_points_IDs.append(rnd_point_ID)
                                    test_points_IDs.append(rnd_point_ID)
                                    flag=False
                                else:
                                    nretry+=1
                                    if (nretry==len(cur_region)):
                                        sys.exit("Could not build validation set")
                                    else:
                                        flag=True
                    
                        j+=1

            else:
                sys.exit("Wrong stratification pattern.\n Please edit sample (regions) sizes)")

            random.shuffle(sampled_test_points_IDs)
            self.test_set = sampled_test_points_IDs
        else:
            sampled_test_points_IDs = self.test_set

        return sampled_test_points_IDs
        
    def check_overlap(self):
        training   = self.training_set
        validation = self.validation_set
        if (training!=None and validation!=None):
            flag=False
            for i in training:
                for j in validation:
                    if (i==j):
                        flag=True
                        break
        else:
            sys.exit("--->Make sure you have generated the training and validation sets before")
            #time.sleep(0.5)
        if flag==False:
            result="--->No overlap between training and validation data sets"

        else:
            result=" There are overlapping between tr and val points. \
                     Please rerun sampling of both"
        
        print(result)
        #time.sleep(0.5)
        
        return result
        
    def build_data_set(self,N,type_="Train", all_prop=False):  
        
        dataset             = []
        input_file_content  = readfile(self.filename,self.prop)
        input_data          = input_file_content[1]
        feats_indices       = input_file_content[3]
        prop_index          = input_file_content[4]
        all_prop_ind        = input_file_content[5]
        
        if (type_.upper()=="TRAIN"):
            points_IDs      = self.get_training_point_IDs(N)
            if (len(points_IDs)>N):
                points_IDs = points_IDs[:N]            
        elif (type_.upper()=="VALID"):
            points_IDs      = self.get_validation_point_IDs(N)
            if (len(points_IDs)>N):
                points_IDs = points_IDs[:N] 
                
        elif (type_.upper()=="TEST"):
            points_IDs      = self.get_test_point_IDs(N)
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
            
    
    def write_data_set(self, output_filename, sample_size, type_="Train", write_all_prop=False):
        input_file        = readfile(self.filename,self.prop)
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
        
    def GetInputFiles(self):
        if self.inputFiles==None:
            self.inputFiles=list()
            for item in os.listdir(self.inputDir):
                if item.endswith(".csv"):
                    self.inputFiles.append(os.path.join(self.inputDir,item))               
         

    def GetAtomLabel(self,filename):
        elements = ["C","O","H","N","P","S","F","Cl"]
        digits   = [str(x) for x in range(10)]
        string   = os.path.split(filename)[1].split(".csv")[0]
        for entry in string.split("_"):
            if len(entry)>=2 and entry.isalnum() and entry[0] in elements and entry[1] in digits:
                atom = entry
                break
        return atom
 
    def Execute(self):
        atoms = list()
        self.GetInputFiles()
        file0         = self.inputFiles[0]
        print(f"POLUS: Processing {file0}")
        self.filename = file0
        atom       = self.GetAtomLabel(file0)
        outfile_TR = os.path.join(self.outputDir,self.systemName+"_"+atom.upper()+"_TRAINING_SET.csv")
        outfile_VL = os.path.join(self.outputDir,self.systemName+"_"+atom.upper()+"_INT_VALIDATION_SET.csv")
        outfile_TS = os.path.join(self.outputDir,self.systemName+"_"+atom.upper()+"_EXT_VALIDATION_SET.csv")
        self.write_data_set(output_filename=outfile_TR, sample_size=self.trainSize, type_="Train", write_all_prop=self.writeAllProps)
        self.write_data_set(output_filename=outfile_VL, sample_size=self.validSize, type_="Valid", write_all_prop=self.writeAllProps)
        self.write_data_set(output_filename=outfile_TS, sample_size=self.testSize, type_="Test", write_all_prop=self.writeAllProps)
        if self.valTest:
            for value in self.test_set:
                 self.validation_set.append(value)
            self.test_set  = self.validation_set.copy()
            self.validSize = len(self.validation_set)
            self.testSize  = len(self.test_set)
        trset = self.training_set
        vlset = self.validation_set
        tsset = self.test_set
        atoms.append(atom)
        for file_ in self.inputFiles[1:]:
            print(f"POLUS: Processing {file_}")
            atom          = self.GetAtomLabel(file_)
            self.filename = file_
            outfile_TR = os.path.join(self.outputDir,self.systemName+"_"+atom.upper()+"_TRAINING_SET.csv")
            outfile_VL = os.path.join(self.outputDir,self.systemName+"_"+atom.upper()+"_INT_VALIDATION_SET.csv")
            outfile_TS = os.path.join(self.outputDir,self.systemName+"_"+atom.upper()+"_EXT_VALIDATION_SET.csv")
            self.write_data_set(output_filename=outfile_TR, sample_size=self.trainSize, type_="Train", write_all_prop=self.writeAllProps)
            self.write_data_set(output_filename=outfile_VL, sample_size=self.validSize, type_="Valid", write_all_prop=self.writeAllProps)
            self.write_data_set(output_filename=outfile_TS, sample_size=self.testSize, type_="Test", write_all_prop=self.writeAllProps)
            if self.training_set!=trset:
                RaiseWarning(message="POLUS: Training set has been altered in stratified sampling")       
            if self.validation_set!=vlset:
                RaiseWarning(message="POLUS: Validation set has been altered in stratified sampling")       
            if self.test_set!=tsset:
                RaiseWarning(message="POLUS: Test set has been altered in stratified sampling")       
            atoms.append(atom)
        # Write job details for ferebus
        if not self.writeAllProps:
            props = [self.prop]
        WriteJobDetails(InDir=self.inputDir,OutDir=self.outputDir,system=self.systemName,atoms=atoms,props=props,overwrite=True)
     
