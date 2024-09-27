import os
import time
import random
import sys
import shutil
from typing import List, Union
import csv
from  polus.utils.read_module import readfile
from  polus.utils.logging import RaiseWarning

listOfLists = List[List[Union[float,int,str]]]
ID_sequence    = List[int]

class RS():  # RS: Random Sampling
    
    def __init__(self,filename, output_prop,*args):
        
        self.filename       = filename
        self.prop           = output_prop
        self.training_set   = None
        self.validation_set = None
        self.test_set = None
        self.pop_size = None
        for ar in args:
            if isinstance(ar,int):
                self.val_set_size = ar
            else:
                RaiseWarning(message=" Unexpected optional argument encountered")
                print("POLUS: Program complains::: Unexpected optional argument encountered!")   

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

    def get_training_point_IDs(self,size_train):
        
        # Get pop size
        if self.pop_size!=None:
            pop_size = self.pop_size
        else:
            ref_data_filename = self.filename
            file_content = readfile(ref_data_filename,self.prop)
            pop_size = len(file_content[1])
        
        # Get train set idx
        if self.training_set!=None:
            result = self.training_set
        else:
            indices = [x for x in range(pop_size)]
            for _ in range(10):
                random.shuffle(indices)
            result = indices[:size_train]
            self.training_set = result
            count_duplicates = 0
            for i in range(len(result)):
                index_i = result[i]
                for j in range(i+1,len(result)):
                    index_j = result[j]
                    if index_i==index_j:
                        count_duplicates+=1
            #print ("POLUS: Number of duplicate points:", count_duplicates)
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
        for _ in range(10):
            random.shuffle(indices)
        complementary = []
        if self.validation_set==None:
            for idx in indices:
                if idx not in training_set_IDs:
                    complementary.append(idx)
            result = complementary[:size_valid]
            self.validation_set = result
        else:
            result = self.validation_set
                
        
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
        
        if self.test_set==None:
            indices = [x for x in range(pop_size)]
            for _ in range(10):
                random.shuffle(indices)
            complementary = []
            for idx in indices:
                if idx not in training_set_IDs and idx not in validation_set_IDs:
                    complementary.append(idx)
            random.shuffle(complementary)
            result = complementary[:size_test]
            self.test_set = result
        else:
            result=self.test_set
        
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
            sys.exit("POLUS| Dataset type can only be 'T' (training) or 'V'(Validation)")
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

    def write_data_set(self, output_filename, sample_size, type_="Train", write_all_prop=False,current_prop=None):
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
                data_row[-1]=data_row[-1][:-2]
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

