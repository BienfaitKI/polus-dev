#!/usr/bin/env python
import os
import shutil
import sys
from utils.read_module import readfile
from files.inputs import get_list_input_files
import time

class recovQ00():
    def __init__(self,list_atoms,system_name,working_directory=None,target_prop=None, geom_IDs=None,input_directory=None):
        self.list_atoms             = list_atoms
        self.natoms                 = len(list_atoms)
        self.system_name            = system_name
        if target_prop==None:
            self.target_prop        = "q00"
        else:
            self.target_prop        = target_prop

        if working_directory==None:
            self.wdir               = os.getcwd()
        else:
            self.wdir               = working_directory
        self.input_dir              = input_directory 
        self.geom_IDs               = geom_IDs
        self.list_wfn_q00           = None
        self.atomic_q00             = None
        self.molecular_q00          = None
        self.recov_q00              = None
        self.abs_recov_q00          = None
        self.list_files             = None


    def generate_list_files(self):
        self.list_files = get_list_input_files(self.input_dir)
        "+++ Print list input files +++"
        print("POLUS| Input files ...")
        for file_ in self.list_files:
            print("POLUS|           {}".format(file_))

    def get_list_wfn_q00(self):
        if self.list_files==None:
            self.generate_list_files()
        input_data = readfile(self.list_files[0], self.target_prop)
        file_body, prop_idx = input_data[1], input_data[4]
        if self.geom_IDs==None:
            list_geom_IDs = [i for i in range(len(file_body))]
        else:
            list_geom_IDs = self.geom_IDs
        self.geom_IDs = list_geom_IDs
        if self.list_wfn_q00==None:
            self.list_wfn_q00 = list()
            for i in range(len(file_body)):
                if i in list_geom_IDs:
                    self.list_wfn_q00.append(0.0)
    
    def get_atomic_q00(self):
        if self.list_files==None:
            self.generate_list_files()
        size_dataset = len(readfile(self.list_files[0],prop=self.target_prop)[1])
        if self.geom_IDs==None:
            list_geom_IDs = [i for i in range(size_dataset)]
        else:
            list_geom_IDs = self.geom_IDs

        if self.atomic_q00==None:
            self.atomic_q00 = {}
            count_atoms = 0
            for i in range(self.natoms):
                atom = self.list_atoms[i]
                self.atomic_q00[atom] = list()
                data = readfile(self.list_files[i],prop=self.target_prop)
                file_content, q00_idx = data[1], data[4]
                for j in range(len(file_content)):
                    if j in list_geom_IDs:
                        line = file_content[j]
                        self.atomic_q00[atom].append(eval(line.split(",")[q00_idx]))
                count_atoms +=1

        return self.atomic_q00

    def get_molecular_q00(self):
        "+++ Check for existence of atomic IQA energies +++"
        if self.atomic_q00==None:
            atomic_q00 = self.get_atomic_q00()
        else:
            atomic_q00 = self.atomic_q00
        
        "+++ Compute molecular IQA energies +++"
        if self.molecular_q00 == None:
            self.molecular_q00 = list()
            for i in range(len(atomic_q00[self.list_atoms[0]])):
                mol_q00 = 0
                for key in atomic_q00.keys():
                    mol_q00 += atomic_q00[key][i]
                self.molecular_q00.append(mol_q00)

    def get_recovery_charge(self):
        self.get_list_wfn_q00()
        self.get_molecular_q00()
        recovery_q00 = list()    # Defined as (wfn - iqa)
        abs_recovery_q00 = list()
        if len(self.list_wfn_q00)==len(self.molecular_q00):
            for i in range(len(self.list_wfn_q00)):
                rec_q00 = self.list_wfn_q00[i] - self.molecular_q00[i]
                recovery_q00.append(rec_q00)
                abs_recovery_q00.append(abs(rec_q00))
        else:
            sys.exit("POLUS| Program complains :: WFN-Q00 and QCT-Q00 not macthing")

        self.recov_q00 = recovery_q00
        self.abs_recov_q00 = abs_recovery_q00
        
        return recovery_q00, abs_recovery_q00

    def write_recov_err_files(self,threshold=1.0,output_filename = None):
        if output_filename!=None:
            filename = output_filename
        else:
            filename = self.system_name.upper()+"-RECOVERY-Q00"

        print("POLUS| Computing recovery charges ...")

        REC_Q00 = self.get_recovery_charge()
        if REC_Q00 == None:
            sys.exit("POLUS| Program complains::: Unable to execute the task due to missing WFN energies")
        recov_q00, abs_recov_q00 = REC_Q00[0], REC_Q00[1]
        geom_IDs = self.geom_IDs
        with open(filename,"w") as rec_q00_file:
            rec_q00_file.write("ID,REC_Q00(e)\t\t ABS_REC_Q00(e)\n")
            for i in range(len(recov_q00)):
                if i!=len(recov_q00)-1:
                    rec_q00_file.write(str(geom_IDs[i])+"\t\t"+str(recov_q00[i])+"\t\t"+str(abs_recov_q00[i])+"\n")
                else:
                    rec_q00_file.write(str(geom_IDs[i])+"\t\t"+str(recov_q00[i])+"\t\t"+str(abs_recov_q00[i]))

        bad_geoms_file = self.system_name.upper()+"-BAD-GEOMETRIES-Q00"
        good_geoms_file = self.system_name.upper()+"-GOOD-GEOMETRIES-Q00"
        count_bad_geoms = 0
        count_good_geoms = 0
        good_geoms_IDs = []
        my_bad_geoms_file = open(bad_geoms_file,"w")
        my_good_geoms_file = open(good_geoms_file,"w")
        my_bad_geoms_file.write("Geom_ID\t\tRecovery_q00\n")
        my_good_geoms_file.write("Geom_ID\t\tRecovery_q00\n")
        for i in range(len(recov_q00)):
            if (abs(recov_q00[i])>threshold):
                count_bad_geoms +=1
                my_bad_geoms_file.write(str(geom_IDs[i])+"\t\t"+str(recov_q00[i])+"\n")   
            else:
                count_good_geoms +=1
                my_good_geoms_file.write(str(geom_IDs[i])+"\t\t"+str(recov_q00[i])+"\n")   
                good_geoms_IDs.append(i)

        print("POLUS| Job completed successfully")
        print(f"POLUS| Total geometries               {len(geom_IDs)}")
        print(f"POLUS| Removed geometries             {count_bad_geoms}")
        print(f"POLUS| Retained geometries            {count_good_geoms}")
        print("POLUS| Writing filtered datasets ...") 
        filtered_dir = os.path.join(self.wdir,"FILTERED-BY-Q00")
        if os.path.isdir(filtered_dir):
            shutil.rmtree(filtered_dir)
        os.mkdir(filtered_dir)
        print("POLUS| Filtered files ...")
        if (self.list_files!=None):
            for file in self.list_files:
                with open(file,"r") as old_file:
                    content = old_file.readlines()
                    header = content[0]
                    body = content[1:]
                #new_filename = os.path.join(filtered_dir,os.path.split(file)[1][:-4]+"_filtered_q00.csv")
                new_filename = os.path.join(filtered_dir,os.path.split(file)[1])
                with open(new_filename,"w") as new_file:
                    new_file.write(header)
                    for i in range(len(body)):
                        if i in good_geoms_IDs:
                            new_file.write(body[i])
                print("POLUS|           {}".format(new_filename))
                time.sleep(0.1)

        

        
        


        
        
        
