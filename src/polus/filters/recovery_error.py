#!/usr/bin/env python
import os
import shutil
import sys
from polus.utils.read_module import readfile
from polus.files.inputs import GetListInputFiles
from polus.utils.logging import RaiseError

import time

class recovEnergy():
    def __init__(self,list_atoms,system_name,list_files=None,working_directory=None,target_prop=None, geom_IDs=None,input_directory=None):
        self.list_atoms             = list_atoms
        self.natoms                 = len(list_atoms)
        self.system_name            = system_name
        if target_prop==None:
            self.target_prop        = "wfn_energy"
        else:
            self.target_prop        = target_prop

        if working_directory==None:
            self.wdir               = os.getcwd()
        else:
            self.wdir               = working_directory
        self.input_dir              = input_directory 
        self.geom_IDs               = geom_IDs
        self.list_wfn_energies      = None
        self.atomic_IQA_energies    = None
        self.molecular_IQA_energies = None
        self.recov_energies         = None
        self.abs_recov_energies     = None
        self.list_files             = list_files


    def generate_list_files(self):
        print("POLUS: Input files reading from recovery_error module ...")
        list_input_dir = os.listdir(self.input_dir)
        if self.list_files==None:
            self.list_files = list()
            for atom in self.list_atoms:
                found = False
                for file_ in list_input_dir:
                    if atom in file_.split("_"):
                        self.list_files.append(os.path.join(self.input_dir,file_))
                        found = True
                if not found:
                    RaiseError(message=f"POLUS: Input file for atom {atom} not found")

        for file_ in self.list_files:
            print("POLUS:           {}".format(file_))

    def get_list_wfn_energies(self):
        if self.list_files==None:
            self.generate_list_files()
        input_data = readfile(self.list_files[0], self.target_prop)
        file_body, prop_idx = input_data[1], input_data[4]
        if self.geom_IDs==None:
            list_geom_IDs = [i for i in range(len(file_body))]
        else:
            list_geom_IDs = self.geom_IDs
        self.geom_IDs = list_geom_IDs

        if self.list_wfn_energies==None:
            self.list_wfn_energies = []
            for i in range(len(file_body)):
                if i in list_geom_IDs:
                    line = file_body[i]
                    self.list_wfn_energies.append(eval(line.split(",")[prop_idx]))

    
    def get_atomic_IQA_energies(self):
        if self.list_files==None:
            self.generate_list_files()
        size_dataset = len(readfile(self.list_files[0],prop="iqa")[1])
        if self.geom_IDs==None:
            list_geom_IDs = [i for i in range(size_dataset)]
        else:
            list_geom_IDs = self.geom_IDs

        if self.atomic_IQA_energies==None:
            self.atomic_IQA_energies = {}
            count_atoms = 0
            for i in range(self.natoms):
                atom = self.list_atoms[i]
                self.atomic_IQA_energies[atom] = []
                data = readfile(self.list_files[i])
                file_content, iqa_idx = data[1], data[4]
                for j in range(len(file_content)):
                    if j in list_geom_IDs:
                        line = file_content[j]
                        self.atomic_IQA_energies[atom].append(eval(line.split(",")[iqa_idx]))
                count_atoms +=1

        return self.atomic_IQA_energies

    def get_molecular_iqa_energies(self):
        "+++ Check for existence of atomic IQA energies +++"
        if self.atomic_IQA_energies==None:
            atomic_IQA_energies = self.get_atomic_IQA_energies()
        else:
            atomic_IQA_energies = self.atomic_IQA_energies
        
        "+++ Compute molecular IQA energies +++"
        if self.molecular_IQA_energies == None:
            self.molecular_IQA_energies = []
            for i in range(len(atomic_IQA_energies[self.list_atoms[0]])):
                mol_energy = 0
                for key in atomic_IQA_energies.keys():
                    mol_energy += atomic_IQA_energies[key][i]
                self.molecular_IQA_energies.append(mol_energy)

    def get_recovery_energies(self):
        self.get_list_wfn_energies()
        self.get_molecular_iqa_energies()
        recovery_energies = list()    # Defined as (wfn - iqa)
        abs_recovery_energies = list()
        if len(self.list_wfn_energies)==len(self.molecular_IQA_energies):
            for i in range(len(self.list_wfn_energies)):
                rec_energy = self.list_wfn_energies[i] - self.molecular_IQA_energies[i]
                recovery_energies.append(rec_energy)
                abs_recovery_energies.append(abs(rec_energy))
        else:
            sys.exit("Wfn energies and IQA energies not matching")

        self.recov_energies = recovery_energies
        self.abs_recov_energies = abs_recovery_energies
        
        return recovery_energies, abs_recovery_energies

    def write_recov_err_files(self,threshold=1.0,output_filename = None, dual_flag=None):
        if output_filename!=None:
            filename = output_filename
        else:
            filename = self.system_name.upper()+"_RECOV_ENERGIES"

        print("POLUS| Computing recovery energies ...")

        REC_ENERGY = self.get_recovery_energies()
        MOL_IQA    = self.molecular_IQA_energies
        MOL_WEn    = self.list_wfn_energies
        if REC_ENERGY == None:
            sys.exit("POLUS| Program complains::: Unable to execute the task due to missing WFN energies")
        recov_energies, abs_recov_energies = REC_ENERGY[0], REC_ENERGY[1]
        geom_IDs = self.geom_IDs
        with open(filename,"w") as rec_energy_file:
            rec_energy_file.write("ID \tIQA_MOL \tWFN_EN\tREC_ENERGY(Ha)\tABS_REC_ENERGY(Ha)\tREC_ENERGY(kJ/mol)\tABS_REC_ENERGY(kJ/mol)\n")
            for i in range(len(recov_energies)):
                if i!=len(recov_energies)-1:
                    rec_energy_file.write(str(geom_IDs[i])+"\t"+str(MOL_IQA[i])+"\t"+str(MOL_WEn[i])+"\t"+str(recov_energies[i])+"\t"+str(abs_recov_energies[i])\
                                          +"\t"+str(2625.5*recov_energies[i])+"\t"+str(2625.5*abs_recov_energies[i])+"\n")
                else:
                    rec_energy_file.write(str(geom_IDs[i])+"\t"+str(MOL_IQA[i])+"\t"+str(MOL_WEn[i])+"\t"+str(recov_energies[i])+"\t"+str(abs_recov_energies[i])\
                                          +"\t"+str(2625.5*recov_energies[i])+"\t"+str(2625.5*abs_recov_energies[i]))

        if (dual_flag!=None or dual_flag==False):
            bad_geoms_file  = self.system_name.upper()+"-BAD-GEOMETRIES-IQA"
            good_geoms_file = self.system_name.upper()+"-GOOD-GEOMETRIES-IQA"
            filtered_dir    = os.path.join(self.wdir,"FILTERED-BY-IQA")
        else:
            bad_geoms_file  = self.system_name.upper()+"-BAD-GEOMETRIES-DUAL"
            good_geoms_file = self.system_name.upper()+"-GOOD-GEOMETRIES-DUAL"
            filtered_dir    = os.path.join(self.wdir,"FILTERED-BY-DUAL")

        count_bad_geoms = 0
        count_good_geoms = 0
        good_geoms_IDs = []
        my_bad_geoms_file = open(bad_geoms_file,"w")
        my_good_geoms_file = open(good_geoms_file,"w")
        my_bad_geoms_file.write("Geom_ID\t\tRecovery_error(kJ/mol)\n")
        my_good_geoms_file.write("Geom_ID\t\tRecovery_error(kJ/mol)\n")
        for i in range(len(recov_energies)):
            if (abs(2625.5*recov_energies[i])>threshold):
                count_bad_geoms +=1
                my_bad_geoms_file.write(str(geom_IDs[i])+"\t\t"+str(2625.5*recov_energies[i])+"\n")   
            else:
                count_good_geoms +=1
                my_good_geoms_file.write(str(geom_IDs[i])+"\t\t"+str(2625.5*recov_energies[i])+"\n")   
                good_geoms_IDs.append(i)

        print("POLUS| Job completed successfully")
        print(f"POLUS| Total geometries               {len(geom_IDs)}")
        print(f"POLUS| Removed geometries             {count_bad_geoms}")
        print(f"POLUS| Retained geometries            {count_good_geoms}")
        print("POLUS| Writing filtered datasets ...") 
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
                #new_filename = os.path.join(filtered_dir,os.path.split(file)[1][:-4]+"_filtered_iqa.csv")
                new_filename = os.path.join(filtered_dir,os.path.split(file)[1])
                with open(new_filename,"w") as new_file:
                    new_file.write(header)
                    for i in range(len(body)):
                        if i in good_geoms_IDs:
                            new_file.write(body[i])
                print("POLUS|           {}".format(new_filename))
                time.sleep(0.1)

        

        
        


        
        
        
