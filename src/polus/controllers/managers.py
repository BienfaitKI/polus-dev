#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 13:53:16 2022

@author: q63762bi
"""

import time
import sys
import os
import shutil
from samplers.RS.randomSampling import RS
from samplers.passiveSampling import PS
from samplers.uncertESS import UESS
from samplers.uncertESS2 import UESS_val_uess
from samplers.uncertESS3 import UESS_revised
from samplers.stratRS import SRS
from utils.printing import *
from utils.read_module import read_job_details, check_property
from filters.recovery_error import recovEnergy
from filters.iqa_correction import iqa_correct
from configparser import ConfigParser


def get_FInFile_IO_Dirs(input_dir,output_dir):
    First_input_filename = None
    """ Find first input file based on system_name and atom """ 
    job_details = read_job_details()
    system_name = job_details[0]
    atom        = job_details[5][0]

    cwd = os.getcwd()
    if input_dir==None or input_dir==os.getcwd():
        input_files_dir = os.getcwd()
    else:
        input_files_dir = os.path.join(cwd,input_dir)
    if output_dir==None or output_dir==os.getcwd():
        output_files_dir = os.getcwd()
    else:
        output_files_dir = os.path.join(cwd,output_dir)
    if input_files_dir!=None and not os.path.isdir(input_files_dir):
        sys.exit("Cannot find input directory")
    if output_files_dir!=os.getcwd() and not os.path.isdir(output_files_dir):
        try:
            os.mkdir(output_files_dir)
        except FileExistsError:
            print(f"{os.path.join(cwd,output_dir)} already exists") 

    OutDir = output_files_dir

    basename = None
    basename1 = system_name.upper()+"_"+atom.upper()+".csv"
    basename2 = atom.upper()+"_features_with_properties.csv"
    if basename1 in os.listdir(default_input_dir) or basename1 in os.listdir(input_files_dir):
        basename = basename1
    elif basename2 in os.listdir(default_input_dir) or basename2 in os.listdir(input_files_dir):
        basename = basename2
    else:
        print(" Cannot build basename. Job halted!")
        sys.exit()

    if basename in os.listdir(input_files_dir):
        First_input_filename = os.path.join(input_files_dir,basename)
        InDir = input_files_dir
    else:
        First_input_filename = os.path.join(default_input_dir,basename)
        InDir = default_input_dir

    return First_input_filename, InDir, OutDir

def get_mol_strat_details(input_dir, system_name, list_atoms,prop):
    most_informative_atom = list_atoms[0]
    input_file = os.path.join(input_dir,system_name+"_"+list_atoms[0]+".csv")
    if not os.path.isfile(input_file):
        print (f" Cannot find input file {input_file}")
        sys.exit()
    else:
        job = SRS(input_file,prop)
        largest_strat_regions = job.get_number_of_bins() 
    for atom in list_atoms[1:]:
        input_file = os.path.join(input_dir,system_name+"_"+atom+".csv")
        if not os.path.isfile(input_file):
            print (f" Cannot find input file {input_file}")
            sys.exit()
        else:
            job = SRS(input_file,prop)
            strat_regions = job.get_number_of_bins() 
            if (strat_regions > largest_strat_regions):
                most_informative_atom = atom
                largest_strat_regions = strat_regions

    return most_informative_atom, largest_strat_regions
          

def run_SRS(input_dir=None,output_dir=None,atoms=None,props=None,train=None,val=None,test=None):
    start = time.time()
    input_filename = None
    all_prop = False

    "+++ Get job details +++"
    job_details = read_job_details()
    system_name = job_details[0]
    strat_method = job_details[1]
    list_atoms = job_details[5]
    list_props = job_details[6]
    tr_set_size = job_details[2]
    ival_set_size = job_details[3]
    eval_set_size = job_details[4]
    if atoms!=None and isinstance(atoms,str):
        list_atoms = [entry for entry in atoms.split("_")]
    if props!=None and isinstance(props,str):
        list_props = [entry for entry in props.split("_")]
    if train!=None and isinstance(train,int):
        tr_set_size = train 
    if val!=None and isinstance(val,int):
        ival_set_size = train 
    if test!=None and isinstance(test,int):
        eval_set_size = train 
    
    """Check I/O directories"""
     
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,output_dir)

    """ Perform sampling"""
    for prop in list_props:
        most_inform_atom, number_of_bins = get_mol_strat_details(input_files_dir, system_name, list_atoms,prop)
        InFile = os.path.join(input_files_dir,system_name+"_"+most_inform_atom+".csv")
        print (f" |>>> Prop: {prop} ")
        print (f" |>>> Most informative atom: {most_inform_atom} ")
        print (f" |>>> Largest number of bins: {number_of_bins} ")
        try:
            print(f" SAR has launched the prob-job for molecular-wise stratified random sampling for {prop}")
            prob_job = SRS(InFile, prop)
            prob_job.get_validation_point_IDs(ival_set_size)
            prob_job.get_training_point_IDs(tr_set_size)
            prob_job.get_test_point_IDs(eval_set_size)
            print(" Outcome: successful")
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(" SAR cannot perform molecular-wise stratified random sampling with the prob-job")
            print(" The following error was caught:\n")
            print(exc_value)
            sys.exit()
    
        """ Save train, val and test sets IDs """
        train_set_IDs = prob_job.get_training_set()
        val_set_IDs   = prob_job.get_validation_set()
        test_set_IDs  = prob_job.get_test_set()

        """ Write out indices of Tr, Val and Test sets"""
        tr_file=open(os.path.join(output_files_dir,"INDICES-TRAIN.idx"),"w")
        for val in train_set_IDs:
            tr_file.write(str(val)+"\n")
        val_file=open(os.path.join(output_files_dir,"INDICES-VAL.idx"),"w")
        for val in val_set_IDs:
            val_file.write(str(val)+"\n")
        test_file=open(os.path.join(output_files_dir,"INDICES-TEST.idx"),"w")
        for val in test_set_IDs:
            test_file.write(str(val)+"\n")
        tr_file.close()
        val_file.close()
        test_file.close()

        """ Let us do the sampling now """
        for atom in list_atoms:
            if (system_name.upper()+"_"+atom+".csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,system_name.upper()+"_"+atom+".csv")
            elif (atom+"_features_with_properties.csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,atom+"_features_with_properties.csv")
            else:
                print(f" SAR cannot find input file for atom {atom}")
                sys.exit()

            print (f" SAR will now perform the task: SRS [{atom}]-{prop}") 
            outdir = os.path.join(output_files_dir,prop)
            tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
            ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
            eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
            if check_property(input_filename,prop):
                job = SRS(input_filename, prop)  
                job.set_training_set(train_set_IDs)
                job.set_validation_set(val_set_IDs)
                job.set_test_set(test_set_IDs)
                job.write_data_set(tr_set_filename,tr_set_size,"Train",all_prop)
                job.write_data_set(ival_set_filename,ival_set_size,"Valid",all_prop)
                job.write_data_set(eval_set_filename,eval_set_size,"Test",all_prop)
            else:
                print(f" >>> SRS [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()

def run_PS(input_dir=None,output_dir=None,atoms=None,props=None,train=None,val=None,test=None):
    start = time.time()
    input_filename = None
    
    all_prop = False
    
    " Get job details"
    job_details = read_job_details()
    system_name = job_details[0]
    strat_method = job_details[1]
    tr_set_size = job_details[2]
    ival_set_size = job_details[3]
    eval_set_size = job_details[4]
    list_atoms = job_details[5]
    list_props = job_details[6]
    PS_init_method = job_details[7]
    PS_uncertainty_method = job_details[8]
    if atoms!=None and isinstance(atoms,str):
        list_atoms = [entry for entry in atoms.split("_")]
    if props!=None and isinstance(props,str):
        list_props = [entry for entry in props.split("_")]
    if train!=None and isinstance(train,int):
        tr_set_size = train 
    if val!=None and isinstance(val,int):
        ival_set_size = train 
    if test!=None and isinstance(test,int):
        eval_set_size = train 
    
        
    """Check I/O directories"""
     
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,output_dir)

    """Run sampling mol-IDs for all atoms and properties"""
    if (len(list_props)>=1 and FInFile!=None):
        try:
            print(" SAR has launched the prob-job for molecular-wise passive sampling")
            prob_job = PS(FInFile, list_props[0],list_atoms[0],PS_init_method,PS_uncertainty_method)
            prob_job.get_validation_point_IDs(ival_set_size)
            prob_job.get_training_point_IDs(tr_set_size)
            prob_job.get_test_point_IDs(eval_set_size)
            print(" Outcome: successful")
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(" SAR cannot perform molecular-wise passive sampling with the prob-job")
            print(" The following error was caught:\n")
            print(exc_value)
            sys.exit()
    
    """ Save train, val and test sets IDs """
    train_set_IDs = prob_job.get_training_set()
    val_set_IDs   = prob_job.get_validation_set()
    test_set_IDs  = prob_job.get_test_set()
    uncertainty_data = prob_job.get_uncertainty_curve_copy()

    """ Write out indices of Tr, Val and Test sets"""
    tr_file=open(os.path.join(output_files_dir,"INDICES-TRAIN.idx"),"w")
    for val in train_set_IDs:
        tr_file.write(str(val)+"\n")
    val_file=open(os.path.join(output_files_dir,"INDICES-VAL.idx"),"w")
    for val in val_set_IDs:
        val_file.write(str(val)+"\n")
    test_file=open(os.path.join(output_files_dir,"INDICES-TEST.idx"),"w")
    for val in test_set_IDs:
        test_file.write(str(val)+"\n")
    tr_file.close()
    val_file.close()
    test_file.close()

    """ Let us do the sampling now """
    for atom in list_atoms:
        for  prop in list_props:
            if (system_name.upper()+"_"+atom+".csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,system_name.upper()+"_"+atom+".csv")
            elif (atom+"_features_with_properties.csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,atom+"_features_with_properties.csv")
            else:
                print(f" SAR cannot find input file for atom {atom}")
                sys.exit()

            print (f" SAR will now perform the task: PS [{atom}]-{prop}") 
            outdir = os.path.join(output_files_dir,prop)
            tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
            ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
            eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
            if check_property(input_filename,prop):
                job = PS(input_filename, prop, atom, PS_init_method, PS_uncertainty_method)  
                job.set_training_set(train_set_IDs)
                job.set_validation_set(val_set_IDs)
                job.set_test_set(test_set_IDs)
                job.set_uncertainty_curve_copy(uncertainty_data)
                job.write_data_set(tr_set_filename,tr_set_size,"Train",system_name,all_prop, prop)
                job.write_data_set(ival_set_filename,ival_set_size,"Valid",system_name,all_prop, prop)
                job.write_data_set(eval_set_filename,eval_set_size,"Test",system_name,all_prop, prop)
            else:
                print(f" >>> PS [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()

def run_UESS(input_dir=None,output_dir=None,atoms=None,props=None,train=None,val=None,test=None):
    start = time.time()
    input_filename = None
    
    all_prop = False
    
    " Get job details"
    job_details = read_job_details()
    system_name = job_details[0]
    strat_method = job_details[1]
    tr_set_size = job_details[2]
    ival_set_size = job_details[3]
    eval_set_size = job_details[4]
    list_atoms = job_details[5]
    list_props = job_details[6]
    PS_init_method = job_details[7]
    PS_uncertainty_method = job_details[8]
    if atoms!=None and isinstance(atoms,str):
        list_atoms = [entry for entry in atoms.split("_")]
    if props!=None and isinstance(props,str):
        list_props = [entry for entry in props.split("_")]
    if train!=None and isinstance(train,int):
        tr_set_size = train 
    if val!=None and isinstance(val,int):
        ival_set_size = train 
    if test!=None and isinstance(test,int):
        eval_set_size = train 
    
    
    """Check I/O directories"""
     
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,output_dir)

    """ Perform sampling"""
    for prop in list_props:
        most_inform_atom, number_of_bins = get_mol_strat_details(input_files_dir, system_name, list_atoms,prop)
        InFile = os.path.join(input_files_dir,system_name+"_"+most_inform_atom+".csv")
        print (f" |>>> Prop: {prop} ")
        print (f" |>>> Most informative atom: {most_inform_atom} ")
        print (f" |>>> Largest number of bins: {number_of_bins} ")
        try:
            print(f" SAR has launched the prob-job for molecular-wise stratified random sampling for {prop}")
            prob_job = UESS(InFile, prop, most_inform_atom, strat_method,PS_init_method, PS_uncertainty_method)  
            prob_job.get_validation_point_IDs(ival_set_size)
            prob_job.get_training_point_IDs(tr_set_size)
            prob_job.get_test_point_IDs(eval_set_size)
            print(" Outcome: successful")
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(" SAR cannot perform molecular-wise stratified random sampling with the prob-job")
            print(" The following error was caught:\n")
            print(exc_value)
            sys.exit()
    
        """ Save train, val and test sets IDs """
        train_set_IDs    = prob_job.get_training_set()
        val_set_IDs      = prob_job.get_validation_set()
        test_set_IDs     = prob_job.get_test_set()
        uncertainty_data = prob_job.get_local_uncertainty_curves()

        """ Write out indices of Tr, Val and Test sets"""
        tr_file=open(os.path.join(output_files_dir,"INDICES-TRAIN.idx"),"w")
        for val in train_set_IDs:
            tr_file.write(str(val)+"\n")
        val_file=open(os.path.join(output_files_dir,"INDICES-VAL.idx"),"w")
        for val in val_set_IDs:
            val_file.write(str(val)+"\n")
        test_file=open(os.path.join(output_files_dir,"INDICES-TEST.idx"),"w")
        for val in test_set_IDs:
            test_file.write(str(val)+"\n")
        tr_file.close()
        val_file.close()
        test_file.close()

        """ Let us do the sampling now """
        for atom in list_atoms:
            if (system_name.upper()+"_"+atom+".csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,system_name.upper()+"_"+atom+".csv")
            elif (atom+"_features_with_properties.csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,atom+"_features_with_properties.csv")
            else:
                print(f" SAR cannot find input file for atom {atom}")
                sys.exit()

            print (f" SAR will now perform the task: UESS [{atom}]-{prop}") 
            outdir = os.path.join(output_files_dir,prop)
            tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
            ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
            eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
            if check_property(input_filename,prop):
                job = UESS(input_filename, prop, atom, strat_method,PS_init_method, PS_uncertainty_method)  
                job.set_training_set(train_set_IDs)
                job.set_validation_set(val_set_IDs)
                job.set_test_set(test_set_IDs)
                job.set_local_uncertainty_curves(uncertainty_data)
                job.write_data_set(tr_set_filename,tr_set_size,"Train",system_name,all_prop,prop)
                job.write_data_set(ival_set_filename,ival_set_size,"Valid",system_name,all_prop,prop)
                job.write_data_set(eval_set_filename,eval_set_size,"Test",system_name,all_prop,prop)
            else:
                print(f" >>> UESS [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()

def run_UESS2(input_dir=None,output_dir=None,atoms=None,props=None,train=None,val=None,test=None):
    start = time.time()
    input_filename = None
    
    all_prop = False
    
    " Get job details"
    job_details = read_job_details()
    system_name = job_details[0]
    strat_method = job_details[1]
    tr_set_size = job_details[2]
    ival_set_size = job_details[3]
    eval_set_size = job_details[4]
    list_atoms = job_details[5]
    list_props = job_details[6]
    PS_init_method = job_details[7]
    PS_uncertainty_method = job_details[8]
    if atoms!=None and isinstance(atoms,str):
        list_atoms = [entry for entry in atoms.split("_")]
    if props!=None and isinstance(props,str):
        list_props = [entry for entry in props.split("_")]
    if train!=None and isinstance(train,int):
        tr_set_size = train 
    if val!=None and isinstance(val,int):
        ival_set_size = train 
    if test!=None and isinstance(test,int):
        eval_set_size = train 
    
    
    """Check I/O directories"""
     
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,output_dir)

    """ Perform sampling"""
    for prop in list_props:
        most_inform_atom, number_of_bins = get_mol_strat_details(input_files_dir, system_name, list_atoms,prop)
        InFile = os.path.join(input_files_dir,system_name+"_"+most_inform_atom+".csv")
        print (f" |>>> Prop: {prop} ")
        print (f" |>>> Most informative atom: {most_inform_atom} ")
        print (f" |>>> Largest number of bins: {number_of_bins} ")
        try:
            print(f" SAR has launched the prob-job for molecular-wise stratified random sampling for {prop}")
            prob_job = UESS_val_uess(InFile, prop, most_inform_atom, strat_method,PS_init_method, PS_uncertainty_method)  
            prob_job.get_training_point_IDs(tr_set_size)
            prob_job.get_validation_point_IDs(ival_set_size)
            prob_job.get_test_point_IDs(eval_set_size)
            print(" Outcome: successful")
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(" SAR cannot perform molecular-wise stratified random sampling with the prob-job")
            print(" The following error was caught:\n")
            print(exc_value)
            sys.exit()
    
        """ Save train, val and test sets IDs """
        train_set_IDs    = prob_job.get_training_set()
        val_set_IDs      = prob_job.get_validation_set()
        test_set_IDs     = prob_job.get_test_set()
        uncertainty_data = prob_job.get_local_uncertainty_curves()
        uncertainty_data2= prob_job.get_uncertainty_curve_val()

        """ Write out indices of Tr, Val and Test sets"""
        tr_file=open(os.path.join(output_files_dir,"INDICES-TRAIN.idx"),"w")
        for val in train_set_IDs:
            tr_file.write(str(val)+"\n")
        val_file=open(os.path.join(output_files_dir,"INDICES-VAL.idx"),"w")
        for val in val_set_IDs:
            val_file.write(str(val)+"\n")
        test_file=open(os.path.join(output_files_dir,"INDICES-TEST.idx"),"w")
        for val in test_set_IDs:
            test_file.write(str(val)+"\n")
        tr_file.close()
        val_file.close()
        test_file.close()

        """ Let us do the sampling now """
        for atom in list_atoms:
            if (system_name.upper()+"_"+atom+".csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,system_name.upper()+"_"+atom+".csv")
            elif (atom+"_features_with_properties.csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,atom+"_features_with_properties.csv")
            else:
                print(f" SAR cannot find input file for atom {atom}")
                sys.exit()

            print (f" SAR will now perform the task: UESS2 [{atom}]-{prop}") 
            outdir = os.path.join(output_files_dir,prop)
            tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
            ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
            eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
            if check_property(input_filename,prop):
                job = UESS_val_uess(input_filename, prop, atom, strat_method,PS_init_method, PS_uncertainty_method)  
                job.set_training_set(train_set_IDs)
                job.set_validation_set(val_set_IDs)
                job.set_test_set(test_set_IDs)
                job.set_local_uncertainty_curves(uncertainty_data)
                job.set_uncertainty_curve_val(uncertainty_data2)
                job.write_data_set(tr_set_filename,tr_set_size,"Train",system_name,all_prop,prop)
                job.write_data_set(ival_set_filename,ival_set_size,"Valid",system_name,all_prop,prop)
                job.write_data_set(eval_set_filename,eval_set_size,"Test",system_name,all_prop,prop)
            else:
                print(f" >>> UESS2 [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()

def run_UESS3(input_dir=None,output_dir=None,atoms=None,props=None,train=None,val=None,test=None):
    start = time.time()
    input_filename = None
    
    all_prop = False
    
    " Get job details"
    job_details = read_job_details()
    system_name = job_details[0]
    strat_method = job_details[1]
    tr_set_size = job_details[2]
    ival_set_size = job_details[3]
    eval_set_size = job_details[4]
    list_atoms = job_details[5]
    list_props = job_details[6]
    PS_init_method = job_details[7]
    PS_uncertainty_method = job_details[8]
    if atoms!=None and isinstance(atoms,str):
        list_atoms = [entry for entry in atoms.split("_")]
    if props!=None and isinstance(props,str):
        list_props = [entry for entry in props.split("_")]
    if train!=None and isinstance(train,int):
        tr_set_size = train 
    if val!=None and isinstance(val,int):
        ival_set_size = train 
    if test!=None and isinstance(test,int):
        eval_set_size = train 
    
    
    """Check I/O directories"""
     
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,output_dir)

    """ Perform sampling"""
    for prop in list_props:
        most_inform_atom, number_of_bins = get_mol_strat_details(input_files_dir, system_name, list_atoms,prop)
        InFile = os.path.join(input_files_dir,system_name+"_"+most_inform_atom+".csv")
        print (f" |>>> Prop: {prop} ")
        print (f" |>>> Most informative atom: {most_inform_atom} ")
        print (f" |>>> Largest number of bins: {number_of_bins} ")
        #try:
        #    print(f" SAR has launched the prob-job for molecular-wise stratified random sampling for {prop}")
        prob_job = UESS_revised(InFile, prop, most_inform_atom, strat_method,PS_init_method, PS_uncertainty_method)  
        prob_job.get_training_point_IDs(tr_set_size)
        prob_job.get_validation_point_IDs(ival_set_size)
        prob_job.get_test_point_IDs(eval_set_size)
        print(" Outcome: successful")
        #except:
        #    exc_type, exc_value, exc_traceback = sys.exc_info()
        #    print(" SAR cannot perform molecular-wise stratified random sampling with the prob-job")
        #    print(" The following error was caught:\n")
        #    print(exc_value)
         #   sys.exit()
    
        """ Save train, val and test sets IDs """
        train_set_IDs    = prob_job.get_training_set()
        val_set_IDs      = prob_job.get_validation_set()
        test_set_IDs     = prob_job.get_test_set()
        uncertainty_data = prob_job.get_uncertainty_curve_train()
        uncertainty_data2= prob_job.get_uncertainty_curve_val()

        """ Write out indices of Tr, Val and Test sets"""
        tr_file=open(os.path.join(output_files_dir,"INDICES-TRAIN.idx"),"w")
        for val in train_set_IDs:
            tr_file.write(str(val)+"\n")
        val_file=open(os.path.join(output_files_dir,"INDICES-VAL.idx"),"w")
        for val in val_set_IDs:
            val_file.write(str(val)+"\n")
        test_file=open(os.path.join(output_files_dir,"INDICES-TEST.idx"),"w")
        for val in test_set_IDs:
            test_file.write(str(val)+"\n")
        tr_file.close()
        val_file.close()
        test_file.close()

        """ Let us do the sampling now """
        for atom in list_atoms:
            if (system_name.upper()+"_"+atom+".csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,system_name.upper()+"_"+atom+".csv")
            elif (atom+"_features_with_properties.csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,atom+"_features_with_properties.csv")
            else:
                print(f" SAR cannot find input file for atom {atom}")
                sys.exit()

            print (f" SAR will now perform the task: UESS3 [{atom}]-{prop}") 
            outdir = os.path.join(output_files_dir,prop)
            tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
            ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
            eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
            if check_property(input_filename,prop):
                job = UESS_revised(input_filename, prop, atom, strat_method,PS_init_method, PS_uncertainty_method)  
                job.set_training_set(train_set_IDs)
                job.set_validation_set(val_set_IDs)
                job.set_test_set(test_set_IDs)
                job.set_uncertainty_curve_train(uncertainty_data)
                job.set_uncertainty_curve_val(uncertainty_data2)
                job.write_data_set(tr_set_filename,tr_set_size,"Train",system_name,all_prop,prop)
                job.write_data_set(ival_set_filename,ival_set_size,"Valid",system_name,all_prop,prop)
                job.write_data_set(eval_set_filename,eval_set_size,"Test",system_name,all_prop,prop)
            else:
                print(f" >>> UESS3 [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()

def run_RS(input_dir=None,output_dir=None,atoms=None,props=None,train=None,val=None,test=None):
    start = time.time()
    input_filename = None
    
    all_prop = False
    
    " Get job details"
    job_details = read_job_details()
    system_name = job_details[0]
    strat_method = job_details[1]
    tr_set_size = job_details[2]
    ival_set_size = job_details[3]
    eval_set_size = job_details[4]
    list_atoms = job_details[5]
    list_props = job_details[6]
    PS_init_method = job_details[7]
    PS_uncertainty_method = job_details[8]
    if atoms!=None and isinstance(atoms,str):
        list_atoms = [entry for entry in atoms.split("_")]
    if props!=None and isinstance(props,str):
        list_props = [entry for entry in props.split("_")]
    if train!=None and isinstance(train,int):
        tr_set_size = train 
    if val!=None and isinstance(val,int):
        ival_set_size = train 
    if test!=None and isinstance(test,int):
        eval_set_size = train 
    
    """Check I/O directories"""
     
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,output_dir)

    """Run sampling mol-IDs for all atoms and properties"""
    if (len(list_props)>=1 and FInFile!=None):
        try:
            print(" SAR has launched the prob-job for molecular-wise random sampling")
            prob_job = RS(FInFile, list_props[0])  # does not matter which prop
            prob_job.get_training_point_IDs(tr_set_size)
            prob_job.get_validation_point_IDs(ival_set_size)
            prob_job.get_test_point_IDs(eval_set_size)
            print(" Outcome: successful")
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(" SAR cannot perform molecular-wise random sampling with the prob-job")
            print(" The following error was caught:\n")
            print(exc_value)
            sys.exit()

    """ Save train, val and test sets IDs """
    train_set_IDs = prob_job.get_training_set()
    val_set_IDs   = prob_job.get_validation_set()
    test_set_IDs  = prob_job.get_test_set()

    """ Write out indices of Tr, Val and Test sets"""
    tr_file=open(os.path.join(output_files_dir,"INDICES-TRAIN.idx"),"w")
    for val in train_set_IDs:
        tr_file.write(str(val)+"\n")
    val_file=open(os.path.join(output_files_dir,"INDICES-VAL.idx"),"w")
    for val in val_set_IDs:
        val_file.write(str(val)+"\n")
    test_file=open(os.path.join(output_files_dir,"INDICES-TEST.idx"),"w")
    for val in test_set_IDs:
        test_file.write(str(val)+"\n")
    tr_file.close()
    val_file.close()
    test_file.close()
    """ Let us do the sampling now """

    for atom in list_atoms:
        for  prop in list_props:
            if (system_name.upper()+"_"+atom+".csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,system_name.upper()+"_"+atom+".csv")
            elif (atom+"_features_with_properties.csv" in os.listdir(input_files_dir)):
                input_filename = os.path.join(input_files_dir,atom+"_features_with_properties.csv")
            else:
                print(f" SAR cannot find input file for atom {atom}")
                sys.exit()

            print (f" SAR will now perform the task: RS [{atom}]-{prop}") 
            outdir = os.path.join(output_files_dir,prop)
            tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
            ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
            eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
            if check_property(input_filename,prop):
                job = RS(input_filename, prop)  # does not matter which prop
                job.set_training_set(train_set_IDs)
                job.set_validation_set(val_set_IDs)
                job.set_test_set(test_set_IDs)
                job.write_data_set(tr_set_filename,tr_set_size,"Train",all_prop,prop)
                job.write_data_set(ival_set_filename,ival_set_size,"Valid",all_prop,prop)
                job.write_data_set(eval_set_filename,eval_set_size,"Test",all_prop,prop)
            #print (f"{atom}-{prop} samples successfully generated!")
            else:
                print(f" >>> RS [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()


def run_recov_calcs():

    job_details = read_job_details()
    recov_filter = job_details[9]
        
    " Get job details"
    job_details = read_job_details()
    system_name =job_details[0]
    list_atoms = job_details[5]
    
    job =  recovEnergy(list_atoms,system_name)
    job.write_recov_err_files(threshold=recov_filter)
    
    print_goodbye_message()

def run_iqa_corr(outdir,input_dir):
    print_welcome_message()
        
    " Get job details"
    job_details = read_job_details()
    system_name = job_details[0]
    list_atoms = job_details[5]

    job= iqa_correct(list_atoms,system_name)
    print (" This job will perform two tasks:")
    print (" Task 1: Generate and write out the raw and corrected IQA energies")
    print (" Task 2: Generate and write out the corrected reference IQA energies")
    print ("")
    print (" Performing Task 1 ...")
    job.write_raw_and_corrected_atomic_iqa_energies(output_dir=outdir)
    print (" Performing Task 2 ...")
    job.write_corrected_reference_data(write_all_prop=False,output_dir=outdir)
    wfn_energy=job.wfn_energies
    iqa_energy=job.molecular_iqa_energies
    min_wfn = min(wfn_energy)
    min_iqa = min(iqa_energy)
    """Check I/O directories"""
    FInFile, input_files_dir, output_files_dir = get_FInFile_IO_Dirs(input_dir,outdir)
    wfn_file = open(os.path.join(output_files_dir,"WFN-ENERGY"), "w")
    iqa_file = open(os.path.join(output_files_dir,"IQA-ENERGY"), "w")
    count_below = 0
    count_beyond = 0
    count_between = 0
    for val in wfn_energy:
        val = 2625.5*(val-min_wfn)
        #val = (val-min_wfn)
        wfn_file.write(str(val)+"\n")
        if val <35.5:
            count_below+=1
        if val >39.0:
            count_beyond+=1
        if 35.5<=val<=39.0:
            count_between+=1
    wfn_file.close() 
    for val in iqa_energy:
        val = 2625.5*(val-min_iqa)
        #val = (val-min_iqa)
        iqa_file.write(str(val)+"\n")
    wfn_file.close() 
    iqa_file.close() 
    print (f"BELOW  35.5 kJ/mol {100*count_below/len(wfn_energy)}")
    print (f"BEYOND 39.0 kJ/mol {100*count_beyond/len(wfn_energy)}")
    print (f"BETWEEN 39.0 kJ/mol {100*count_between/len(wfn_energy)}")
    print_goodbye_message()

             
def create_output_directories(list_props):
    cwd = os.getcwd()
    for prop in list_props:
        prop_dir = os.path.join(cwd,prop)
        if not os.path.isdir(prop_dir):
            os.mkdir(prop_dir)
            #print (f"{prop_dir} successfully created")
        else:
            shutil.rmtree(prop_dir)
            os.mkdir(prop_dir)
            #print (f"{prop_dir} removed and recreated!")
            

            
