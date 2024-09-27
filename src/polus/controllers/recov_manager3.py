#!/usr/bin/env python
import os
import sys
from config.user_inputs import read_config, read_cmd_args
from utils.io import get_FInFile_IO_Dirs
from filters.recovery_error import recovEnergy
from filters.recovery_error2 import recovQ00
from utils.printing import print_goodbye_message

global root
root = os.getcwd()

def run_recov_dual_calcs():
    "+++ Gather job details +++"
    cmd_args          = read_cmd_args()
    config_           = read_config()
    system_name       = config_[0]
    input_dir         = cmd_args[0]
    output_dir        = cmd_args[1]
    iqa_filter        = config_[7]
    q00_filter        = config_[8]
    list_atoms        = cmd_args[2]
    query             = get_FInFile_IO_Dirs(root,input_dir,output_dir)
    FInFile           = query[0]
    input_files_dir   = query[1]
    output_files_dir  = query[2]
    if (list_atoms==None):
        list_atoms    = config_[1]
        if (list_atoms==None):
            list_atoms    = query[3]
            if (list_atoms == None):
                sys.exit("POLUS| Program complains::: Cannot determine list of atoms")
    print("POLUS| IQA-FILTER threshold  {}".format(iqa_filter))    
    print("POLUS| Q00-FILTER threshold  {}".format(q00_filter))    
    "+++ Perform RECOV-Q00 task +++" 
    job =  recovQ00(list_atoms=list_atoms, \
                       system_name=system_name, \
                       working_directory=root, \
                       target_prop="q00", \
                       input_directory=input_files_dir)

    job.write_recov_err_files(threshold=q00_filter, \
                              output_filename = system_name.upper()+"-RECOVERY-Q00")
    
    "+++ Perform RECOV-IQA task +++" 
    job2 =  recovEnergy(list_atoms=list_atoms, \
                       system_name=system_name, \
                       working_directory=root, \
                       target_prop="wfn_energy", \
                       input_directory=input_files_dir)

    job2.write_recov_err_files(threshold=iqa_filter, \
                              output_filename = system_name.upper()+"-RECOVERY-IQA")

    "+++ Perform RECOV-IQA task +++" 
    input_dir_dual = os.path.join(root,"FILTERED-BY-Q00")
    print("POLUS| Performing recov task via DUAL-FILTER")
    print("POLUS| Input directory {}".format(input_dir_dual))
    job3 =  recovEnergy(list_atoms=list_atoms, \
                       system_name=system_name, \
                       working_directory=root, \
                       target_prop="wfn_energy", \
                       input_directory=input_dir_dual)

    job3.write_recov_err_files(threshold=iqa_filter, \
                              output_filename = system_name.upper()+"-RECOVERY-DUAL", \
                              dual_flag=True)
    print_goodbye_message()

