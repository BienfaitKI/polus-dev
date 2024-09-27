#!/usr/bin/env python
import os
import sys
import numpy as np
from polus.config.user_inputs import read_config
from polus.config.user_inputs import read_atom_names
from polus.utils.logging import RaiseError, RaiseWarning
from polus.utils.read_module import readfile

def write_job_details(input_directory,root):
    config_details = read_config()
    system_name    = config_details[0]
    atoms          = config_details[1]
    props          = config_details[2]
    files          = os.listdir(input_directory)
    natoms         = len(files)

    if (atoms==None):
        atoms = read_atom_names(input_directory)

    "+++ Print list of input files +++"
    print("POLUS| Input files ...")
    for file_ in files:
        print("POLUS|           {}".format(os.path.join(input_directory,file_)))

    if ("alf" not in files[0]):
        status = 1
        print(" POLUS cannot write job-details file for FEREBUS")
    else:
        jd_file    = os.path.join(root,"job-details") # TODO: should be .ini file
        if not  os.path.isfile(jd_file):
            with open(jd_file,"w") as jd:
                jd.write("[system]\n")
                jd.write("system_name {} \n".format(system_name))            
                jd.write("natoms {} \n".format(natoms))            
                jd.write("atoms ")            
                for i in range(len(atoms)):
                    atom = atoms[i]
                    if (i!=natoms-1):
                        jd.write("{} ".format(atom))            
                    else:
                        jd.write("{} \n".format(atom))            
                jd.write("props ")            
                for j in range(len(props)):
                    prop = props[j]
                    if (j!=len(props)-1):
                        jd.write("{} ".format(prop))            
                    else:
                        jd.write("{} \n\n".format(prop))            
                jd.write("[ALF]\n")            
                has_zero_idx = check_alf(files)
                for file_ in files:
                    atom = file_.split(".")[0].split("_")[0]
                    if has_zero_idx:
                        alf  = [str(int(x)+1) for x in file_.split(".")[0].split("_alf_")[1].split("_")]
                    else:
                        alf  = [str(int(x)) for x in file_.split(".")[0].split("_alf_")[1].split("_")]
                    jd.write("{} {} {} {} \n".format(atom,alf[0],alf[1],alf[2]))
        status = 0
    return status


def check_alf(files):
    files_ = list()
    for file_ in files:
        files_.append(os.path.split(file_)[1])
    has_zero_idx = False
    for file_ in files_:
        alf  = [int(x)+1 for x in file_.split(".")[0].split("_alf_")[1].split("_")]
        if 0 in alf:
            has_zero_idx = True
            break

    return has_zero_idx
    
        
 
def WriteJobDetails(InDir=None,OutDir=None,system=None,atoms=None,props=None,overwrite=True):
    if (InDir==None):
        RaiseError(message="Program cannot find input directory")
    if (OutDir==None):
        RaiseError(message="Program cannot find output directory")
    if (system==None):
        RaiseError(message="Invalid system name")
    if (atoms==None):
        RaiseError(message="Invalid atom labels")
    if (props==None):
        RaiseError(message="Invalid target properties")
    files          = os.listdir(InDir)
    natoms         = len(atoms)

    if ("alf" not in files[0]):
        RaiseWarning(message=" Program cannot write job-details file for FEREBUS")
    else:
        jd_file    = os.path.join(OutDir,"job-details") # TODO: should be .ini file
        if overwrite or not os.path.isfile(jd_file):
            with open(jd_file,"w") as jd:
                jd.write("[system]\n")
                jd.write("system_name {} \n".format(system))            
                jd.write("natoms {} \n".format(natoms))            
                jd.write("atoms ")            
                for i in range(len(atoms)):
                    atom = atoms[i]
                    if (i!=natoms-1):
                        jd.write("{} ".format(atom))            
                    else:
                        jd.write("{} \n".format(atom))            
                jd.write("props ")            
                for j in range(len(props)):
                    prop = props[j]
                    if (j!=len(props)-1):
                        jd.write("{} ".format(prop))            
                    else:
                        jd.write("{} \n\n".format(prop))            
                jd.write("[ALF]\n")            
                has_zero_idx = check_alf(files)
                for file_ in files:
                    atom = file_.split(".")[0].split("_")[0]
                    try:
                        if has_zero_idx:
                            alf  = [str(int(x)+1) for x in file_.split(".")[0].split("_alf_")[1].split("_")]
                        else:
                            alf  = [str(int(x)) for x in file_.split(".")[0].split("_alf_")[1].split("_")]
                    except:
                        alf = None
                    if alf  == None:
                        alf  = ["1","2","3"]
                    jd.write("{} {} {} {} \n".format(atom,alf[0],alf[1],alf[2]))

                # Write prop descriptions
                jd.write("\n[PROPERTIES]\n")     
                jd.write(f"{'key':<10} {'Min':>12} {'Max':>12} {'Range':>12} {'Mean':>12} {'Median':>12} {'Std':>12} {'CV':>12}\n")
                for prop in props:
                    data = getPropDescription(atoms,prop,InDir)
                    for key in data.keys():
                        task = prop+"-"+key
                        cv   = data[key][5]/abs(data[key][3])
                        jd.write(f"{task:<10} {data[key][0]:>12.6f} {data[key][1]:>12.6f} {data[key][2]:>12.6f} {data[key][3]:>12.6f} {data[key][4]:>12.6f} {data[key][5]:>12.6f} {cv:>12.6f}\n")

def check_alf(files):
    has_zero_idx = False
    for file_ in files:
        alf  = [int(x) for x in file_.split(".")[0].split("_alf_")[1].split("_")]
        if 0 in alf:
            has_zero_idx = True
            break

    return has_zero_idx
    
        
def getPropDescription(atoms,prop,input_directory):
    data  = dict()
    files = list()
    items = os.listdir(input_directory)
    for atom in atoms:
        filename = None
        for item in items:
            if item.endswith(".csv") and atom in item.split("_"):
                 filename = os.path.join(input_directory,item)
        prop_vector = readfile(filename,prop) [2]
        data[atom]  = list()
        data[atom].append(min(prop_vector)) 
        data[atom].append(max(prop_vector)) 
        data[atom].append(max(prop_vector)-min(prop_vector)) 
        data[atom].append(np.mean(np.array(prop_vector)))
        data[atom].append(np.median(np.array(prop_vector)))
        data[atom].append(np.std(np.array(prop_vector)))
    return data
