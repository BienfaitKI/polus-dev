# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 18:24:26 2023

@author: Bienfait
"""
import time
import os
from polus.utils.logging import RaiseError
from polus.utils.read_module import check_input_files, read_job_details
from datetime import datetime

def PrintJobDetails():
    now       = datetime.now()
    date      = now. strftime("%A %d %B %Y")
    time      = now. strftime("%H:%M:%S")
    hostname  = str(os.popen("hostname").read()).strip().strip("\n")
    osName    = os.name
    try:
        user      = str(os.getlogin())
    except:
        user      = "UNKNOWN"
    pid       = str(os.getpid())
    #print("GENERAL JOB DETAILS")
    #print("*"*82)
    print(f"POLUS: {'Job started on':<30} {date:>44}")
    print(f"POLUS: {'Job started at':<30} {time:>44}")
    print(f"POLUS: {'Job running on':<30} {hostname:>44}")
    #print(f"POLUS: {'Operating System':<30} {osName:>44}")
    #print(f"POLUS: {'UserName':<30} {user:>44}")
    #print(f"POLUS: {'Process ID':<30} {pid:>44}")
    #print("*"*82)
    #print("")
    #print("PROGRESS REPORT")
    #print("*"*82)
    #print(f"{'       Task'} {'Duration(s)':>70}")

def PrintOnTerminal(msg=None,duration=None,msgLength=None):
    if isinstance(msg,str):
        print(f"POLUS: {msg}",end="")
    elif isinstance(duration,float) and isinstance(msgLength,int):
        fieldSize = 75 - msgLength
        print(f"{duration:>{fieldSize}.6f}")
    else:
        RaiseError(message=" Program cannot print information on standard output terminal ")

def print_welcome_message():
    print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print ("     @@@@@@@@@@@@        &        Welcome  to  POLUS        &        @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@       &&&       Welcome  to  POLUS       &&&       @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@      &&&&&      Welcome  to  POLUS      &&&&&      @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@     &&&&&&&     Welcome  to  POLUS     &&&&&&&     @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@    &&&&&&&&&    Welcome  to  POLUS    &&&&&&&&&    @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@     &&&&&&&     Welcome  to  POLUS     &&&&&&&     @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@      &&&&&      Welcome  to  POLUS      &&&&&      @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@       &&&       Welcome  to  POLUS       &&&       @@@@@@@@@@@@     ")
    print ("     @@@@@@@@@@@@        &        Welcome  to  POLUS        &        @@@@@@@@@@@@     ")
    print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print (" ")
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(" Date and time =", dt_string)
    print (" ")
    time.sleep(2)

def PrintGBMessage():
    print ("")
    print ("THANKS FOR USING POLUS")
    print ("*"*81)
    print (" For more details: bienfait.isamura@postgrad.manchester.ac.uk")
    print ("                   paul.popelier@manchester.ac.uk")
    print ("*"*81)
    print (" ")

def print_goodbye_message():
    print ("")
    print ("*"*81)
    print ("GOODBYE")
    print ("*"*81)
    print (" Thanks for using POLUS")
    print (" For more details: bienfait.isamura@postgrad.manchester.ac.uk")
    print ("*"*81)
    print (" ")
def print_help():
    print (" ")
    print ("    This code performs various data sampling tasks (from any reference data) into training and validation sets.")
    print ("    Available sampling techniques include:")
    print ("        1. Random sampling (RS)")
    print ("        2. Stratified Random sampling (SRS)")
    print ("        3. Passive sampling (PS)")
    print ("        4. Uncertainty-enhanced stratified sampling (UESS = [ SRS + PS])")
    print (" ")
    print ("    The package depends on two python libraries: statistics and numerical python (numpy).")
    print ("    Make sure they are installed before running the main.py script with the appropriate job type argument.")
    print ("    ===> Assuming PIP is already available on your system, you can install the above libraries by typing the ")
    print ("         following commands on your terminal:")
    print ("         Numpy: pip install numpy")
    print ("         Statistics: pip install statistics")
    print (" ")
    print ("    The stratified random sampling (SRS) algorithm involves two main steps:")
    print ("        1. Stratification of the dataset")
    print ("        2. Random sampling within each stratification bin or region")
    print (" ")
    print ("    The stratification process is guided by the distribution of the target property.")
    print ("    The following stratification rules are implemented:")
    print ("        1. Sturges' rule")
    print ("        2. Scott's rule")
    print ("        3. Rice's rule")
    print ("        4. Doane's rule")
    print ("        5. Freedman-Diaconis rule (Default)")
    print ("        6. Equiprobable bins")
    print ("")
    print ("    Note that the user is also allowed to manually specify the desired number of stratification regions.")
    print ("    However, this is not recommended.")
    print ("")
    print ("    The following input files are required:")
    print ("        1. JOB_DETAILS.txt: contains the details of the sampling")
    print ("        2. Reference data files [system_atom.csv]: contain the raw data for sampling")
    print ("")
    print ("    Usage options")
    print ("    -------------")
    print ("    -h, --help                      print this help message")
    print ("    -v, --version                   print the version of the package")
    print ("    -c, --check                     check availability of input files")
    print ("    -C, --clean                     remove existing output folders")
    print ("    -r, --run      [job_type]       run job_type")
    print ("")
    print ("     [job types]: RS, SRS, PS, UESS")
    print (" ")
    print ("    For more details: bienfait.isamura@postgrad.manchester.ac.uk")
    print ("")

def print_version():
    print ("    SRS version 1.0")
    
def print_input_details():
    list_files = check_input_files()
    print ("    File                                   Status")
    print ("    ----                                   ------")
    if not list_files[0]:
        print ("    Job file                                 NOT FOUND")
        print ("")
        print ("    Nota: SRS is unable to check data files")
    else:
        print ("    Job file                               FOUND")
        job_details = read_job_details()
        atoms = job_details[5]
        for i in range(1,len(list_files)):
            if list_files[i]:
                print (f"    {atoms[i-1]}_file                                FOUND")
            else:
                print (f"    {atoms[i-1]}_file                                NOT FOUND")

def print_job_details(system_name,strat_method,tr_set_size,ival_set_size,\
                      eval_set_size,list_atoms,list_props):
    
    print("JOB DETAILS")
    print("-----------")
    print("")
    print("Parameter               Value")
    print("---------               -----")
    print(f"System name             {system_name}")
    print(f"Strat method            {strat_method}")
    print(f"Tr_set_size             {tr_set_size}")
    print(f"Int.val_set_size        {ival_set_size}")
    print(f"Ext.val_set_size        {eval_set_size}")
    print(f"Atoms                   {list_atoms[0]}  ",end="")
    for i in range(1,len(list_atoms)):
        if i!=len(list_atoms)-1:
            print(f"{list_atoms[i]}  ",end="")
        else:
            print(f"{list_atoms[i]}  ",end="\n")
    print(f"Properties              {list_props[0]}  ",end="")
    for i in range(1,len(list_props)):
        print(f"{list_props[i]}  ",end="")
    print("\n")
