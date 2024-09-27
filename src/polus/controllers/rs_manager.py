import time
import sys
import os
import shutil
import random
from polus.samplers.randomSampling import RS
from polus.utils.printing import *
from polus.utils.read_module import check_property
from polus.utils.io import get_FInFile_IO_Dirs
from polus.config.user_inputs import read_config, read_cmd_args
from polus.files.outputs import write_job_details

global root
root = os.getcwd()

def run_RS(input_dir=None,output_dir=None,atoms=None,props=None,trainSize=None,valSize=None,testSize=None):

    "+++ Set defaults +++"
    start = time.time()
    input_filename = None 
    all_prop = False

    "+++ Gather job details +++"
    cmd_args      = read_cmd_args()
    config_       = read_config()
    system_name   = config_[0]
    #input_dir     = cmd_args[0]
    #output_dir    = cmd_args[1]
    list_atoms    = cmd_args[2]
    list_props    = cmd_args[3]
    train         = cmd_args[6]
    val           = cmd_args[7]
    test          = cmd_args[8]

    "+++ Set user inputs +++"
    if (atoms !=None):
        list_atoms = atoms
    if (props !=None):
        list_props = props
    if (trainSize !=None):
        train = trainSize
    if (valSize !=None):
        val   = valSize
    if (testSize !=None):
        test = testSize

    "+++ Double check inputs +++"
    if (output_dir==None):
        input_dir     = root
    if (output_dir==None):
        output_dir    = root
    if (list_atoms==None):
        list_atoms    = config_[1]
    if (list_props==None):
        list_props    = config_[2]
    if (train==None):
        train         = config_[4]
    if (val==None):
        val           = config_[5]
    if (test==None):
        test          = config_[6]
    if (system_name==None):
        system_name = cmd_args[10]
        if (system_name==None):
            sys.exit("POLUS| program aborted due to missing system's names")
    
    "+++ Set I/O directories +++"
    query             = get_FInFile_IO_Dirs(root,input_dir,output_dir)
    FInFile           = query[0]
    input_files_dir   = query[1]
    output_files_dir  = query[2]
    if (list_atoms==None):
        list_atoms    = query[3]

    "+++ Write JD file for FEREBUS users +++!"
    status = write_job_details(input_files_dir,root)
    if (status==0):
        print("POLUS| Job-details file ...")
        print("POLUS|           {}".format(os.path.join(input_files_dir,"job-details")))
    else:
        print("POLUS! Program cannot write job-details file")

    "+++ Mol-based probing +++!"
    lg_train_set = max(train)
    lg_val_set   = max(val)
    lg_test_set  = max(test)
    try:
        prob_job = RS(FInFile, list_props[0])  # does not matter which prop
        prob_job.get_training_point_IDs(lg_train_set)
        prob_job.get_validation_point_IDs(lg_val_set)
        prob_job.get_test_point_IDs(lg_test_set)
        BIG_TRAIN = prob_job.get_training_set()
        BIG_VAL   = prob_job.get_validation_set()
        BIG_TEST  = prob_job.get_test_set()
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        print("POLUS| Program cannot perform molecular-wise random sampling with the prob-job")
        print("POLUS| The following error was caught:\n")
        print(exc_value)
        sys.exit()
    samples = {}

    "+++ Shufling point IDs +++"
    for _ in range(100):
        random.shuffle(BIG_TRAIN)    
        random.shuffle(BIG_VAL)    
        random.shuffle(BIG_TEST)    

    "+++ Selecting points +++"
    if (len(list_props)>=1 and FInFile!=None):
        count_tasks = 0
        for tr_set_size in train:
            for val_set_size in val:
                for test_set_size in test:
                    samples[count_tasks] = []
                    samples[count_tasks].append(BIG_TRAIN[:tr_set_size])
                    samples[count_tasks].append(BIG_VAL[:val_set_size])
                    samples[count_tasks].append(BIG_TEST[:test_set_size])
                    count_tasks+=1

    "+++ Write files +++"
    print("POLUS| Program will now start writing individual files")
    for key in samples.keys():
        train     = samples[key][0]
        val       = samples[key][1]
        test      = samples[key][2]
        destin    = os.path.join(output_files_dir,"RS-"+str(len(train))+"-"+str(len(val))+"-"+str(len(test)))
        if not os.path.isdir(destin):
            os.mkdir(destin)
        tr_file   = open(os.path.join(output_files_dir,destin+"/INDICES-TRAIN.idx"),"w")
        val_file  = open(os.path.join(output_files_dir,destin+"/INDICES-VAL.idx"),"w")
        test_file = open(os.path.join(output_files_dir,destin+"/INDICES-TEST.idx"),"w")
        for value in train:
            tr_file.write(str(value)+"\n")
        for value in val:
            val_file.write(str(value)+"\n")
        for value in test:
            test_file.write(str(value)+"\n")
        tr_file.close()
        val_file.close()
        test_file.close()
   
        InFiles = os.listdir(input_files_dir) 
        for atom in list_atoms:
            for  prop in list_props:
                input_filename = None
                for file_ in InFiles:
                    if (atom in file_.split("_")):
                            input_filename = os.path.join(input_files_dir,file_)
                print (f"POLUS| Program is performing task: RS [{atom}]-{prop}") 
                outdir = os.path.join(destin,prop)
                tr_set_filename = os.path.join(outdir,system_name+"_"+atom+"_TRAINING_SET.csv")
                ival_set_filename = os.path.join(outdir,system_name+"_"+atom+"_INT_VALIDATION_SET.csv")
                eval_set_filename = os.path.join(outdir,system_name+"_"+atom+"_EXT_VALIDATION_SET.csv")   
                if check_property(input_filename,prop):
                    job = RS(input_filename, prop)  # does not matter which prop
                    job.set_training_set(train)
                    job.set_validation_set(val)
                    job.set_test_set(test)
                    if (len(train)>=1):
                        job.write_data_set(tr_set_filename,len(train),"Train",all_prop,prop)
                    if (len(val)>=1):
                        job.write_data_set(ival_set_filename,len(val),"Valid",all_prop,prop)
                    if (len(test)>=1):
                        job.write_data_set(eval_set_filename,len(test),"Test",all_prop,prop)
                else:
                    print(f"POLUS| Program complains::: RS [{atom}]-{prop} skipped due to missing/unrequired target property")

    wall_time = round(time.time() - start,3)
    print (f" Wall-time(s): {wall_time}")
    print_goodbye_message()


