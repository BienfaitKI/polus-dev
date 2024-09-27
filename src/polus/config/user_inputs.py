import sys
import os
from configparser import ConfigParser
from polus.utils.printing import *

global config_file, cwd, elements, prop_itls
config_file = os.path.join(os.getcwd(),"polus.ini")
cwd         = os.getcwd()
elements    = ["C","H","O","N","P","S"]
prop_itls   = ["i","q"]

def read_atom_names(input_dir):
    if os.path.isdir(input_dir):
        files = os.listdir(input_dir)
        atoms = []
        for file in files:
            if file.endswith(".csv"):
                file = file.split("_")
                for entry in file:
                    if (entry[0] in elements and isinstance(eval(entry[1:]),int)):
                        atoms.append(entry)
    else:
        print("POLUS| Program complains::: Input directory not found")
    
    if (len(atoms)==0):
        print("POLUS| Program complains::: Unable to define list of atoms")
        sys.exit()

    return atoms

def read_config():
    system_name  = None
    atoms        = None
    props        = None
    method       = None
    train        = None
    val          = None
    test         = None
    iqa_filt     = None
    q00_filt     = None
    try:
        config=ConfigParser()
        config.read(os.path.join(config_file))
        system_name = config["system"]["name"]
        method      = config["sampling"]["method"]
        try:
            iqa_filt    = eval(config["filtering"]["iqa_filt"])
            q00_filt    = eval(config["filtering"]["q00_filt"])
        except:
            iqa_filt    = 1.0
            q00_filt    = 0.0005
        if (len(config["system"]["atoms"].split())>0): 
            atoms   = [entry for entry in config["system"]["atoms"].split()]
        if (len(config["system"]["props"].split())>0):
            props   = [entry for entry in config["system"]["props"].split()]
        train       = [int(entry) for entry in config["sampling"]["train"].split()]
        val         = [int(entry) for entry in config["sampling"]["val"].split()]
        test        = [int(entry) for entry in config["sampling"]["test"].split()]
    except BaseException as e:    
        sys.exit('Possible error in the polus.ini file.\n\n The execution of SELECTOR has been aborted')

    return system_name, atoms, props, method, train, val, test, iqa_filt, q00_filt


def read_cmd_args():
    list_cmd_args = sys.argv
    nbr_args      = len(list_cmd_args)
    sampler       = None
    filter_       = None
    input_dir     = None
    output_dir    = None
    atoms         = None
    props         = None
    method        = None
    train         = None
    val           = None
    test          = None
    system        = None
    #+++ Read cmd_args +++#
    for i in range(len(list_cmd_args)):
        arg=list_cmd_args[i]
        if arg=="-I" or arg=="--input-dir":
            input_dir = list_cmd_args[i+1]
            if input_dir == "." or input_dir=="cwd":
                input_dir== cwd
        if arg=="-O" or arg=="--output-dir":
            output_dir = list_cmd_args[i+1]
            if output_dir == "." or output_dir=="cwd":
                output_dir== cwd
        if arg=="-A" or arg=="--atom":
            atoms   = [entry for entry in list_cmd_args[i+1].split("_")]
        if arg=="-P" or arg=="--prop":
            props   = [entry for entry in list_cmd_args[i+1].split("_")]
        if arg=="-R" or arg=="--run":
            sampler = list_cmd_args[i+1]
        if arg=="-S" or arg=="--system":
            system = list_cmd_args[i+1]
        if arg=="--recovery":
            if (list_cmd_args[i+1].upper()=="IQA"):
                filter_ = "IQA-FILTER"
            elif (list_cmd_args[i+1].upper()=="Q00"):
                filter_ = "Q00-FILTER"
            elif (list_cmd_args[i+1].upper()=="DUAL"):
                filter_ = "DUAL-FILTER"
            else:
                print("POLUS| Program complains::: Recov. filter not defined. Value set to IQA-FILTER")
                filter_ = "IQA-FILTER"
        if arg=="-ciqa" or arg=="":
            filter_ = "iqa-correction"
        if arg=="-T" or arg=="--train":
            train   = [int(entry) for entry in list_cmd_args[i+1].split("_")]
        if arg=="-V" or arg=="--val":
            val     = [int(entry) for entry in list_cmd_args[i+1].split("_")]
        if arg=="-t" or arg=="--test":
            test    = [int(entry) for entry in list_cmd_args[i+1].split("_")]
        if arg=="--help" or arg=="-h":
            print_help()
            sys.exit()
        if arg=="--version" or arg=="-v":
            print_version()
            sys.exit()
        if arg=="--check" or arg=="-c":
            print_input_details()
            if nbr_args>=i+1 and list_cmd_args[i+1]=="kill":
                sys.exit()
            if arg=="--clean" or arg=="-C":
                flag = False
                if nbr_args>=i+1:
                    if isinstance(list_cmd_args[i+1],str) and list_cmd_args[i+1].upper()=="ALL":
                        flag= True
                remove_output_folders(flag)

    return input_dir, output_dir, atoms, props, sampler, filter_, train, val, test, nbr_args, system

