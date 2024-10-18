import sys
import os
from configparser import ConfigParser
from polus.utils.printing import *

global config_file, cwd, elements, prop_itls
config_file = os.path.join(os.getcwd(),"polus.ini")
cwd         = os.getcwd()
elements    = ["C","H","O","N","P","S","F","Cl"]
prop_itls   = ["i","q"]

def ReadAtomLabels(input_dir: str) -> list[object]:
    """
    This function reads the atom labels from the input .cvs filenames.
    
    Parameters:
    - input_dir: str -> Path to the directory containing the input .csv files
    """
    
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
        print("POLUS: Program complains::: Input directory not found")
    
    if (len(atoms)==0):
        print("POLUS: Program complains::: Unable to define list of atoms")
        sys.exit()

    return atoms

def GetProps(allProps: bool = False, momentRank: int = 4) -> list[str]:
    """
    This function generates the list of possible target property names.
    
    Parameters:
    - allProps:   Boolean -> Boolean variable indicating whether all the properties should be considered or not
    - momentRank  integer -> Maximum rank for the atomic multipole moments
    """

    if (allProps):
        props = ["iqa"]
        for i in range(momentRank+1):
            for j in range(0,i+1):
                if (j==0):
                    props.append("q"+str(i)+str(j))
                else:
                    props.append("q"+str(i)+str(j)+"c")
                    props.append("q"+str(i)+str(j)+"s")
    else:
        props = ["iqa"]
    return props

def read_config() -> tuple[str, list[object], list[object], str, list[int], list[int], list[int], float, float]:
    """
    This function (depricated) used to read the content of a user-supplied configuration file.
    
    """
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
        # TODO: Replace sys.exit() by polus.utils.logging.RaiseError()
        sys.exit('Possible error in the polus.ini file.\n\n The execution of SELECTOR has been aborted')

    return system_name, atoms, props, method, train, val, test, iqa_filt, q00_filt


def read_cmd_args() -> tuple[str | None, str | None, list[str] | None, list[str] | None, str | None, str | None, list[int] | None, list[int] | None, list[int] | None, int, str | None]:
    """
    This function (depricated) used to read command line arguments.
    
    """
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

