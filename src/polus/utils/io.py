import os
import sys
from polus.config.user_inputs import read_config
from polus.config.user_inputs import read_cmd_args
from polus.config.user_inputs import read_atom_names
from polus.files.inputs import GetInputBasename


def GetRandInFile(InDir):
    """
    This function reads a random file from the input Directory
    
    Parameters:
    - InDir: str -> Path to the directory containing input .csv files
    """
    if (os.path.isdir(InDir)):
        files = [os.path.join(InDir,x) for x in os.listdir(InDir)]
        for file_ in files:
            if IsCSV(file_):
                InFile = file_
                break
    else:
        RaiseError(message="Cannot find input directory")
    return InFile

def IsCSV(file_):
    """
    This function checks if a file has comma separated content
    
    Parameters:
    - file_: str -> Path to the file to be checked
    """
    return file_.endswith(".csv")

def get_FInFile_IO_Dirs(root,input_dir,output_dir):
    default_input_dir = os.path.join(root,"input_files")
    cmd_args          = read_cmd_args()
    config_details    = read_config()
    input_dir         = cmd_args[0]
    output_dir        = cmd_args[1]
    FIFile            = None
    "+++ Find first input file based on system_name+++" 
    config_details = read_config()
    system_name    = config_details[0]
        

    if input_dir==None or input_dir==root:
        input_files_dir = root
    else:
        input_files_dir = os.path.join(root,input_dir)
    if output_dir==None or output_dir==root:
        output_files_dir = root
    else:
        output_files_dir = os.path.join(root,output_dir)
    if input_files_dir!=None and not os.path.isdir(input_files_dir):
        sys.exit("POLUS| Program complains::: Cannot find input directory")
    if output_files_dir!=root and not os.path.isdir(output_files_dir):
        try:
            os.mkdir(output_files_dir)
        except FileExistsError:
            print(f"{os.path.join(cwd,output_dir)} already exists") 

    OutDir = output_files_dir

    "+++ Find atom names +++"
    atoms = None
    if (config_details[1]!=None):
        atom_1     = config_details[1][0]
    else:
        atoms      = read_atom_names(input_files_dir)
        atom_1     = atoms[0]

    "+++ Set base filename +++"
    basename = GetInputBasename(system_name,atom_1,input_files_dir)
    if (basename!=None):
        First_input_filename = os.path.join(input_files_dir,basename)
        InDir = input_files_dir
    else:
        basename = GeInputBasename(system_name,atom_1,default_input_dir)
        if (basename!=None):
            First_input_filename = os.path.join(default_input_dir,basename)
            InDir = default_input_dir
        else:
            print(" Cannot build basename. Job halted!")
            sys.exit()

    return First_input_filename, InDir, OutDir, atoms
