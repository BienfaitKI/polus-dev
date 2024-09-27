import os
import sys


def GetInputBasename(system_name, atom, input_dir):
    basename    = None
    basename1   = system_name+"_"+atom.upper() 
    basename2   = atom.upper()+"_features_with_properties"
    basename3   = atom.upper()+"_processed_data"
    input_files = os.listdir(input_dir)
    for file in input_files:
        if file.endswith(".csv"):
            tail = os.path.split(file)[1]
            if basename1 in file.strip():
                basename = file
                break
            if basename2 in file.strip():
                basename = file
                break
            if basename3 in file.strip():
                basename = file
                break
            if atom.upper() in tail.split("_"):
                baseline = file
    return basename

def GetListInputFiles(input_dir):
    files = list()
    if os.path.isdir(input_dir):
        for file_ in os.listdir(input_dir):
            if file_.endswith(".csv"):
                files.append(os.path.join(input_dir,file_))
    else:
        sys.exit("POLUS: Program complains::: {} is not a directory ".format(input_dir))

    return files
