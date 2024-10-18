import os
import sys


def check_file(path: str) -> bool:
    """
    This function checks the existence of a file
    
    Parameters:
    - path: str -> Path of the file to check
    """
    
    if os.path.isfile(path):
        exist=True
    else:
        exist=False
        
    return exist

def get_filesize(filename: str, unit: str = "bytes") -> None:
    """
    This function evaluates and returns the size of a file
    
    Parameters:
    - filename: str -> Path to the file whose size must be determined
    - unit:     str -> unit to use to express the file's size. Supported units include bytes, Kb and MB
    """
    
    if check_file(filename):
        size = os.path.getsize(filename)
        if (unit=="bytes"):
            size=size
        elif (unit=="Kb" or unit=="Kbytes"):
            size=round(size/1024,2)
        elif (unit=="MB" or unit=="MegaBytes"):
            size=round(size/1048576,2)
        else:
            size=size
    else:
        print("POLUS| Program complains::: File {} does not exist".format(filename))
        print_goodbye_message()
        sys.exit()
        
    return size

def check_property(filename: str, prop: str = "iqa") -> bool:
    """
    This function checks the existence of the data associated with the target physical properties. It basically reads 
    the file's header and locate a column that matches the property's name.
    
    Parameters:
    - filename: str -> Path to the file to be checked.
    - prop:     str -> Name of the target physical property.
    """

    if check_file(filename):    
        with open(filename,'r') as myfile:
            file_content = myfile.readlines()
            header       = file_content[0]
        if (prop.upper() in header.upper()):
            found = True
        else:
            found = False
    else:
        found = False
        print("Input file not found")

    return found
        
def readfile(filename: str, prop: str = "iqa"):
    """
    This function reads the content of a file and returns many components including the header, the file body, a vector of target property values,
    feature column indices, property column indicie and the indices of the remaining columns (possibly hosting data for other physical properties). 
    
    Parameters:
    - filename: str -> Path to the file to read.
    - prop:     str -> Name of the target physical property.
    """
    prop_vector = []
    features_indices = []
    prop_index = None
    if check_file(filename):   
        is_open_elsewhere = True
        i = 0
        while is_open_elsewhere: 
            try:
                myfile = open(filename,"r")
                is_open_elsewhere = False
            except Exception as e:
                i+=1
                if (i<100):
                    print("POLUS| Program could not open the file: Attempt ({})".format(i))
                else:
                    print("POLUS| Program could not open the file {}".format(filename))
                    print_goodbye_message()
                    sys.exit()

        file_content = myfile.readlines()
        myfile.close()
        header       = file_content[0]
        file_body    = file_content[1:]
        if (prop.upper() in header.upper()):
            for i in range(len(header.split(','))):
                cur_header_term = header.split(',')[i]
                if (prop.upper() in cur_header_term.upper()):
                    prop_index=i
            if prop_index == None:
                print("POLUS! Program could not find property of interest")
                SAR_goodbye()
                sys.exit()
                
        for i in range(len(file_body)):
            line = file_body[i]
            if prop_index==None:
                print("POLUS| Program is unable to compute prop-index. Job halted! ")
                SAR_goodbye()
                sys.exit()
            else:
                prop_vector.append(eval(line.split(',')[prop_index]))
            
        header_list=header.split(',')
        for j in range(len(header_list)):
            if (header_list[j][0]=='f'):
                features_indices.append(j)
        all_prop_indices = []       
        for value in range(len(header.split(','))):
            if value not in features_indices:
                all_prop_indices.append(value)

    else:
        sys.exit("POLUS| Program complains::: File {} does not exist".format(filename))
            
    return header, file_body, prop_vector, features_indices, prop_index, all_prop_indices

def read_job_details() -> tuple[str, str, int, int, int, list[str], list[str], str, str, object]:
    system_name = None
    strat_method = None
    tr_set_size = None
    ival_set_size = None
    eval_set_size = None
    job_detail_file = None
    atoms = None
    props = None
    PS_init_method = None
    PS_uncertainty_method = None
    recov_filter = 1.0
    cwd = os.getcwd()
    
    for file in os.listdir(cwd):
        if file.endswith(".txt") and "JOB" in file:
            job_detail_file = os.path.join(cwd,file)
    
    if job_detail_file == None:
        sys.exit("Job file not found")
        
    is_open_elsewhere = True
    i = 0
    while is_open_elsewhere:
        try:
            job_file = open(job_detail_file,"r")
            is_open_elsewhere = False
        except Exception as e:
            i+=1
            if (i<100):
                print(f" Could not open the file: Attempt ({i})")
            else:
                print(f" Could not open the file {filename}")
                SAR_goodbye()
                sys.exit()
         
    content = job_file.readlines()
    for line in content:
        line = line.split()
        if (len(line)<=1):
            continue
        else:
            if "system_name" == line[0]:
                system_name = line[1]
            if "strat_method" == line[0]:
                strat_method = line[1]
            if "tr_set_size" == line[0]:
                tr_set_size = int(line[1])
            if "ival_set_size" == line[0]:
                ival_set_size = int(line[1])
            if "eval_set_size" == line[0]:
                eval_set_size = int(line[1])
            if "atoms"== line[0]:
                atoms = line[1:]
            if "props" == line[0]:
                props = line[1:]
            if "PS_init_method"==line[0]:
                PS_init_method = line[1]
            if "PS_uncertainty_method"==line[0]:
                PS_uncertainty_method = line[1]
            if "recov_filter"==line[0]:
                recov_filter = eval(line[1])
       
    result = [system_name,strat_method,tr_set_size,ival_set_size,\
              eval_set_size, atoms, props,PS_init_method, PS_uncertainty_method,\
              recov_filter]
    if None in result:
        sys.exit("Could not read job details. Check file for  missing entries!")        

    return result


def check_input_files() -> list[bool]:
    #TODO: should take as input the actual input dir
    list_exists = []
    job_detail_file = "JOB_DETAILS.txt"
    if os.path.isfile(job_detail_file):
        list_exists.append(True)
        job_details = read_job_details()
        system_name = job_details [0]
        atoms = job_details[5]
        list_cwd = os.listdir(os.getcwd())
        list_input_dir = os.listdir(os.path.join(os.getcwd(),"input_files"))
        for atom in atoms:
            filename1 = system_name.upper()+"_"+atom.upper()+".csv"
            filename2 = atom.upper()+"_features_with_properties.csv"
            if (filename1 in list_input_dir or filename2 in list_input_dir):
                list_exists.append(True)
            elif (filename1 in list_cwd or filename2 in list_cwd):
                list_exists.append(True)
            else:
                list_exists.append(False)
    else:
        list_exists.append(False)

    return list_exists
def SAR_goodbye() -> None:
    print (" ")
    print (" Thanks for using POLUS (Sampling Algorithms for Regression).")
    print (" For more details: bienfait.isamura@postgrad.manchester.ac.uk")
    print (" ")

