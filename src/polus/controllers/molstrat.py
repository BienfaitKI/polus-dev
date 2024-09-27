#!/usr/bin/env python

import os
import sys
from samplers.stratRS import SRS



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
