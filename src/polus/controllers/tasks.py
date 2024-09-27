import os
import sys
from polus.controllers.managers import run_SRS, run_PS, run_UESS, run_UESS2, run_UESS3, run_iqa_corr
from polus.controllers.rs_manager import run_RS
from polus.controllers.recov_manager import run_recov_calcs
from polus.controllers.recov_manager2 import run_recov_q00_calcs
from polus.controllers.recov_manager3 import run_recov_dual_calcs



def task_controller(ncmd_args,input_dir,output_dir,sampler,filter_,atoms,props,train,val,test):
    if ncmd_args==1:
        run_RS(input_dir,output_dir,atoms,props,train,val,test)
    else:
        if (sampler=="SRS"): 
            run_SRS(input_dir,output_dir,atoms,props,train,val,test)
        elif (sampler=="PS"):
            run_PS(input_dir,output_dir,atoms,props,train,val,test)
        elif (sampler=="UESS"):
            run_UESS(input_dir,output_dir,atoms,props,train,val,test)
        elif (sampler=="UESS2"):
            run_UESS2(input_dir,output_dir,atoms,props,train,val,test)
        elif (sampler=="UESS3"):
            run_UESS3(input_dir,output_dir,atoms,props,train,val,test)
        elif (sampler=="RS"):
            run_RS(input_dir,output_dir,atoms,props,train,val,test)
        else:
            if (filter_=="IQA-FILTER"):
                run_recov_calcs()
            if (filter_=="Q00-FILTER"):
                run_recov_q00_calcs()
            if (filter_=="DUAL-FILTER"):
                run_recov_dual_calcs()
            elif (filter_=="iqa-correction"):
                run_iqa_corr(output_dir,input_dir)
            else:
                sys.exit("POLUS| Program aborted without performing any job")

