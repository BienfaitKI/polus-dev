import os
from polus.utils.logging import PrintInfo


IQA = "iqa"
Q00 = "q00"

class constants():
    def __init__(self):
        self.root               = os.getcwd()
        self.default_indir      = os.path.join(self.root,"input_files")
        self.default_outdir     = os.path.join(self.root,"datasets") 
        self.default_filtdir    = os.path.join(self.root,"filtered") 
        self.default_sampler    = "RS"
        self.available_samplers = ["RS","PS","SRS","UESS","UESS2","UESS3"]
        self.supportedElements  = ["C","H","O","N","F","Cl","P","S","AG","AT"]
        self.default_train      = [2000] 
        self.default_valid      = [500]
        self.default_test       = [1000] 
        self.default_ncmd_args  = 10
        
        PrintInfo(message=' Default parameters loaded')


 

