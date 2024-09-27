import time
import sys
import os
import shutil
import random
import logging
from polus.samplers.PS.passiveSampling import PS
from polus.utils.printing import *
from polus.utils.read_module import check_property
from polus.utils.io import GetRandInFile
from polus.utils.defaults import constants
from polus.utils.userInputs import ReadAtomLabels, GetProps
from polus.utils.logging import RaiseError, RaiseWarning
from polus.config.user_inputs import read_config, read_cmd_args
from polus.files.outputs import WriteJobDetails


global logger
logger = logging.getLogger(__name__)

class PSampler():
    def __init__(self,InDir=None,OutDir=None,atoms=None,props=None,trainSize=None,validSize=None,testSize=None,systemName=None,allProp=False):
        self.inputDir    = InDir    
        self.outputDir   = OutDir    
        self.atoms       = atoms
        self.props       = props
        self.trainSize   = trainSize
        self.validSize   = validSize
        self.testSize    = testSize
        self.systemName  = systemName
        self.allProp     = allProp
        self.defaults    = constants()
        self.start       = time.time()
    
    def SetInputFileDir(self):
        if self.inputDir   == None:
            self.inputDir   = vars(self.defaults)["default_indir"] 
        else:
            self.inputDir   = self.inputDir

    def SetOutputFileDir(self):
        if self.outputDir  == None:
            self.outputDir  = vars(self.defaults)["default_outdir"] 
        else:
            self.outputDir  = self.outputDir

        if os.path.isdir(self.outputDir):
            shutil.rmtree(self.outputDir)
        os.mkdir(self.outputDir)

    def SetListOfAtoms(self):
        if self.atoms      == None:
            self.atoms      = ReadAtomLabels(self.inputDir)
            if (self.atoms == None):
               RaiseError(message="missing atoms")
        else:
            self.atoms      = self.atoms

    def SetListOfProps(self):
        if self.props      == None:
            self.props      = GetProps(allProps=True,momentRank=4)
        else:
            self.props      = self.props

    def SetTrainSetSize(self):
        if self.trainSize  == None:
            self.trainSize  = vars(self.defaults)["default_train"]
        else:
            if (isinstance(self.trainSize,list)):
                self.trainSize = self.trainSize
            else: 
                RaiseWarning(message="Invalid training set size. Value set to default.")
                self.trainSize  = vars(self.defaults)["default_train"]

    def SetValidSetSize(self):
        if self.validSize  == None:
            self.validSize  = vars(self.defaults)["default_valid"]
        else:
            if (isinstance(self.validSize,list)):
                self.validSize = self.validSize
            else: 
                RaiseWarning(message="Invalid Validation set size. Value set to default")
                self.validSize  = vars(self.defaults)["default_valid"]

    def SetTestSetSize(self):
        if self.testSize  == None:
            self.testSize  = vars(self.defaults)["default_test"]
        else:
            if (isinstance(self.testSize,list)):
                self.testSize = self.testSize
            else: 
                RaiseWarning(message="Invalid test set size. Value set to default")
                self.testSize  = vars(self.defaults)["default_test"]

    def SetSystemName(self):
        if self.systemName == None:
            self.systemName = "MyMOL"
        else:
            self.systemName = self.systemName

    def CheckUserInputs(self):
        self.SetInputFileDir()
        self.SetOutputFileDir()
        self.SetSystemName()
        self.SetListOfAtoms()
        self.SetListOfProps()
        self.SetListOfProps()
        self.SetTrainSetSize()
        self.SetValidSetSize()
        self.SetTestSetSize()

    def PrintUserInputs(self):
        print(f"{'INDIR':<15} {self.inputDir}")
        print(f"{'OUTDIR':<15} {self.outputDir}")
        print(f"{'SYSTEM':<15} {self.systemName}")
        print(f"{'TRAIN':<15} {self.trainSize}")
        print(f"{'VALID':<15} {self.validSize}")
        print(f"{'TEST':<15} {self.testSize}")

    def GenerateJobDetailsFile(self):
        WriteJobDetails(InDir=self.inputDir,\
                        OutDir=self.outputDir, \
                        system=self.systemName, \
                        atoms=self.atoms, \
                        props=self.props)

    def MolBasedSampling(self):
        lg_train_set = max(self.trainSize)
        lg_val_set   = max(self.validSize)
        lg_test_set  = max(self.testSize)
        FInFile      = GetRandInFile(self.inputDir)
        #try:
        prob_job = RS(filename=FInFile,output_prop=self.props[0])  
        prob_job.get_training_point_IDs(lg_train_set)
        prob_job.get_validation_point_IDs(lg_val_set)
        prob_job.get_test_point_IDs(lg_test_set)
        BIG_TRAIN = prob_job.get_training_set()
        BIG_VAL   = prob_job.get_validation_set()
        BIG_TEST  = prob_job.get_test_set()
        #except:
        #    RaiseError(message="Program cannot perform molecular-wise random sampling with the prob-job")

        samples = {}
        for _ in range(100):
            random.shuffle(BIG_TRAIN)    
            random.shuffle(BIG_VAL)    
            random.shuffle(BIG_TEST)    

        if (len(self.props)>=1 and FInFile!=None):
            count_tasks = 0
            for tr_set_size in self.trainSize:
                for val_set_size in self.validSize:
                    for test_set_size in self.testSize:
                        samples[count_tasks] = []
                        samples[count_tasks].append(BIG_TRAIN[:tr_set_size])
                        samples[count_tasks].append(BIG_VAL[:val_set_size])
                        samples[count_tasks].append(BIG_TEST[:test_set_size])
                        count_tasks+=1
        return samples

    def Execute(self):
        self.CheckUserInputs()
        self.PrintUserInputs()
        self.GenerateJobDetailsFile()
        samples = self.MolBasedSampling()
        for key in samples.keys():
            train     = samples[key][0]
            val       = samples[key][1]
            test      = samples[key][2]
            destin    = os.path.join(self.outputDir,"RS-"+str(len(train))+"-"+str(len(val))+"-"+str(len(test)))
            if not os.path.isdir(destin):
                os.mkdir(destin)
            tr_file   = open(os.path.join(destin,"INDICES-TRAIN.idx"),"w")
            val_file  = open(os.path.join(destin,"INDICES-VAL.idx"),"w")
            test_file = open(os.path.join(destin,"INDICES-TEST.idx"),"w")
            for value in train:
                tr_file.write(str(value)+"\n")
            for value in val:
                val_file.write(str(value)+"\n")
            for value in test:
                test_file.write(str(value)+"\n")
            tr_file.close()
            val_file.close()
            test_file.close()
   
            InFiles = os.listdir(self.inputDir) 
            for atom in self.atoms:
                for  prop in self.props:
                    input_filename = None
                    for file_ in InFiles:
                        if (atom in file_.split("_")):
                            input_filename = os.path.join(self.inputDir,file_)
                    print (f"POLUS| Program is performing task: RS [{atom}]-{prop}") 
                    outdir = os.path.join(destin,prop)
                    print(outdir)
                    tr_set_filename = os.path.join(outdir,self.systemName+"_"+atom+"_TRAINING_SET.csv")
                    ival_set_filename = os.path.join(outdir,self.systemName+"_"+atom+"_INT_VALIDATION_SET.csv")
                    eval_set_filename = os.path.join(outdir,self.systemName+"_"+atom+"_EXT_VALIDATION_SET.csv")   
                    if check_property(input_filename,prop):
                        job = RS(input_filename, prop)  # does not matter which prop
                        job.set_training_set(train)
                        job.set_validation_set(val)
                        job.set_test_set(test)
                        if (len(train)>=1):
                            job.write_data_set(tr_set_filename,len(train),"Train",self.allProp,prop)
                        if (len(val)>=1):
                            job.write_data_set(ival_set_filename,len(val),"Valid",self.allProp,prop)
                        if (len(test)>=1):
                            job.write_data_set(eval_set_filename,len(test),"Test",self.allProp,prop)
                    else:
                        print(f"POLUS| Program complains::: RS [{atom}]-{prop} skipped due to missing/unrequired target property")

        wall_time = round(time.time() - self.start,3)
        print (f" Wall-time(s): {wall_time}")
        print_goodbye_message()
        logger.info('Job terminated')


