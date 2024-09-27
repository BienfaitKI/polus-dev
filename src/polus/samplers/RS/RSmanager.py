import time
import sys
import os
import shutil
import random
import logging
from polus.samplers.RS.randomSampling import RS
from polus.utils.printing import *
from polus.utils.read_module import check_property
from polus.utils.io import GetRandInFile
from polus.utils.defaults import constants
from polus.utils.userInputs import ReadAtomLabels, GetProps
from polus.utils.logging import RaiseError, RaiseWarning, PrintInfo
from polus.files.outputs import WriteJobDetails


global logger
logger = logging.getLogger(__name__)

class RSampler():
    def __init__(self,inputDir=None,outputDir=None,atoms=None,props=None,valTest=False,trainSize=None,validSize=None,testSize=None,systemName=None,allProp=False):
        self.inputDir    = inputDir    
        self.outputDir   = outputDir    
        self.atoms       = atoms
        self.props       = props
        self.trainSize   = trainSize
        self.validSize   = validSize
        self.testSize    = testSize
        self.systemName  = systemName
        self.allProp     = allProp
        self.defaults    = constants()
        self.start       = time.time()
        self.valTest     = valTest
        print (f"POLUS: RSampler initialized") 
    
    def SetInputFileDir(self):
        """
        Checks the user input directory and sets it to default if not defined yet.
        """
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
                RaiseWarning(message=" Invalid training set size. Value set to default.")
                self.trainSize  = vars(self.defaults)["default_train"]

    def SetValidSetSize(self):
        if self.validSize  == None:
            self.validSize  = vars(self.defaults)["default_valid"]
        else:
            if (isinstance(self.validSize,list)):
                self.validSize = self.validSize
            else: 
                RaiseWarning(message=" Invalid Validation set size. Value set to default")
                self.validSize  = vars(self.defaults)["default_valid"]

    def SetTestSetSize(self):
        if self.testSize  == None:
            self.testSize  = vars(self.defaults)["default_test"]
        else:
            if (isinstance(self.testSize,list)):
                self.testSize = self.testSize
            else: 
                RaiseWarning(message=" Invalid test set size. Value set to default")
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
        self.SetTrainSetSize()
        self.SetValidSetSize()
        self.SetTestSetSize()

    def PrintUserInputs(self): 
        PrintInfo(message=" INDIR  "+self.inputDir)
        PrintInfo(message=" OUTDIR "+self.outputDir)
        PrintInfo(message=" SYSTEM "+self.systemName)
        PrintInfo(message=" ATOMS  "+ " ".join([str(x) for x in self.atoms]))
        PrintInfo(message=" TRAIN  "+ " ".join([str(x) for x in self.trainSize]))
        PrintInfo(message=" VALID  "+ " ".join([str(x) for x in self.validSize]))
        PrintInfo(message=" TEST   "+ " ".join([str(x) for x in self.testSize]))

    def GenerateJobDetailsFile(self):
        WriteJobDetails(InDir=self.inputDir,\
                        OutDir=self.outputDir, \
                        system=self.systemName, \
                        atoms=self.atoms, \
                        props=self.props)

    def CheckDuplicates(self,train,valid,test):
        ntrain = len(train)
        nvalid = len(valid)
        ntest  = len(test)
        points = train.copy()
        for i in range(nvalid):
            points.append(valid[i])
        for i in range(ntest):
            points.append(test[i])
        dupl   = 0 
        for i in range(len(points)):
            for j in range(i+1,len(points)):
                if (points[i]==points[j]):
                    dupl+=1

        return dupl, 100*(dupl/len(points))

    def MolBasedSampling(self):
        lg_train_set = max(self.trainSize)
        lg_val_set   = max(self.validSize)
        lg_test_set  = max(self.testSize)
        FInFile      = GetRandInFile(self.inputDir)
        try:
            prob_job = RS(filename=FInFile,output_prop=self.props[0])  
            prob_job.get_training_point_IDs(lg_train_set)
            prob_job.get_validation_point_IDs(lg_val_set)
            prob_job.get_test_point_IDs(lg_test_set)
            BIG_TRAIN = prob_job.get_training_set()
            BIG_VAL   = prob_job.get_validation_set()
            BIG_TEST  = prob_job.get_test_set()
        except:
            RaiseError(message=" Program cannot perform molecular-wise random sampling with the prob-job")
        print (f"POLUS: Molecular-based sampling completed") 
        dupl,dupl_ratio=self.CheckDuplicates(BIG_TRAIN,BIG_VAL,BIG_TEST)
        print (f"POLUS: {'Dupl. points':<18} {dupl:>10} ") 
        print (f"POLUS: {'Dupl. Ratio(%)':<18} {dupl_ratio:>10.3f} ") 

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
                        if self.valTest:
                            for i in samples[count_tasks][2]:
                                samples[count_tasks][1].append(i)
                            samples[count_tasks][2] = samples[count_tasks][1].copy()
                             
                        count_tasks+=1
        return samples

    def Execute(self):
        self.CheckUserInputs()
        self.PrintUserInputs()
        self.GenerateJobDetailsFile()
        samples = self.MolBasedSampling()
        print (f"POLUS: Writing files (in progress)") 
        count_pools = 0
        for key in samples.keys():
            count_pools +=1
            train     = samples[key][0]
            val       = samples[key][1]
            test      = samples[key][2]
            print (f"POLUS: {'POOL':<10} {count_pools:>10}") 
            print (f"POLUS: {'TRAIN':<10} {len(train):>10}") 
            print (f"POLUS: {'VALID':<10} {len(val):>10}") 
            print (f"POLUS: {'TEST':<10} {len(test):>10}") 
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
            print (f"POLUS: {'Task':<10} {'Duration(s)':>10}") 
            for atom in self.atoms:
                for  prop in self.props:
                    input_filename = None
                    for file_ in InFiles:
                        if (atom in file_.split("_")):
                            input_filename = os.path.join(self.inputDir,file_)
                    monitor_task = time.time()
                    current_task = atom+"-"+prop 
                    outdir = os.path.join(destin,prop)
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
                        print(f"POLUS: Program complains::: RS [{atom}]-{prop} skipped due to missing/unrequired target property")
                    duration = time.time() - monitor_task
                    print (f"POLUS: {current_task:<10} {duration: >10.3f}") 

        wall_time = round(time.time() - self.start,3)
        PrintInfo(message=" Job terminated")
        PrintInfo(message=" Wall-time(s): "+str(wall_time))
        PrintInfo(message=" Thanks for using POLUS")
        PrintInfo(message=" Please report any bugs/unexpected behaviour to: bienfait.isamura@postgrad.manchester.ac.uk")


