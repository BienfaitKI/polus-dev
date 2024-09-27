import os
import time
import random
import copy
from polus.utils.logging import RaiseError, PrintInfo
from polus.utils.printing import PrintOnTerminal
from polus.samplers.RS.randomSampling import RS
from polus.utils.userInputs import ReadAtomLabels, GetProps
from polus.utils.read_module import check_property
from polus.utils.logging import RaiseWarning
from polus.files.outputs import WriteJobDetails


class SELECT():
    def __init__(self,atoms=None, props=None, valTest=True, allProp=False,indexFile=None,trainSize=1000,validSize=500,testSize=1000,systemName=None,excludedIndexFile=None,inputDir=None,outputDir=None,externalSet=None,considerExtSet=True,randomSeed=10):
        self.indexFile         = indexFile
        self.excludedIndexFile = excludedIndexFile
        self.inputDir          = inputDir
        self.outputDir         = outputDir
        self.trainSize         = trainSize
        self.validSize         = validSize
        self.testSize          = testSize
        self.systemName        = systemName
        self.atoms             = atoms
        self.props             = props
        self.allProp           = allProp
        self.externalSet       = externalSet
        self.considerExtSet    = considerExtSet
        self.randomSeed        = randomSeed
        self.valTest           = valTest
        self.supportedElements = ["C","H","O","N","S","F","P","Cl"]
        self.digits            = ["0","1","2","3","4","5","6","7","8","9"]
        self.IDs               = None
        self.allPointIDs       = None
        self.excludedIDs       = None
        self.filteredUserIDs   = None
        self.filteredPointIDs  = None
        self.filteredTrain     = None
        self.filteredValid     = None
        self.filteredTest      = None
        self.inputFiles        = None

    def SetExternalSet(self):
        if self.externalSet == None:
            self.externalSet = [-1]
        else:
            tempExternalSet  = list()
            if self.excludedIDs == None:
                self.ReadExcludedIndexFile()
            for i in self.externalSet:
                if i not in self.excludedIDs:
                    tempExternalSet.append(i)
            self.externalSet = copy.deepcopy(tempExternalSet)
            

    def ReadIndexFile(self):
        if not os.path.isfile(self.indexFile):
            RaiseError(message=" Program cannot find index file")
        else:
            if self.IDs  == None:
                with open(self.indexFile,"r") as myfile:
                     content = myfile.readlines()
                     if "#" in content[0]:
                         self.IDs = [int(x.split()[0]) for x in content[1:]]
                     else:
                         self.IDs = [int(x.split()[0]) for x in content]
            else:
                if isinstance(self.IDs[0],int):
                    self.IDs  = self.IDs
                else:
                    RaiseError(message=" Invalid Index Values")


    def ReadExcludedIndexFile(self):
        if not os.path.isfile(self.excludedIndexFile):
            RaiseWarning(message=" Program cannot find excluded index file")
            self.excludedIDs  = list()
        else:
            if self.excludedIDs  == None:
                with open(self.excludedIndexFile,"r") as myfile:
                     content = [line.split("\n")[0] for line in myfile.readlines()]
                     self.excludedIDs = [int(x) for x in content]
            else:
                if isinstance(self.excludedIDs[0],int):
                    self.excludedIDs  = self.excludedIDs
                else:
                    RaiseError(message=" Invalid Excluded Index Values")

    def SetSystemName(self):
        if self.systemName == None:
            self.systemName = "MOL"
        else:
            if isinstance(self.systemName,str):
                self.systemName = self.systemName.upper()
            else:
                RaiseError(message=" Invalid System Name")

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

    def SetInputDir(self):
        if self.inputDir == None:
            self.inputDir = os.getcwd()
        else:
            self.inputDir = self.inputDir
            if not os.path.isdir(self.inputDir):
                os.makedirs(self.inputDir)

        
    def SetOutputDir(self):
        if self.outputDir == None:
            self.outputDir = os.getcwd()
        else:
            self.oututDir = self.outputDir
            if not os.path.isdir(self.outputDir):
                os.makedirs(self.outputDir)

    def SetInputFiles(self):
        if self.inputFiles == None:
            self.inputFiles = list()
            if self.inputDir == None:
                self.SetInputDir()
            csvFiles = list()
            for file in os.listdir(self.inputDir):
                if file.endswith(".csv"):
                    csvFiles.append(file)
            if self.atoms == None:
                self.SetListOfAtoms()
            for atom in self.atoms:
                for file in csvFiles:
                    if atom.upper()==self.GetAtomName(file):
                        self.inputFiles.append(os.path.join(self.inputDir,file))

    def SetListPoints(self):
        if self.inputFiles == None:
            self.SetInputFiles()
        randomFile       = self.inputFiles[0]
        self.allPointIDs = list()
        with open(randomFile,"r") as myfile:
            ngeoms = len(myfile.readlines())-1
        for i in range(ngeoms):
            self.allPointIDs.append(i)
       
    def SetFilteredPointIDs(self):
        if self.allPointIDs == None:
            self.SetListPoints()
        if self.excludedIDs == None:
            self.ReadExcludedIndexFile()
        if self.filteredPointIDs == None: 
            self.filteredPointIDs = list()
            for ID in self.allPointIDs:
                if ID not in self.excludedIDs:
                    self.filteredPointIDs.append(ID)

    def SetFilteredTrain(self):
        if self.filteredPointIDs  == None:
            self.SetFilteredPointIDs()
        if self.IDs == None:
            self.ReadIndexFile()
        if self.filteredUserIDs == None:
            self.filteredUserIDs = list()
            for value in self.IDs:
                if value in self.filteredPointIDs:
                    self.filteredUserIDs.append(value)

        self.filteredTrain = copy.deepcopy(self.filteredUserIDs)
        if len(self.filteredTrain)<self.trainSize:
            trainComplement    = list()
            complementSize     = self.trainSize - len(self.filteredTrain)
            for value in self.filteredPointIDs:
                if value not in self.filteredTrain and value not in self.externalSet:
                    trainComplement.append(value)
            random.seed(self.randomSeed)
            for _ in range(100):
                random.shuffle(trainComplement)
            if len(trainComplement) >= complementSize:
                self.filteredTrain = self.filteredTrain + trainComplement[:complementSize]
            else:
                self.filteredTrain = self.filteredTrain + trainComplement[:]
        else:
             self.filteredTrain    = self.filteredTrain[:self.trainSize]
        
           
    def SetFilteredValid(self):
        if self.filteredPointIDs  == None:
            self.SetFilteredPointIDs()
        if self.filteredTrain == None:
            self.SetFilteredTrain()
        if not self.considerExtSet or self.externalSet==None:
            complementIDs = list()
            for value in self.filteredPointIDs:
                if value not in self.filteredTrain:
                    complementIDs.append(value)
            random.seed(self.randomSeed)
            for _ in range(100):
                random.shuffle(complementIDs)
            self.filteredValid = complementIDs[:self.validSize]
        else:
            if self.externalSet == [-1]:
                RaiseError(message="OOP! External set is missing ")
            if len(self.externalSet)>=self.validSize:
                self.filteredValid = self.externalSet[:self.validSize]
            else:
                self.filteredValid = self.externalSet[:]
                

    def SetFilteredTest(self):
        if self.filteredPointIDs  == None:
            self.SetFilteredPointIDs()
        if self.filteredTrain == None:
            self.SetFilteredTrain()
        if self.filteredValid == None:
            self.SetFilteredValid()

        complementIDs = list()
        if not self.considerExtSet or self.externalSet==None:
            if (self.trainSize + self.validSize + self.testSize) <= len(self.filteredPointIDs):
                for value in self.filteredPointIDs:
                    if value not in self.filteredTrain and value not in self.filteredValid:
                        complementIDs.append(value)
            else:
                for value in self.filteredPointIDs:
                    if value not in self.filteredTrain:
                        complementIDs.append(value)
            random.seed(self.randomSeed)
            for _ in range(100):
                random.shuffle(complementIDs)
            self.filteredTest = complementIDs[:self.testSize]
        else:
            if self.externalSet == [-1]:
                RaiseError(message="OOP! External set is missing ")
            if len(self.externalSet)>=(self.validSize+self.testSize):
                self.filteredTest = self.externalSet[self.validSize:self.validSize+self.testSize]
            elif len(self.externalSet)>(self.testSize):
                lowLimit = len(self.externalSet)-self.testSize
                self.filteredTest = self.externalSet[lowLimit:]
            else:
                RaiseError(message="OOP! Insufficient out-of-sample data")
            

    def GetAtomName(self,file_):
        atomName = None
        file_    = os.path.split(file_)[1]
        for entry in file_.split("_"):
            if len(entry)<=5:
                char1 = entry[0]
                if char1.upper() in self.supportedElements:
                    if (len(entry)==2) and entry[1] in self.digits: 
                       atomName = entry
                    if (len(entry)==3) and entry[2] in self.digits: 
                       atomName = entry
        if atomName == None:
            RaiseError(message="Program cannot define atom name ")

        return atomName

    def WriteFiles(self):
        print(f"POLUS: Writing files upon index-based sampling ...")
        if self.inputFiles == None:
            self.SetInputFiles()
        if self.atoms      == None:
            self.SetListOfAtoms()
        if self.props      == None:
            self.SetListOfProps()

        # Write FEREBUS job Details File
        WriteJobDetails(InDir=self.inputDir,OutDir=self.outputDir,system=self.systemName,atoms=self.atoms,props=self.props)

        # Check valTest
        if self.valTest:
            for value in self.filteredTest:
                self.filteredValid.append(value)
            self.filteredTest = self.filteredValid.copy()

        # Write Training, Validation and Test set files
        if isinstance(self.trainSize,int) and isinstance(self.validSize,int) and isinstance(self.testSize,int):
            for atom in self.atoms:
                #msg_   = atom.upper()
                #start  = time.time()
                #PrintInfo(message = msg_)
                #PrintOnTerminal(msg = msg_)
                for  prop in self.props:
                    input_filename = None
                    for file_ in self.inputFiles:
                        testAtom  = self.GetAtomName(file_)
                        if atom == testAtom:
                            input_filename = file_
                    current_task = atom+"-"+prop 
                    outdir = os.path.join(self.outputDir,prop)
                    tr_set_filename = os.path.join(outdir,self.systemName+"_"+atom+"_TRAINING_SET.csv")
                    ival_set_filename = os.path.join(outdir,self.systemName+"_"+atom+"_INT_VALIDATION_SET.csv")
                    eval_set_filename = os.path.join(outdir,self.systemName+"_"+atom+"_EXT_VALIDATION_SET.csv")   
                    if check_property(input_filename,prop):
                        job = RS(input_filename, prop)  # does not matter which prop
                        job.set_training_set(self.filteredTrain)
                        job.set_validation_set(self.filteredValid)
                        job.set_test_set(self.filteredTest)
                        if (len(self.filteredTrain)>=1):
                            job.write_data_set(tr_set_filename,len(self.filteredTrain),"Train",self.allProp,prop)
                        if (len(self.filteredValid)>=1):
                            job.write_data_set(ival_set_filename,len(self.filteredValid),"Valid",self.allProp,prop)
                        if (len(self.filteredTest)>=1):
                            job.write_data_set(eval_set_filename,len(self.filteredTest),"Test",self.allProp,prop)
                    #else:
                    #    RaiseWarning(message=f"Program complains::: RS [{atom}]-{prop} skipped due to missing/unrequired target property")
                #PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def Execute(self):
        self.ReadIndexFile()
        self.ReadExcludedIndexFile()
        self.SetExternalSet()
        self.SetSystemName()
        self.SetInputDir()
        self.SetOutputDir()
        self.SetInputFiles()
        self.SetListOfAtoms()
        self.SetListOfProps()
        self.SetFilteredTrain()
        self.SetFilteredValid()
        self.SetFilteredTest()
        self.WriteFiles()
