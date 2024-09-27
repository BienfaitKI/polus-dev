import os
import time
import numpy as np
from polus.utils.logging import RaiseError, PrintInfo
from polus.utils.userInputs import ReadAtomLabels
from polus.utils.read_module import readfile
from polus.utils.printing import PrintOnTerminal
from scipy.stats import gaussian_kde



class Files():
    def __init__(self,inputDir=None,atoms=None,prop=None,systemName=None,smoothDistr=False):
        self.inputDir          = inputDir
        self.prop              = prop
        self.systemName        = systemName
        self.atoms             = atoms
        self.smoothD           = smoothDistr
        self.supportedElements = ["C","H","O","N","F","Cl","P","S"]
        self.digits            = ["0","1","2","3","4","5","6","7","8","9"]
        self.inputFiles        = None
        self.atomGroups        = None
        self.distributions     = None
        self.ngeoms            = None
        self.nAtomGroups       = None

    def SetInputDir(self):
        if self.inputDir  == None:
            self.inputDir  = os.path.join(os.getcwd(),"input_files")
        else:
            if isinstance(self.inputDir,str):
                self.inputDir = os.path.join(os.getcwd(),self.inputDir)
            else:
                RaiseError(message=" Invalid input directory ")
     
    def SetProp(self):
        if self.prop == None:
            self.prop = "iqa"
        else:
            if isinstance(self.prop,str):
                self.prop = self.prop
            else: 
               RaiseError(message=" Invalid target property")

    def SetListOfAtoms(self):
        if self.atoms      == None:
            self.atoms      = ReadAtomLabels(self.inputDir)
            if (self.atoms == None):
               RaiseError(message=" Missing atoms")
        else:
            if isinstance(self.atoms,list) and isinstance(self.atoms[0],str):
                self.atoms      = self.atoms
            else:
               RaiseError(message=" Invalid atoms")

    def SetInputFiles(self):
        if self.inputFiles == None:
            if self.inputDir == None:
                self.SetInputDir()
            csvFiles = list()
            for file in os.listdir(self.inputDir):
                if file.endswith(".csv"):
                    csvFiles.append(file)
            if self.atoms == None:
                self.SetListOfAtoms()
            self.inputFiles = list()
            for atom in self.atoms:
                for file in csvFiles:
                    if atom.upper()==self.GetAtomName(file):
                        self.inputFiles.append(os.path.join(self.inputDir,file))
    
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

    def GetElementName(self,label):
        element = ""
        for char in label:
            if ord(char.upper())>=65 and ord(char.upper())<=90:
                element = element+char.upper()
            else:
                break
        return element

    def SetAtomGroups(self):
        if self.atomGroups == None:
            self.atomGroups = dict()
            if self.atoms == None:
                self.SetListOfAtoms()
            for atom in self.atoms:
                element   = self.GetElementName(atom)   
                if element not in list(self.atomGroups.keys()):
                    self.atomGroups[element] = [atom]
                else:
                    self.atomGroups[element].append(atom)
        else:
            if isinstance(self.atomGroups,dict):
                self.atomGroups = self.atomGroups
            else:
                RaiseError(message="Program cannot define atom groups ")
        self.nAtomGroups = len(list(self.atomGroups.keys()))

    def SetDistributions(self):
        msg_         = "Setting distributions of target property"
        start        = time.time()
        PrintInfo(message = msg_)
        PrintOnTerminal(msg = msg_)
        if self.inputFiles == None:
            self.SetInputFiles()
        if self.atomGroups == None:
            self.SetAtomGroups()
        if self.distributions == None:
            self.distributions = dict()
            for element in list(self.atomGroups.keys()):
                self.distributions[element] = list()
                for atom in self.atomGroups[element]:
                    D = self.ReadTargetProperty(atom)
                    if self.smoothD:
                        self.distributions[element].append(self.SmoothDistribution(D))
                    else:
                        self.distributions[element].append(D)
            
        else:
            if isinstance(self.distributions,dict):
                self.distributions = self.distributions
            else:
                RaiseError(message=" Invalid distributions")
        PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))
        return self.distributions

    def ReadTargetProperty(self,atom):
        if self.inputFiles == None:
            self.SetInputFiles()
        fileOfInterest = None
        for file_ in self.inputFiles:
            testAtom = self.GetAtomName(file_)
            if atom.upper()==testAtom.upper():
                fileOfInterest = file_
                break
        if fileOfInterest == None:
            RaiseError(message="Program Cannot Determine File of Interest")
        else:
            fileContent = readfile(fileOfInterest, prop=self.prop)
            propVector  = fileContent[2]
        self.ngeoms     = len(propVector)
        return propVector



    def SmoothDistribution(self,X):
        if not isinstance(X,np.ndarray):
            X     = np.array(X)
        kde       = gaussian_kde(X,bw_method="silverman")
        positions = np.linspace(start=np.min(X),stop=np.max(X),num=5000)
        Xapprox   = kde.evaluate(positions)
         
        return Xapprox
        
