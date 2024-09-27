import os
import shutil
import random
import math
import copy
import time
import numpy as np
from polus.distributions.commons import Files
from polus.utils.logging import PrintInfo, RaiseError
from polus.utils.printing import PrintOnTerminal, PrintGBMessage
from scipy.stats import kstest, cramervonmises_2samp, epps_singleton_2samp, anderson_ksamp
from scipy.stats import PermutationMethod
from scipy.special import kl_div


class ANALYSER(Files):
    """
    This class analyses distributions of target physical properties for every atom in the system.

    Attributes:
    - inputDir:          str       -> Path to the directory containing input .cvs files
    - atoms:             list(str) -> List of atom names, each written as a str object 
    - prop:              str       -> Target property
    - compareMethod:     str       -> Name of the statistical method to use to compare pairs of distributions
    - systemNam:         str       -> Name of the system (molecule or cluster)
    - outputDir:         str       -> Path to the output directory
    - xyzTraj:           str       -> Path to an xyz trajectory whose structures led to the distribution of atomic properties
    - combineMethod:     str       -> Method to be used when creating composite datasets for equivalent atoms
    - pValue:            float     -> p-value to be used in the decision making 

    """
    def __init__(self,inputDir=None,atoms=None,prop=None,statComputeMethod="auto",systemName=None,outputDir=None,xyzTraj=None,combineMethod=None,pValue=None,compareMethod=None,eqAtomsSM=None,smoothDistr=False):
        super().__init__(inputDir,atoms,prop,systemName,smoothDistr)
        self.inputDir         = inputDir
        self.outputDir        = outputDir
        self.pValue           = pValue
        self.xyzTraj          = xyzTraj
        self.smoothD          = smoothDistr
        self.combineMethod    = combineMethod
        self.statComputeMethod= statComputeMethod
        self.supportedMethods = ["KS","CVM","ES","KL"]
        self.method           = compareMethod
        self.eqAtomsSM        = eqAtomsSM
        self.distros          = None
        self.equivalentAtoms  = None
        self.remainingAtoms   = None
        self.atomTypes        = None
        self.testResults      = None
        self.eqAtomsSampleIDs = None
        self.atomTypeEquiv    = None
        self.atomTypesFile    = None

    def SetCombineMethod(self):
        if self.combineMethod==None:
            self.combineMethod="random"
        else:
            if isinstance(self.combineMethod,str):
                self.combineMethod = self.combineMethod
            else:
                RaiseError(message="Invalid combination method for equivalent atoms")

    def SetUserPValue(self):
        msg_         = "Setting User-defined pValue"
        start        = time.time()
        PrintInfo(message = msg_)
        PrintOnTerminal(msg = msg_)
        if self.pValue == None:
            self.pValue = 0.05
        PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def SetOutputDir(self):
        msg_         = "Setting Output Directory"
        start        = time.time()
        PrintInfo(message = msg_)
        PrintOnTerminal(msg = msg_)
        if self.outputDir == None:
            self.outputDir = os.path.join(os.getcwd(),"DISTROS")
        else:
            if isinstance(self.outputDir,str):
                self.outputDir = os.path.join(os.getcwd(),self.outputDir)
            else:
                RaiseError(message="Invalid Output Directory")
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)   

        PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))
    def SetEqAtomsSM(self):
        msg_         = "Setting EqAtoms Sampling Method"
        start        = time.time()
        PrintInfo(message = msg_)
        PrintOnTerminal(msg = msg_)
        if self.eqAtomsSM == None:
            self.eqAtomsSM == "RS"
        else:
            if isinstance(self.eqAtomsSM,str):
                self.eqAtomsSM = self.eqAtomsSM
            else:
                RaiseError(message="Invalid sampling method for equivalent atoms")
        PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def SetDistros(self):
        self.distros = self.SetDistributions()

    def SetMethod(self):
        if self.method == None:
            self.method = "KS"
        else:
            if isinstance(self.method,str):
                self.method = self.method
            else:
                RaiseError(message=" Invalid Comparison Method")
     

    def FindEquivalentAtoms(self):
        msg_         = "Finding Equivalent Atoms"
        start        = time.time()
        print("POLUS: Finding equivalent atoms")
        print(f"{'Atom1':>7} {'Atom2':>7}   {'p-value':>12}   {'Stat':>12}  {'Decision':>12}")
        #PrintInfo(message = msg_)
        #PrintOnTerminal(msg = msg_)
        if self.distros == None:
            self.SetDistros()
        if self.equivalentAtoms == None:
            self.testResults      = dict()
            self.equivalentAtoms  = dict()
            for element in list(self.atomGroups.keys()):
                atomGroup = self.atomGroups[element]
                for i in range(len(atomGroup)):
                    atom1     = atomGroup[i]
                    D1 = self.distros[element][i]
                    #e1 = np.max(D1) - np.min(D1)
                    D1       = [(X-np.mean(D1))/np.std(D1) for X in D1]
                    #D1       = [(X-np.min(D1))/e1 for X in D1]
                    for j in range(i+1,len(atomGroup)):
                        atom2    = atomGroup[j]
                        D2       = self.distros[element][j]
                        #e2       = np.max(D2) - np.min(D1)
                        # Transform D2
                        D2       = [(Y-np.mean(D2))/np.std(D2) for Y in D2]
                        #D2       = [(Y-np.min(D2)) for Y in D2]
                        #if self.prop == "iqa":
                        #    D1        = [X*2625.5 for X in D1]
                        #    D2        = [Y*2625.5 for Y in D2]
                        if self.statComputeMethod=="exact" and self.method!="AD":
                            D1,D2 = self.GetCompactDistributions(D1,D2)
                        testPair = self.ComparePairOfAtoms(sorted(D1),sorted(D2))
                        atomPair = atom1+"-"+atom2
                        self.testResults[atomPair] = list()
                        self.testResults[atomPair].append(testPair[0])
                        self.testResults[atomPair].append(testPair[1])
                        self.testResults[atomPair].append(testPair[2])
                        print(f"{atom1:<7} {atom2:<7}   {testPair[1]:>12.6e}   {testPair[2]:>12.6e}   {testPair[0]:>12}")
                        if testPair[0] == "E":
                            if atom1 in list(self.equivalentAtoms.keys()):
                                self.equivalentAtoms[atom1].append(atom2)
                            else:
                                self.equivalentAtoms[atom1] = [atom2]
        #PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def GetCompactDistributions(self,D1,D2):
        D1.sort()
        D2.sort()
        D1_, D2_ = list(), list()
        n1,n2=len(D1),len(D2)
        if n1>250 and n2>250:
            step = math.floor(n1/250)
            for i in range(250):
                D1_.append(D1[i*step])
                D2_.append(D2[i*step])
        return D1_,D2_

    def ComparePairOfAtoms(self,distro1,distro2):
        if self.method == None:
            self.SetMethod()
            if self.method not in self.supportedMethods:
                self.method = "KS"
        if isinstance(distro1,list) and isinstance(distro2,list):
            D1 = np.array(distro1)
            D2 = np.array(distro2)
            # Perform Test
            if self.method == "KS":
                test = kstest(D1,D2)
            elif self.method == "CVM":
                test = cramervonmises_2samp(D1,D2,method=self.statComputeMethod)
            elif self.method == "ES":
                test = epps_singleton_2samp(D1,D2)
            elif self.method == "AD":
                rng = np.random.default_rng()
                permMethod = PermutationMethod(n_resamples=999, random_state=rng)
                test = anderson_ksamp([D1,D2], method=permMethod)
            else:
                test = kl_div(D1,D2)
            # Set Test Result
            if not self.method == "KL":
                if test.pvalue >= self.pValue:
                    result = "E"
                else:
                    result = "NE"
            else:
                if np.sum(test) < 0.1 :
                    result = "NE"
                else:
                    result = "NE"
        else:
            RaiseError(message=" Invalid Distributions")
        if not self.method == "KL":
            outcome         = [result, test.pvalue, test.statistic]
        else: 
            outcome         = [result, 0.0, np.sum(test)]
        return outcome

    def WriteTestResults(self):
        if self.testResults == None:
            self.FindEquivalentAtoms()
        if self.outputDir   == None:
            self.SetOutputDir()
        outputFilename = os.path.join(self.outputDir,self.systemName.upper()+"-TEST-RESULTS-"+self.method.upper()+"-"+self.prop.upper()+".dat")
        outFile        = open(outputFilename,"w")
        print(f"POLUS: Test Results ...")
        print(f"{'Atom1':>7} {'Atom2':>7}   {'p-value':>12}   {'Stat':>12} {'Decision':>12}")
        outFile.write(f"{'Atom1':>7} {'Atom2':>7}   {'p-value':>12}   {'Stat':>12} {'Decision':>12}\n")
        #outFile.write("*"*45+"\n")
        for key in list(self.testResults.keys()):
            testSummary = self.testResults[key]
            atom1,atom2 = key.split("-")
            decision    = testSummary[0]
            pvalue      = testSummary[1]
            statistic   = testSummary[2]
            print(f"{atom1:>7} {atom2:>7}   {pvalue:>12.6e}   {statistic:>12.6e} {decision:>12}")
            outFile.write(f"{atom1:>7} {atom2:>7}   {pvalue:>12.6e}   {statistic:<12.6e} {decision:>12}\n")
        
        # TODO: Write list of equivalent atoms
        outFile.close()
           

    def SampleEquivalentAtoms(self):
        if self.eqAtomsSampleIDs == None:
            self.eqAtomSampleIDs = dict()
            if self.equivalentAtoms == None:
                self.FindEquivalentAtoms()
            if self.inputFiles  == None:
                self.SetInputFiles()
            if self.eqAtomsSM == None:
                self.SetEqAtomsSM()
            for element in list(self.equivalentAtoms.keys()):
                atomGroup         = self.equivalentAtoms[element]
                groupSize         = len(atomGroup)
                eqAtomsSampleSize = int(math.ceil(float(self.ngeoms)/float(groupSize)))
                self.eqAtomsSampleIDs[element] = list()
                if self.eqAtomsSM == "RS":
                    allGeomIDs     = [i for i in range(self.ngeoms)]
                    for _ in range(10):
                        random.shuffle(allGeomIDs)
                    self.eqAtomsSampleIDs[element] = allGeomIDs[:eqAtomsSampleSize]
                else:
                    RaiseError(message="Job aborted [Only RS is implemented for eq. atoms for the moment")
                
        else:
            if isinstance(self.eqAtomsSampleIDs,dict):
                self.eqAtomsSampleIDs = self.eqAtomsSampleIDs
            else:
                RaiseError(message="Invalid Sample IDs for equivalent atoms")

    def SetRemainingAtoms(self):
        if self.atoms == None:
            self.SetListOfAtoms()
        if self.equivalentAtoms == None:
            self.FindEquivalentAtoms()
        allEquivalentAtoms = list()
        for element in list(self.equivalentAtoms.keys()):
            for atom in self.equivalentAtoms[element]:
                allEquivalentAtoms.append(atom)
        self.remainingAtoms = list()
        for atom in self.atoms:
            if atom not in allEquivalentAtoms:
                self.remainingAtoms.append(atom)


    def SetAtomTypes(self):
        if self.equivalentAtoms == None:
            self.FindEquivalentAtoms()
        if self.remainingAtoms == None:
            self.SetRemainingAtoms()
        self.atomTypes = list()
        for atom in self.remainingAtoms:
            self.atomTypes.append(atom)
        for leadAtom in list(self.equivalentAtoms.keys()):
            atomGroup = leadAtom
            for atom in self.equivalentAtoms[leadAtom]:
                atomGroup = atomGroup + "-" + atom.upper()
            self.atomTypes.append(atomGroup)
        self.atomTypes    = self.RefineAtomTypes(self.atomTypes)
        print(f"POLUS: Printing atom types ...")
        print(f"POLUS: {'ID':>10} {'AtomType':>20}")
        count = 0
        for atomType in self.atomTypes:
            count = count + 1
            print(f"POLUS: {count:>10} {atomType:>20}")

    def RefineAtomTypes(self,atomTypes):
        realAtomGroups = list()
        singleAtoms    = list()
        for atomGroup in atomTypes:
            if len(atomGroup.split("-"))>1:
                realAtomGroups.append(atomGroup)
            else:
                singleAtoms.append(atomGroup)
        copyRealAtomGroups = list()
        for i in range(len(realAtomGroups)):
            rAG1    = realAtomGroups[i]
            atoms1  = rAG1.split("-")
            rAG1c   = rAG1
            for j in range(i+1,len(realAtomGroups)):
                identical = False
                rAG2      = realAtomGroups[j]
                atoms2    = rAG2.split("-")
                for atom in atoms1:
                    if atom in atoms2:
                        identical = True
                if identical:
                    rAG1c = rAG1c + "-" + rAG2
            copyRealAtomGroups.append(rAG1c)
        copyRealAtomGroups = self.CheckSubGroups(copyRealAtomGroups)
        if self.atoms == None:
            self.SetListOfAtoms()
        uniqueAtoms    = list()
        for atom in self.atoms:
            isUnique   = True
            for atomType in copyRealAtomGroups:
                if atom in atomType.split("-"):
                    isUnique = False
            if isUnique:
                uniqueAtoms.append(atom)
        AtomTypes = copyRealAtomGroups + uniqueAtoms
                    
        return AtomTypes
   
    def CheckSubGroups(self,atomGroups):
        finalAtomGroups_ = list()
        for i in range(len(atomGroups)):
            group1     = set(atomGroups[i].split("-"))
            isSubset_  = False
            for j in range(len(atomGroups)):
                if i!=j:
                    group2 = set(atomGroups[j].split("-"))
                    if group1.issubset(group2):
                        isSubset_ = True
            if not isSubset_:
                finalAtomGroups_.append(atomGroups[i])
        finalAtomGroups = list()
        for atomGroup in finalAtomGroups_:
            atomGroup = "-".join(list(set(atomGroup.split("-"))))
            finalAtomGroups.append(atomGroup)
        
        return finalAtomGroups

    def ReadAtomTypes(self,file_):
        self.atomTypes = list()
        with open(file_,"r") as myfile:
            content = myfile.readlines()
        for line in content:
            line = line[:-1]
            self.atomTypes.append(line)

    def CreateCompositeFiles(self,atomTypesFile=None,nGeoms=None): 
        if self.atomTypes  == None:
            if atomTypesFile== None:
                self.SetAtomTypes()
            else:
                self.ReadAtomTypes(atomTypesFile)
        if self.inputFiles == None:
            self.SetInputFiles()
        if self.combineMethod ==None:
            self.SetCombineMethod()
        if self.ngeoms == None:
           if nGeoms != None and isinstance(nGeoms,int):
               self.ngeoms = ngeoms
           else:
               self.SetListOfAtoms()
               self.ReadTargetProperty(atom=self.atoms[0])
        self.SetATEquivalence(self.atomTypes)
        newDir = os.path.join(self.outputDir,"NEWFILES") 
        if not os.path.isdir(newDir):
            os.mkdir(newDir)
        """ Write/copy unchanged files """ 
        countAGs = 0
        AGFile = open(os.path.join(newDir,self.systemName+"-AGs.dat"),"w")
        for key in self.atomTypeEquiv:
            if "AG" not in key:
                origFile = self.FindOrigFile(key,self.inputFiles)
                if origFile == None:
                    print("POLUS: Unable to find original file for atom "+key)
                    RaiseError(message="Unable to find original file for atom "+key)
                else:   
                    fileTail    = os.path.split(origFile)[1]
                    destination = os.path.join(newDir,fileTail)
                    shutil.copy(origFile,destination)
            else:
                countAGs  = countAGs + 1
                AG        =  self.atomTypeEquiv[key]
                kw        = "AG"+str(countAGs)
                atoms     = " ".join(AG)
                AGFile.write(f"{kw:<10} {atoms:<30}\n")
                nAtoms    = len(AG)
                shareSize = math.floor(self.ngeoms/nAtoms)
                complSize = self.ngeoms - nAtoms*shareSize
                if self.combineMethod=="random":
                    tempFile_   = self.PerformRandomCombination(AG,shareSize,complSize,self.inputFiles)
                    #fileTail    = os.path.split(tempFile_)[1].split(".csv")[0]
                    #newfileTail = fileTail.split("TEMP-FILE-")[1]+"_processed_data.csv" 
                    #destination = os.path.join(newDir,newfileTail)
                    destination = os.path.join(newDir,"AG"+str(countAGs)+"_processed_data.csv")
                    shutil.copy(tempFile_,destination)
                    os.remove(tempFile_)
        AGFile.close()

    def PerformRandomCombination(self,AG,shareSize,complSize,origFiles):
        random.seed(1)
        allGeoms = [i for i in range(self.ngeoms)]
        for _ in range(10):
            random.shuffle(allGeoms)
        CommonGeomIDs = allGeoms[:shareSize]
        OtherGeomIDs  = allGeoms[shareSize:shareSize+complSize]
        fileBodies    = list()
        combinedBody  = list()
        for atom in AG:
            origFile = self.FindOrigFile(atom,origFiles)
            with open(origFile,"r") as myfile:
                content = myfile.readlines()
            header     = content[0]
            body       = content[1:]
            fileBodies.append(body)
            for r in range(len(body)):
                if r in CommonGeomIDs:
                    combinedBody.append(body[r])
        for s in range(len(fileBodies[0])):
            if s in OtherGeomIDs:
                combinedBody.append(fileBodies[0][s])    
        tempFilename = os.path.join(os.getcwd(),"TEMP-FILE-"+"-".join(AG)+".csv")
        with open(tempFilename,"w") as tempFile:
            tempFile.write(header)
            for line in combinedBody:
                tempFile.write(line)
        return tempFilename
         

    def FindOrigFile(self,keyword,files):
        origFile = None
        for file_ in files:
            tail = os.path.split(file_)[1]
            if keyword in tail.split("_"):
                origFile = file_
                break
        return origFile

    def SetATEquivalence(self,ATypes):
        self.atomTypeEquiv = dict()
        countATGroups      = 0
        for AT in ATypes:
            if len(AT.split("-"))<=1:
                self.atomTypeEquiv[AT] = list(AT)
            else:
                countATGroups +=1
                self.atomTypeEquiv["AG"+str(countATGroups)] = AT.split("-")

    def WriteAtomTypes(self):
        if self.atomTypes == None:
            self.SetAtomTypes()
        self.atomTypesFile  = os.path.join(self.outputDir,self.systemName+"-ATOM-TYPES.dat")
        with open(self.atomTypesFile,"w") as myfile:
            for AT in self.atomTypes:
                myfile.write(AT+"\n")
        
    def Execute(self):
        self.SetCombineMethod()
        self.SetUserPValue()
        self.SetOutputDir()
        self.SetEqAtomsSM()
        self.SetMethod()
        self.SetDistros()
        self.FindEquivalentAtoms()
        self.SetRemainingAtoms()
        self.SetAtomTypes()
        self.WriteTestResults()
        self.WriteAtomTypes()
