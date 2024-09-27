import os
import shutil
import sys
import time
from polus.utils.logging import RaiseError, RaiseWarning, PrintInfo
from polus.utils.read_module import readfile
from polus.files.inputs import GetListInputFiles
from polus.utils.defaults import constants, IQA, Q00
from polus.utils.userInputs import ReadAtomLabels

class recovEnergy():
    def __init__(self,atoms,systemName,workingDir=None,targetProp=None, geomIDs=None,inputDir=None,outputDir=None,dualFilter=None,threshold=None):
        self.targetProp             = targetProp 
        self.workingDir             = workingDir 
        self.inputDir               = inputDir 
        self.outputDir              = outputDir 
        self.geomIDs                = geomIDs
        self.atoms                  = atoms
        self.systemName             = systemName
        self.dualFilter             = dualFilter
        self.threshold              = threshold
        self.atomicProp             = IQA
        self.defaults               = constants()
        self.listWfnEnergies        = None
        self.atomicIqaEnergies      = None
        self.molecularIqaEnergies   = None
        self.recoveryEnergies       = None
        self.absRecoveryEnergies    = None
        self.listFiles              = None
        self.natoms                 = None
        self.filteredDir            = None
        print (f"POLUS: RecovEnergy initialized") 

    def SetInputFileDir(self):
        if self.inputDir   == None:
            self.inputDir   = vars(self.defaults)["default_indir"] 
        else:
            self.inputDir   = self.inputDir

    def SetOutputFileDir(self):
        if self.outputDir  == None:
            self.outputDir  = vars(self.defaults)["default_filtdir"] 
        else:
            self.outputDir  = self.outputDir

        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

    def SetListOfAtoms(self):
        if self.atoms      == None:
            self.atoms      = ReadAtomLabels(self.inputDir)
            if (self.atoms == None):
               RaiseError(message=" missing atoms")
        else:
            self.atoms      = self.atoms

    def SetSystemName(self):
        if self.systemName == None:
            self.systemName = "MyMOL"
        else:
            self.systemName = self.systemName

    def SetTargetProp(self):
        if self.targetProp == None:
            self.targetProp = "wfn_energy"
        else:
            self.targetProp = self.targetProp

    def SetWorkingDir(self):
        if self.workingDir == None:
            self.workingDir = os.getcwd()
        else:
            self.workingDir = self.workingDir

    def SetListFiles(self):
        if self.listFiles  == None:
            self.listFiles = GetListInputFiles(self.inputDir)
        else:
            self.listFiles = self.listFiles

    def SetGeomIDs(self):
        if self.geomIDs   == None:
            randomFile     = self.listFiles[0]
            datasetSize    = len(readfile(randomFile,prop="iqa")[1])
            self.geomIDs   = [i for i in range(datasetSize)]
        else:
            self.geomIDs   = self.geomIDs

    def SetNatoms(self):
        if self.natoms == None:
            self.natoms = len(self.atoms)
        else:
            self.natoms = self.natoms

    def SetDualFilter(self):
        if self.dualFilter == None:
            self.dualFilter = False
        else:
            if (isinstance(self.dualFilter,bool)):
                self.dualFilter = self.dualFilter
            else:
                RaiseWarning(message= " Invalid Dual filter parameter. Value set to False")
                self.dualFilter = False
               
    def SetThreshold(self):
        if (self.threshold == None):
            self.threshold = 1.0  # Expressed in kJ/mol
        else: 
            if (isinstance(self.threshold,float) or isinstance(self.threshold,int)):
                self.threshold = self.threshold
            else:
                RaiseWarning(message= " Invalid Threshold parameter. Value set to default (1.0 kJ/mol)")
                self.threshold = 1.0

    def CheckUserInputs(self):
        self.SetInputFileDir()
        self.SetOutputFileDir()
        self.SetSystemName()
        self.SetTargetProp()
        self.SetWorkingDir()
        self.SetListOfAtoms()
        self.SetListFiles()
        self.SetGeomIDs()
        self.SetNatoms()
        self.SetDualFilter()
        self.SetThreshold()

    def PrintUserInputs(self): 
        PrintInfo(message=" INDIR   "+self.inputDir)
        PrintInfo(message=" OUTDIR  "+self.outputDir)
        PrintInfo(message=" WORKDIR "+self.workingDir)
        PrintInfo(message=" SYSTEM  "+self.systemName)
        PrintInfo(message=" ATOMS   "+ " ".join([str(x) for x in self.atoms]))

    def ReadWfnEnergies(self):
        randomFile          = self.listFiles[0]
        inputData           = readfile(randomFile, self.targetProp)
        file_body, prop_idx = inputData[1], inputData[4]

        if self.listWfnEnergies == None:
            self.listWfnEnergies = list()
            for i in range(len(file_body)):
                if i in self.geomIDs:
                    line = file_body[i]
                    self.listWfnEnergies.append(eval(line.split(",")[prop_idx]))

    
    def ReadAtomicIqaEnergies(self):
        randomFile          = self.listFiles[0]
        datasetSize         = len(readfile(randomFile,prop=self.atomicProp)[1])

        if self.atomicIqaEnergies   == None:
            self.atomicIqaEnergies   = dict()
            count_atoms = 0
            for i in range(self.natoms):
                atom = self.atoms[i]
                self.atomicIqaEnergies[atom] = []
                data = readfile(self.listFiles[i],prop=self.atomicProp)
                file_content, iqa_idx = data[1], data[4]
                for j in range(len(file_content)):
                    if j in self.geomIDs:
                        line = file_content[j]
                        self.atomicIqaEnergies[atom].append(eval(line.split(",")[iqa_idx]))
                count_atoms +=1

        return self.atomicIqaEnergies

    def SetMolecularIqaEnergies(self):
        if self.atomicIqaEnergies==None:
            atomicIqaEnergies = self.ReadAtomicIqaEnergies()
        else:
            atomicIqaEnergies = self.atomicIqaEnergies
        
        "+++ Compute molecular IQA energies +++"
        if self.molecularIqaEnergies   == None:
            self.molecularIqaEnergies   = list() 
            for i in range(len(atomicIqaEnergies[self.atoms[0]])):
                mol_energy = 0
                for key in atomicIqaEnergies.keys():
                    mol_energy += atomicIqaEnergies[key][i]
                self.molecularIqaEnergies.append(mol_energy)
        else:
            self.molecularIqaEnergies   = self.molecularIqaEnergies() 


    def GenerateRecoveryEnergies(self):
        if self.listWfnEnergies      == None:
            self.ReadWfnEnergies()
        if self.molecularIqaEnergies == None:
            self.SetMolecularIqaEnergies()
        self.recoveryEnergies         = list()  # Defined as (wfn - iqa)
        self.absRecoveryEnergies      = list()
        if len(self.listWfnEnergies)==len(self.molecularIqaEnergies):
            for i in range(len(self.listWfnEnergies)):
                recov_energy = self.listWfnEnergies[i] - self.molecularIqaEnergies[i]
                self.recoveryEnergies.append(recov_energy)
                self.absRecoveryEnergies.append(abs(recov_energy))
        else:
            RaiseError(message=" number of WFN and IQA energy terms not matching")

        
        return self.recoveryEnergies, self.absRecoveryEnergies

    def writeFiles(self,outputFilename = None):
        geomIDs      = self.geomIDs
        dual_flag    = self.dualFilter
        if outputFilename!=None:
            filename = outputFilename
        else:
            filename = os.path.join(self.outputDir,self.systemName.upper()+"_RECOV_IQA")

        if self.recoveryEnergies == None:
            recov_energies, abs_recov_energies = self.GenerateRecoveryEnergies()
        else:
            recov_energies     = self.recoveryEnergies
            abs_recov_energies = [abs(x) for x in recov_energies]
        MOL_IQA    = self.molecularIqaEnergies
        MOL_WEn    = self.listWfnEnergies
        if recov_energies == None:
            RaiseError(message=" Program unable to complete task due to missing WFN energies ")
        with open(filename,"w") as rec_energy_file:
            #rec_energy_file.write(f"{'ID':>10} {'IQA-MOL':>15} {'WFN':>15} {'RECOV-En(Ha)':>15} {'ABS-RecEn(Ha)':>15} {'RECOV-En(kJ/mol)':>20} {'ABS-RecEn(kJ/mol)':>20}\n")
            rec_energy_file.write(f"#FIELDS: Geom-ID IQA-MOL WFN-En RECOV-En(Ha) ABS-RECOV-En(Ha) RECOV-En(kJ/mol) ABS-RECOV-En(kJ/mol)\n")
            for i in range(len(recov_energies)):
                rec_energy_file.write(f"{str(geomIDs[i]):<7} {MOL_IQA[i]:>12.6f} {MOL_WEn[i]:>12.6f} {recov_energies[i]:>10.6f} {abs_recov_energies[i]:>10.6f} {2625.5*recov_energies[i]:>10.6f} {2625.5*abs_recov_energies[i]:>10.6f} \n")

        if not dual_flag:
            badGeomsFile  = os.path.join(self.outputDir,self.systemName.upper()+"-BAD-GEOMETRIES-IQA")
            goodGeomsFile = os.path.join(self.outputDir,self.systemName.upper()+"-GOOD-GEOMETRIES-IQA")
            filteredDir   = os.path.join(self.outputDir,"FILTERED-BY-IQA")
        else:
            badGeomsFile  = os.path.join(self.outputDir,self.systemName.upper()+"-BAD-GEOMETRIES-DUAL")
            goodGeomsFile = os.path.join(self.outputDir,self.systemName.upper()+"-GOOD-GEOMETRIES-DUAL")
            filteredDir   = os.path.join(self.outputDir,"FILTERED-BY-DUAL")

        countBadGeoms      = 0
        countGoodGeoms     = 0
        goodGeomIDs        = list()
        myBadGeomsFile     = open(badGeomsFile,"w")
        myGoodGeomsFile    = open(goodGeomsFile,"w")
        myBadGeomsFile.write(f"#FIELDS: Geom-ID RECOV-En(kJ/mol)\n")
        myGoodGeomsFile.write(f"#FIELDS: Geom-ID RECOV-En(kJ/mol)\n")
        for i in range(len(recov_energies)):
            if (abs(2625.5*recov_energies[i])>self.threshold):
                countBadGeoms +=1
                myBadGeomsFile.write(f"{geomIDs[i]:<7} {2625.5*recov_energies[i]:>10.6f} \n")   
            else:
                countGoodGeoms +=1
                myGoodGeomsFile.write(f"{geomIDs[i]:<7} {2625.5*recov_energies[i]:>10.6f} \n")   
                goodGeomIDs.append(i)

        print("POLUS: Job completed successfully")
        print(f"POLUS: Total geometries               {len(geomIDs)}")
        print(f"POLUS: Removed geometries             {countBadGeoms}")
        print(f"POLUS: Retained geometries            {countGoodGeoms}")
        if not os.path.isdir(filteredDir):
            os.mkdir(filteredDir)
        if (self.filteredDir==None):
            self.filteredDir = filteredDir
        if (self.listFiles!=None):
            for file in self.listFiles:
                with open(file,"r") as oldFile:
                    content = oldFile.readlines()
                    header = content[0]
                    body = content[1:]
                newFilename = os.path.join(filteredDir,os.path.split(file)[1])
                with open(newFilename,"w") as newFile:
                    newFile.write(header)
                    for i in range(len(body)):
                        if i in goodGeomIDs:
                            newFile.write(body[i])

    def Execute(self):
        "+++ Check user inputs +++"
        self.CheckUserInputs()
        "+++ Print user inputs +++"
        self.PrintUserInputs()
        "+++ Read user inputs +++"
        print(f"POLUS: Reading Wfn Energies")
        self.ReadWfnEnergies()
        "+++ Read Atomic IQA energies +++"
        print(f"POLUS: Reading Atomic IQA nergies")
        self.ReadAtomicIqaEnergies()
        "+++ Set Molecular IQA energies +++"
        print(f"POLUS: Computing Molecular IQA energies")
        self.SetMolecularIqaEnergies()
        "+++ Generate Recovery errors +++"
        print(f"POLUS: Computing recovery (IQA) errors")
        self.GenerateRecoveryEnergies()
        "+++ Write files +++"
        self.writeFiles()

class recovQ00():
    def __init__(self,atoms=None, systemName=None,workingDir=None,targetProp=None, molCharge=0,geomIDs=None,inputDir=None,outputDir=None,threshold=None):
        self.targetProp             = targetProp 
        self.workingDir             = workingDir 
        self.inputDir               = inputDir 
        self.outputDir              = outputDir 
        self.geomIDs                = geomIDs
        self.atoms                  = atoms
        self.systemName             = systemName
        self.threshold              = threshold
        self.defaults               = constants()
        self.atomicProp             = Q00
        self.molCharge              = float(molCharge)
        self.listWfnQ00             = None
        self.atomicQ00              = None
        self.molecularQ00           = None
        self.recoveryQ00            = None
        self.absRecoveryQ00         = None
        self.listFiles              = None
        self.natoms                 = None
        self.filteredDir            = None
        print (f"POLUS: RecovQ00 initialized") 


    def SetInputFileDir(self):
        if self.inputDir   == None:
            self.inputDir   = vars(self.defaults)["default_indir"] 
        else:
            self.inputDir   = self.inputDir

    def SetOutputFileDir(self):
        if self.outputDir  == None:
            self.outputDir  = vars(self.defaults)["default_filtdir"] 
        else:
            self.outputDir  = self.outputDir

        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

    def SetListOfAtoms(self):
        if self.atoms      == None:
            self.atoms      = ReadAtomLabels(self.inputDir)
            if (self.atoms == None):
               RaiseError(message=" missing atoms")
        else:
            self.atoms      = self.atoms

    def SetSystemName(self):
        if self.systemName == None:
            self.systemName = "MyMOL"
        else:
            self.systemName = self.systemName

    def SetTargetProp(self):
        if self.targetProp == None:
            self.targetProp = "q00"
        else:
            self.targetProp = self.targetProp

    def SetWorkingDir(self):
        if self.workingDir == None:
            self.workingDir = os.getcwd()
        else:
            self.workingDir = self.workingDir

    def SetListFiles(self):
        if self.listFiles  == None:
            self.listFiles = GetListInputFiles(self.inputDir)
        else:
            self.listFiles = self.listFiles

    def SetGeomIDs(self):
        if self.geomIDs   == None:
            randomFile     = self.listFiles[0]
            datasetSize    = len(readfile(randomFile,prop=self.atomicProp)[1])
            self.geomIDs   = [i for i in range(datasetSize)]
        else:
            self.geomIDs   = self.geomIDs

    def SetNatoms(self):
        if self.natoms == None:
            self.natoms = len(self.atoms)
        else:
            self.natoms = self.natoms

    def SetMolCharge(self):
        if self.molCharge == None:
            self.molCharge = 0.0
        else:
            self.molCharge = self.molCharge

    def SetThreshold(self):
        if (self.threshold == None):
            self.threshold = 0.001  # Expressed in e
        else: 
            if (isinstance(self.threshold,float) or isinstance(self.threshold,int)):
                self.threshold = self.threshold
            else:
                RaiseWarning(message= " Invalid Threshold parameter. Value set to default (1 me)")
                self.threshold = 0.001

    def CheckUserInputs(self):
        self.SetInputFileDir()
        self.SetOutputFileDir()
        self.SetSystemName()
        self.SetTargetProp()
        self.SetWorkingDir()
        self.SetListOfAtoms()
        self.SetListFiles()
        self.SetGeomIDs()
        self.SetNatoms()
        self.SetThreshold()

    def PrintUserInputs(self): 
        PrintInfo(message=" INDIR   "+self.inputDir)
        PrintInfo(message=" OUTDIR  "+self.outputDir)
        PrintInfo(message=" WORKDIR "+self.workingDir)
        PrintInfo(message=" SYSTEM  "+self.systemName)
        PrintInfo(message=" ATOMS   "+ " ".join([str(x) for x in self.atoms]))

    def SetWfnQ00(self):
        randomFile          = self.listFiles[0]
        with open(randomFile,"r") as file_:
            content = file_.readlines()[1:]
        npoints     = len(content)
        if (self.molCharge==None or not isinstance(self.molCharge,float)):
            self.SetMolCharge()
        if self.listWfnQ00 == None:
            self.listWfnQ00 = list()
            for i in range(npoints):
                if i in self.geomIDs:
                    self.listWfnQ00.append(self.molCharge)

    def ReadAtomicQ00(self):
        randomFile          = self.listFiles[0]
        datasetSize         = len(readfile(randomFile,prop=self.atomicProp)[1])

        if self.atomicQ00      == None:
            self.atomicQ00      = dict()
            count_atoms = 0
            for i in range(self.natoms):
                atom = self.atoms[i]
                self.atomicQ00[atom] = []
                data = readfile(self.listFiles[i],prop=self.atomicProp)
                file_content, q00_idx = data[1], data[4]
                for j in range(len(file_content)):
                    if j in self.geomIDs:
                        line = file_content[j]
                        self.atomicQ00[atom].append(eval(line.split(",")[q00_idx]))
                count_atoms +=1

        return self.atomicQ00

    def SetMolecularQ00(self):
        if self.atomicQ00==None:
            atomicQ00 = self.ReadAtomicQ00()
        else:
            atomicQ00 = self.atomicQ00
        
        "+++ Compute molecular Q00 +++"
        if self.molecularQ00   == None:
            self.molecularQ00   = list() 
            for i in range(len(atomicQ00[self.atoms[0]])):
                mol_q00 = 0
                for key in atomicQ00.keys():
                    mol_q00 += atomicQ00[key][i]
                self.molecularQ00.append(mol_q00)
        else:
            self.molecularQ00   = self.molecularQ00() 
        

    def GenerateRecoveryQ00(self):
        if self.listWfnQ00      == None:
            self.SetWfnQ00()
        if self.molecularQ00 == None:
            self.SetMolecularQ00()
        self.recoveryQ00         = list()  # Defined as (wfn - q00)
        self.absRecoveryQ00      = list()
        if len(self.listWfnQ00)==len(self.molecularQ00):
            for i in range(len(self.listWfnQ00)):
                recov_q00 = self.listWfnQ00[i] - self.molecularQ00[i]
                self.recoveryQ00.append(recov_q00)
                self.absRecoveryQ00.append(abs(recov_q00))
        else:
            RaiseError(message=" number of WFN and Q00 terms not matching")

        
        return self.recoveryQ00, self.absRecoveryQ00


    def writeFiles(self,outputFilename = None):
        geomIDs      = self.geomIDs
        if outputFilename!=None:
            filename = outputFilename
        else:
            filename = os.path.join(self.outputDir,self.systemName.upper()+"_RECOV_Q00")

        if self.recoveryQ00 == None:
            recov_q00, abs_recov_q00 = self.GenerateRecoveryQ00()
        else:
            recov_q00     = self.recoveryQ00
            abs_recov_q00 = [abs(x) for x in recov_q00]
        MOL_Q00     = self.molecularQ00
        MOL_WQ00    = self.listWfnQ00
        if recov_q00 == None:
            RaiseError(message=" Program unable to complete task due to missing WFN Q00")
        with open(filename,"w") as rec_q00_file:
            rec_q00_file.write(f"#FIELDS: Geom-ID Q00-MOL WFN-Q00 RECOV-Q00(e) ABS-RECOV-Q00(e)\n")
            for i in range(len(recov_q00)):
                rec_q00_file.write(f"{str(geomIDs[i]):<7} {MOL_Q00[i]:>12.6f} {MOL_WQ00[i]:>12.6f} {recov_q00[i]:>10.6f} {abs_recov_q00[i]:>10.6f}\n")

        badGeomsFile  = os.path.join(self.outputDir,self.systemName.upper()+"-BAD-GEOMETRIES-Q00")
        goodGeomsFile = os.path.join(self.outputDir,self.systemName.upper()+"-GOOD-GEOMETRIES-Q00")
        filteredDir   = os.path.join(self.outputDir,"FILTERED-BY-Q00")

        countBadGeoms      = 0
        countGoodGeoms     = 0
        goodGeomIDs        = list()
        myBadGeomsFile     = open(badGeomsFile,"w")
        myGoodGeomsFile    = open(goodGeomsFile,"w")
        myBadGeomsFile.write(f"#FIELDS: Geom-ID RECOV-Q00(e)\n")
        myGoodGeomsFile.write(f"#FIELDS: Geom-ID RECOV-Q00(e)\n")
        for i in range(len(recov_q00)):
            if (abs(recov_q00[i])>self.threshold):
                countBadGeoms +=1
                myBadGeomsFile.write(f"{geomIDs[i]:<7} {recov_q00[i]:>10.6f} \n")   
            else:
                countGoodGeoms +=1
                myGoodGeomsFile.write(f"{geomIDs[i]:<7} {recov_q00[i]:>10.6f} \n")   
                goodGeomIDs.append(i)

        print("POLUS: Job completed successfully")
        print(f"POLUS: Total geometries               {len(geomIDs)}")
        print(f"POLUS: Removed geometries             {countBadGeoms}")
        print(f"POLUS: Retained geometries            {countGoodGeoms}")
        if not os.path.isdir(filteredDir):
            os.mkdir(filteredDir)
        if (self.filteredDir==None):
            self.filteredDir = filteredDir
        if (self.listFiles!=None):
            for file in self.listFiles:
                with open(file,"r") as oldFile:
                    content = oldFile.readlines()
                    header = content[0]
                    body = content[1:]
                newFilename = os.path.join(filteredDir,os.path.split(file)[1])
                with open(newFilename,"w") as newFile:
                    newFile.write(header)
                    for i in range(len(body)):
                        if i in goodGeomIDs:
                            newFile.write(body[i])

    def Execute(self):
        "+++ Check user inputs +++"
        self.CheckUserInputs()
        "+++ Print user inputs +++"
        self.PrintUserInputs()
        "+++ Read user inputs +++"
        print(f"POLUS: Setting reference molecular charge")
        self.SetWfnQ00()
        "+++ Read Atomic Q00 +++"
        print(f"POLUS: Reading Atomic charges")
        self.ReadAtomicQ00()
        "+++ Set molecular Q00 +++"
        print(f"POLUS: Computing Molecular charges from atomic Q00")
        self.SetMolecularQ00()
        "+++ Generate Recovery errors +++"
        print(f"POLUS: Computing recovery (Q00) errors")
        self.GenerateRecoveryQ00()
        "+++ Write files +++"
        print(f"POLUS: Writing files")
        self.writeFiles()

        
        


        
        
        
