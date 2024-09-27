import os
import math
import sys
import time
import numpy as np
from polus.utils.logging import RaiseError, PrintInfo
from polus.utils.printing import PrintOnTerminal
from polus.trajectories.calculators import RotateGeometry, ComputeRMSD
from polus.trajectories.globals import weightsVect


class File():
    def __init__(self,filename,refGeomFilename=None,natoms=None,atoms=None):
        self.filename         = filename
        self.natoms           = natoms
        self.atoms            = atoms
        self.refGeomFilename  = refGeomFilename
        self.userWeights      = None
        self.weightDef        = None
        self.content          = None
        self.ngeoms           = None
        self.nME              = None
        self.history          = None
        self.jump             = None
        self.rotTraj          = None
        self.matrRMSD         = None
        self.centroid         = None
        self.argsRMSD         = None
        self.argcount         = None
        self.refGeom          = None
        self.rotated          = None
        

    def SetWeights(self,weightDef):
        self.weightDef = weightDef
        if self.userWeights == None:
            if isinstance(weightDef,list):
                if self.atoms == None:
                    self.ProcessFile()
                if len(weightDef)==self.natoms and (isinstance(weightDef[0],int) or isinstance(weightDef[0],float)):
                    self.userWeights = weightDef
                else:
                    RaiseError(message="Program unable to set weights")
            elif isinstance(weightDef,str):
                if "HL" in weightDef:
                    if self.atoms == None:
                        self.ProcessFile()
                    heavy = eval(weightDef.split("L")[1].split(":")[0]) 
                    light = eval(weightDef.split("L")[1].split(":")[1]) 
                    self.userWeights = [1.0]*self.natoms
                    for i in range(self.natoms):
                        if self.atoms[i][0]=="H":
                            self.userWeights[i] = float(light)
                        else:
                            self.userWeights[i] = float(heavy)
            else:
                RaiseError(message="Program unable to set weights")
        weightsVect = self.userWeights
        W           = self.userWeights

    def ProcessFile(self):
        if (self.natoms == None):
            self.natoms = 0
            self.atoms  = list()
            with open(self.filename,"r") as xyz_file:
                content = xyz_file.readlines()
            self.content = content
            for line in content[2:]:
                line = line.split()
                if (len(line)==4):
                   self.natoms +=1
                   self.atoms.append(line[0]+str(self.natoms))
                else:
                    break
            self.ngeoms = int(len(self.content)/(self.natoms+2))
            self.nME    = self.ngeoms*self.ngeoms
            self.jump   = self.natoms + 2

        else:
            if (isinstance(self.natoms,int)):
                self.natoms = self.natoms
            else:
                RaiseError(message=" Program cannot set the number of atoms")

    def SetArgsRMSD(self,kind="Half"):
        count = 0
        self.argcount = 0
        if kind=="Half":
            self.argsRMSD   = list()
            for i in range(self.ngeoms):
                for j in range(i+1,self.ngeoms):
                    self.argsRMSD.append([i,j])
                    count+=1
                if (i+1)%1000==0:
                    print(f"Size {sys.getsizeof(self.argsRMSD)}")
        else:
            self.argsRMSD   = list()
            for i in range(self.ngeoms):
                for j in range(self.ngeoms):
                    self.argsRMSD.append([i,j])
                    count+=1
                if (i+1)%1000==0:
                    print(f"Size {sys.getsizeof(self.argsRMSD)}")
        self.argcount = count

    def GetArgRMSD(self,count):
        temp_count = count + 1
        if temp_count<=self.nME:
            i = math.floor(temp_count/self.ngeoms)
            j = temp_count%self.ngeoms
            if j!=0:
                j = j - 1
            else:
                j = self.ngeoms - 1 
                i = i - 1 
        else:
            sys.exit("Cannot access beyond max Elements")
        return i,j

    def SetHistory(self):
        self.history = dict()
        if (self.content==None):
            self.ProcessFile()
        for i in range(self.ngeoms):
            self.history[i] = list()
            geom            = self.content[i*self.jump:(i+1)*self.jump]
            for j in range(2,self.jump):
                atom = geom[j].split()[1:]
                self.history[i].append([eval(x) for x in atom])

    def RotateTrajectory(self,Ref_Geom=None,rotateTraj=True,rotMethod="KU"):
        print("POLUS: Rotating trajectory")
        if self.history==None:
            self.SetHistory()
        if Ref_Geom==None:
            RefGeom = self.history[0]
        else:
            RefGeom = self.ReadRefGeom(Ref_Geom)
        if self.rotTraj==None and rotateTraj:
            self.rotated    = True
            #self.rotTraj    = dict()
            self.rotTraj    = list()
            for i in range(self.ngeoms):
                Test_Geom       = self.history[i]
                #self.rotTraj[i] = RotateGeometry(RefGeom,Test_Geom,rotMethod)
                self.rotTraj.append(RotateGeometry(RefGeom,Test_Geom,rotMethod))
        else:
            self.rotated    = False
            if not rotateTraj:
                #self.rotTraj    = dict()
                self.rotTraj    = list()
                for i in range(self.ngeoms):
                    #self.rotTraj[i] = np.array(self.history[i])
                    self.rotTraj.append(np.array(self.history[i]))
            #elif isinstance(self.rotTraj,dict):
            elif isinstance(self.rotTraj,list):
                self.rotTraj = self.rotTraj
            else:
                RaiseError(message= " Invalid rotated trajectory")

    def ComputeRMSDMatrix(self,Ref_Geom=None,rotateTraj=True,rotMethod="KU"):
        #msg_         = "Computing RMSD matrix"
        #start        = time.time()
        print("POLUS: Computing RMSD matrix")
        if self.matrRMSD==None:
            if self.rotTraj == None:
                self.RotateTrajectory(Ref_Geom,rotateTraj,rotMethod)
            self.matrRMSD = np.zeros(shape=(self.ngeoms,self.ngeoms))
            #PrintInfo(message = msg_)
            #PrintOnTerminal(msg = msg_)
            for i in range(self.ngeoms):
                self.matrRMSD[i,i] = 0.0
                for j in range(i+1,self.ngeoms):
                    self.matrRMSD[i,j] = ComputeRMSD(self.rotTraj[i],self.rotTraj[j],self.rotated,self.userWeights)
                    self.matrRMSD[j,i] = self.matrRMSD[i,j]
        else:
            if isinstance(self.matrRMSD,np.ndarray):
                self.matrRMSD = self.matrRMSD
            else:
                RaiseError(message= " Invalid RMSD matrix")
        #PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def FillRMSDMatrix(self, argID,rotMethod="K"):
        i,j  = self.GetArgRMSD(argID)
        rmsd = ComputeRMSD(self.rotTraj[i],self.rotTraj[j],rotMethod,self.rotated,self.userWeights)

        return rmsd
     
    def ComputeCentroid(self,Ref_Geom=None):
        print("POLUS: Computing Virtual Traj. Centroid")
        if self.rotTraj==None:
            self.RotateTrajectory(Ref_Geom)
        self.centroid = list()
        # TODO: Replace dict with list
        #atoms    = dict()
        atoms    = list()
        for i in range(self.natoms):
            atoms.append([])
            atoms[i] = [0.0,0.0,0.0]
            for j in range(self.ngeoms):
                atoms[i][0] = atoms[i][0] + self.rotTraj[j][i,0]
                atoms[i][1] = atoms[i][1] + self.rotTraj[j][i,1]
                atoms[i][2] = atoms[i][2] + self.rotTraj[j][i,2]
            atoms[i] = [x/self.ngeoms for x in atoms[i]]
            self.centroid.append(atoms[i])

        self.centroid = np.array(self.centroid)

    def GetSeedGeometry(self,seedFilename=None):
        if seedFilename==None:
            self.ComputeCentroid()
        else:
            with open(seedFilename,"r") as f:
                content = f.readlines()[2:]
            self.seed = list()
            for i in range(len(content)):
                self.seed.append([eval(x) for x in content[i].split()[1:]])
             
            self.centroid = self.seed.copy()
            self.centroid = np.array(self.centroid)

    def GetAtomSymbol(self,label):
        symbol = ""
        for char in label:
            if ord(char.upper())>=65 and ord(char.upper())<=90:
                symbol = symbol+char.upper()
            else:
                break
        return symbol

    def ReadRefGeom(self,filename):
        msg_         = "Reading reference geometry"
        start        = time.time()
        PrintInfo(message = msg_)
        PrintOnTerminal(msg = msg_)
        if isinstance(filename,str) and os.path.isfile(filename):
            if self.refGeom == None:
                self.refGeom = list()
                with open(filename,"r") as myfile:
                    content = myfile.readlines()
                for atom in content[2:]:
                    atom = atom.split()[1:]
                    if len(atom)==3:
                        self.refGeom.append([eval(x) for x in atom])
            else:
                if isinstance(self.refGeom,list) and isinstance(self.refGeom[0][0],float):
                    self.refGeom = self.refGeom
                else:
                    RaiseError(message="Invalid Reference Geometry")
        else:
            RaiseWarning(message=f"Program unable to find file {filename}")
            
        PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))
        return self.refGeom


    def RotateTrajectoryThroughALF(self,xyzTraj=None,alf=None,systemName=None,centerAtom=None,listAtoms=None,outputFilename=None,outputDir=None):
        from ichor.core.files import XYZ, Trajectory
        from ichor.core.calculators import default_alf_calculator, calculate_alf_features,calculate_alf_cahn_ingold_prelog
        
        # Set SystemName
        if systemName==None:
            systemName="MOL"
        else:
            systemName=systemName
        # Set XYZ Trajectory
        if xyzTraj==None:
            xyzTraj=self.filename
        # Set Atom List
        if isinstance(listAtoms,list):
            atoms = listAtoms
        else:
            atoms = self.GetSortedAtomList() 
        # Set Center Atom
        if centerAtom==None:
            centerAtom = atoms[0]+str(1)
        # Set Output Directory
        if outputDir==None:
            outDir   = os.getcwd()
        else:
            if isinstance(outputDir,str):
                outDir   = outputDir
            else:
                RaiseError(message="Invalid output directory")
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        # Set Output Filename
        if isinstance(outputFilename,str):
            outFilename = os.path.join(outDir,outputFilename)
        else:
            tail        = systemName.upper()+"-ROTATED-TRAJECTORY-ALF-"+centerAtom.upper()+".xyz"
            outFilename = os.path.join(outDir,tail)
        # Perform rotation
        traj = Trajectory(outFilename)
        xyzObject  = XYZ(xyzTraj)
        if alf==None:
            alf        = xyzObject[centerAtom].alf(default_alf_calculator) 
        else:
            alf        = alf
        # Get ALF atom IDs
        alfAtoms = self.GetAlfIDs(alf)
        atoms    = self.NewListAtoms(atoms,alfAtoms)
        #alf       = xyzObject[centerAtom].alf(calculate_alf_cahn_ingold_prelog) 
        Traj       = Trajectory(xyzTraj)
        for T in Traj:
            feat = calculate_alf_features(T[centerAtom],alf)
            atomsInstance = Trajectory.np_array_to_trajectory(feat,f"dummy",atoms)[0]
            traj.add(atomsInstance)

        traj.write()

    def GetAlfIDs(self,alf):
        return [int(atom) for atom in alf]

    def NewListAtoms(self,atoms,alfIDs):
        listAtoms  = list()
        for atomID in alfIDs:
            listAtoms.append(atoms[atomID])
        for i in range(len(atoms)):
            if i not in alfIDs:
                listAtoms.append(atoms[i])
        return listAtoms

    def GetSortedAtomList(self):
        if self.atoms==None:
            self.ProcessFile()
        atoms = []
        count = 0
        for atom in self.atoms:
            atoms.append("")
            for char in atom:
                if char.isalpha():
                    atoms[count]+=char 
            count+=1
        return atoms
   
    def ExtractSubTrajectory(self,atoms=None,geomRange=None,outputFilename=None):
        if outputFilename == None:
            outputFilename = os.path.join(os.getcwd(),"SUBTRAJECTORY.xyz")
        if self.history   == None:
            self.SetHistory()
        if atoms          == None:
            atoms          = [1]
        f = open(outputFilename,"w")
        if geomRange      == None:
            for i in range(self.ngeoms):
                geom = self.history[i]
                f.write(f"{len(atoms)}\n\n")
                for j in atoms:
                    f.write(f"{self.atoms[j-1]:<6} {geom[j-1][0]:12.6f} {geom[j-1][1]:>12.6f} {geom[j-1][2]:>12.6f}\n")

        f.close()
