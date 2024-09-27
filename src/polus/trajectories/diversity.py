import os
import tracemalloc
import time
import sys
import math
import copy
from polus.trajectories.commons import File
from polus.utils.logging import RaiseError, PrintInfo
from polus.utils.printing import PrintOnTerminal, PrintGBMessage
from polus.trajectories.calculators import ComputeRMSD
import  multiprocessing as mp
from multiprocessing.pool import ThreadPool
import numpy as np



class Sampler(File):
    def __init__(self,filename,ncores=16,printPace=10,performSelection=True,writeFerebusInputs=True,nbatches=None,chunkSize=500,groupAverage=False,weightsVector=None,rotateTraj=False,refGeom=None,seedGeom=None,parallel=False,natoms=None,atoms=None,sampleSize=100,systemName="MOL",outputDir=None,mpSM=None,autoStop=False,threshold=None,rotMethod="KU"):
        super().__init__(filename,natoms,atoms)
        self.divIndexFiles   = None
        self.samplePool      = None
        self.largestSubSample= None
        self.externalSet     = None
        self.selectedGeoms   = None
        self.selectStartTime = None
        self.selectCurrentTime = None
        self.refGeomFilename = refGeom
        self.performSelection  = performSelection
        self.selectPrintPace = printPace
        self.ncores          = ncores
        self.chunkSize       = chunkSize
        self.nbatches        = nbatches
        self.weightsVector   = weightsVector
        self.groupAverage    = groupAverage
        self.writeFBSInputs  = writeFerebusInputs
        self.seedFilename    = seedGeom
        self.tmpcentroid     = None
        if isinstance(sampleSize,int):
            self.sampleSize  = [sampleSize]
        elif isinstance(sampleSize,list) and isinstance(sampleSize[0],int):
            self.sampleSize  = sampleSize
        else:
            RaiseError(message="Invalid Sample Size")
        self.selectedXYZ     = None
        self.outputDir       = outputDir
        self.parallel        = parallel
        self.autoStop        = autoStop
        self.nbatches        = None
        self.smallestRMSDs   = list()
        self.systemName      = systemName
        self.rotateTraj      = rotateTraj
        self.rotMethod       = rotMethod
        if mpSM == None:
            self.mpSM        = "spawn"
        else:
            self.mpSM        = mpSM

        if isinstance(threshold,float):
            self.threshold = [threshold]
        elif threshold==None:
            self.threshold = [0.05]
        else:
            self.threshold = threshold

        tracemalloc.start()

    def UpdateSampleSize(self,condition):
        if condition:
            self.sampleSize.append(len(self.samplePool))
        else:
            self.sampleSize=self.sampleSize

    def SetRMSDMatrix(self):
        print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
        # Rotate Trajectory
        self.RotateTrajectory(self.refGeomFilename,self.rotateTraj,self.rotMethod)
        # Set sample pool
        #if not isinstance(self.rotTraj,dict):
        if not isinstance(self.rotTraj,list):
            RaiseError(message=" Unable to execute diversity-based sampling")
        else:
            #self.samplePool = [x for x in self.rotTraj.keys()]
            self.samplePool = [x for x in range(len(self.rotTraj))]
        # Compute Seed/Centroid
        if self.seedFilename == None:
            self.ComputeCentroid()
        else:
            self.GetSeedGeometry(self.seedFilename)

        # Compute RMSD matrix
        print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
        if self.parallel or self.ngeoms>1000:
            # Checking user-defined chunk size
            if self.chunkSize > self.ngeoms:
                self.chunkSize = self.ngeoms
            # Init RMSD matrix
            self.matrRMSD  = np.zeros(shape=(self.ngeoms,self.ngeoms),dtype=np.float32)
            ngeometries    = self.matrRMSD.shape[0]
            nelements      = ngeometries*ngeometries
            nSubBlocks     = math.ceil(ngeometries/float(self.chunkSize))
            print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
            if self.nbatches == None:
                self.nbatches = math.ceil(ngeometries/float(self.chunkSize))
            print(f"POLUS:{' # Geometries in dataset':<30} {ngeometries:>45}")
            print(f"POLUS:{' # Entries of RMSD matrix':<30} {nelements:>45}")
            print(f"POLUS:{' # Sub-Blocks of RMSD matrix':<30} {nSubBlocks:>45}")
            print(f"POLUS:{' # Batches for computing RMSD matrix':<40} {self.nbatches:>35}")
            print(f"POLUS: Filling in the RMSD matrix. This step can take a very long time!!!")
            print(f"POLUS: Filling in the RMSD matrix. This step can take a very long time!!!")
            print(f"POLUS: Filling in the RMSD matrix. This step can take a very long time!!!")
            # Set list of arguments (indices)
            #self.SetArgsRMSD(kind="Full")
            # Fill in sub-blocks of RMSD Matrix
            lwbound     = 0
            print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
            for n in range(self.nbatches):
                upbound = ngeometries*self.chunkSize*(n+1)
                dim_max = (n+1)*self.chunkSize
                if dim_max<=ngeometries:
                    list_args           = [x for x in range(lwbound,upbound)]
                    self.tempRMSDMatrix = np.zeros(shape=(self.chunkSize*self.ngeoms),dtype=np.float32)
                else:
                    remGeoms            = ngeometries - n*self.chunkSize
                    list_args           = [x for x in range(lwbound,nelements)]
                    self.tempRMSDMatrix = np.zeros(shape=(remGeoms*self.ngeoms),dtype=np.float32)
                print(f"POLUS: Creating new pool of processes for batch # {n+1}")
                print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
                with mp.Pool(processes=self.ncores) as p:
                    results=p.map(func=self.FillRMSDMatrix,iterable=list_args) 
                # Collecting results
                for i, result in enumerate(results,start=0):
                    self.tempRMSDMatrix[i] = result
                lwbound     = upbound
                # Transfer data from block to big matrix
                for j in range(len(list_args)):
                    arg1,arg2 = self.GetArgRMSD(list_args[j])
                    self.matrRMSD[arg1,arg2] = self.tempRMSDMatrix[j] 
                self.tempRMSDMatrix            = None
                print(f"POLUS: Destroying processes")
                print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
        else:
            Ref_Geom=None
            self.ComputeRMSDMatrix(Ref_Geom,self.rotateTraj,self.rotMethod)
        print(f"POLUS: RMSD matrix successfully built !!!")
        print(f"POLUS: RMSD matrix successfully built !!!")
        print(f"POLUS: RMSD matrix successfully built !!!")


    def GenerateSample(self):
        if self.selectedXYZ == None and self.selectedGeoms!=None:
            self.selectedXYZ     = dict()
            #if self.rotateTraj:
            #    for i in range(len(self.selectedGeoms)):
            #        self.selectedXYZ[i] = self.rotTraj[self.selectedGeoms[i]]
            #else:
            #    for i in range(len(self.selectedGeoms)):
            #        self.selectedXYZ[i] = np.array(self.history[self.selectedGeoms[i]])
            for i in range(len(self.selectedGeoms)):
                self.selectedXYZ[i] = np.array(self.history[self.selectedGeoms[i]])

    def WriteRotatedXYZTraj(self,outputFilename=None):
        # Set outputDir
        if self.outputDir==None:
            outDir        = os.getcwd()
        else:
            if not os.path.isdir(self.outputDir):
                os.mkdir(self.outputDir)
            outDir    = self.outputDir
        # Set systemName
        if self.systemName == None:
            self.systemName = "MOL"
        # Set Output Filename
        if outputFilename == None:
            outputFilename = os.path.join(outDir,"ROT-"+self.systemName.upper()+"-TRAJ.xyz")
        # Rotate Trajectory
        self.RotateTrajectory(self.refGeomFilename)
        # Write Output Filename
        with open(outputFilename,"w") as myfile:
            #for key in self.rotTraj.keys():
            for key in range(len(self.rotTraj)):
                myfile.write(str(self.natoms)+"\n")
                myfile.write("GEOM-"+str(key)+"\n")
                geom = self.rotTraj[key].tolist()
                for i in range(self.natoms):
                    atom = geom[i]
                    myfile.write(f"{self.atoms[i][0]:<5} {atom[0]:>12.8f} {atom[1]:>12.8f} {atom[2]:>12.8f} \n")


    def WriteSampleXYZ(self,sample_filename=None):
        if sample_filename==None:
            if self.outputDir==None:
                outDir        = os.getcwd()
            else:
                if not os.path.isdir(self.outputDir):
                    os.mkdir(self.outputDir)
                outDir    = self.outputDir
            sample_filename = os.path.join(outDir,self.systemName.upper()+"-SAMPLE-"+str(len(self.selectedGeoms))+".xyz")
            if os.path.isfile(sample_filename):
                os.remove(sample_filename)
        with open(sample_filename,"w") as myfile:
            for key in self.selectedXYZ.keys():
                myfile.write(str(self.natoms)+"\n")
                myfile.write("GEOM-"+str(key)+"\n")
                geom = self.selectedXYZ[key].tolist()
                for i in range(self.natoms):
                    atom = geom[i]
                    myfile.write(f"{self.atoms[i][0]:<5} {atom[0]:>12.8f} {atom[1]:>12.8f} {atom[2]:>12.8f} \n")

    def WriteDiversityMetrics(self,div_filename=None):
        if div_filename==None:
            if self.outputDir==None:
                outDir        = os.getcwd()
            else:
                if not os.path.isdir(self.outputDir):
                    os.mkdir(self.outputDir)
                outDir    = self.outputDir
            div_filename = os.path.join(outDir,self.systemName.upper()+"-DIVERSITY-"+str(len(self.selectedGeoms))+".dat")
            if os.path.isfile(div_filename):
                os.remove(div_filename)
        with open(div_filename,"w") as myfile:
            myfile.write("#FIELDS GEOM-ID sRMSD(A)\n")
            for i in range(len(self.smallestRMSDs)):
                myfile.write(f"{i:<7} {self.smallestRMSDs[i]:>12.8f} \n")

    def WriteIndexFile(self,index_filename=None):
        if index_filename==None:
            if self.outputDir==None:
                outDir        = os.getcwd()
            else:
                if not os.path.isdir(self.outputDir):
                    os.mkdir(self.outputDir)
                outDir    = self.outputDir
            index_filename = os.path.join(outDir,self.systemName.upper()+"-INDEX-"+str(len(self.selectedGeoms))+".dat")
            if os.path.isfile(index_filename):
                os.remove(index_filename)
        if self.divIndexFiles == None:
            self.divIndexFiles = dict()
        self.divIndexFiles[len(self.selectedGeoms)] = index_filename
        with open(index_filename,"w") as myfile:
            for i in self.selectedGeoms:
                myfile.write(f"{i:<7}\n")

    def WriteCentroid(self,centroid_filename=None):
        if not isinstance(self.centroid,np.ndarray):
            self.ComputeCentroid()
        if centroid_filename==None:
            if self.outputDir==None:
                outDir        = os.getcwd()
            else:
                if not os.path.isdir(self.outputDir):
                    os.mkdir(self.outputDir)
                outDir    = self.outputDir
            centroid_filename = os.path.join(outDir,self.systemName.upper()+"-CENTROID.xyz")
            if os.path.isfile(centroid_filename):
                os.remove(centroid_filename)
        with open(centroid_filename,"w") as myfile:
            myfile.write(f"{self.natoms}\n")
            myfile.write(f"{self.systemName.upper()}-CENTROID\n")
            for i in range(self.natoms):
                symbol = self.GetAtomSymbol(self.atoms[i])
                coords = self.centroid[i]
                myfile.write(f"{symbol:<7} {coords[0]:>12.8f} {coords[1]:>12.8f} {coords[2]:>12.8f} \n")

    def WriteDistanceMatrix(self,outputFilename=None):
        if outputFilename == None:
            if self.outputDir == None:
                outputFilename = os.path.join(os.getcwd(),"DISTANCE-MATRIX.dat")
            else:
                if not os.path.isdir(self.outputDir):
                    os.mkdir(self.outputDir)
                outputFilename = os.path.join(self.outputDir,"DISTANCE-MATRIX.dat")
        if self.matrRMSD == None:
            self.SetRMSDMatrix()
        data = self.matrRMSD.tolist()
        with open(outputFilename,"w") as f:
            f.write("#FIELDS GEOM-ID1 GEOM-ID2 RMSD(angstrom)\n")
            for i in range(self.ngeoms):
                for j in range(self.ngeoms):
                    f.write(f"{i+1:>10} {j+1:>10} {self.matrRMSD[i,j]:>12.6f}\n")


    def GetSmallestRMSD(self,Test_Geom_ID):
        if self.selectedGeoms == None:
            RMSD_array = ComputeRMSD(self.rotTraj[Test_Geom_ID],self.centroid,self.rotMethod,False,self.userWeights)
        else:
            RMSD_array = self.matrRMSD[Test_Geom_ID,self.selectedGeoms]
        return np.min(RMSD_array)

    def SetExternalSet(self):
        if self.externalSet == None:
            self.externalSet = list()
            if self.writeFBSInputs:
                maxTrainSize     = sorted(self.sampleSize,reverse=True)[1]
                self.externalSet = self.selectedGeoms[maxTrainSize:]
            else:
                maxTrainSize     = sorted(self.sampleSize,reverse=True)[0]
                self.externalSet = list()
                for ID in self.samplePool:
                    if ID not in self.selectedGeoms:
                        self.externalSet.append(ID)
        print(f"POLUS: External Set Size {len(self.externalSet)}")

    def SelectGeoms(self):
        args = self.samplePool.copy()
        if not self.groupAverage:
            #args = self.samplePool[1:]
            #if len(self.samplePool)>=20000:
                #nbatches  = math.ceil(len(self.samplePool)/self.ncores)
                #for n in nbatches:
                #    if (n+1)*self.ncores:
                #        args_ = args[n*self.ncores:(n+1)*self.ncores]
                #    else:
                #        args_ = args[n*self.ncores:]
                #    with mp.Pool(processes=self.ncores) as p:
                #    results   = p.map(func=self.GetSmallestRMSD,iterable=args) 
                #    best_rmsd = max(results)
                #    best      = args[results.index(best_rmsd)]

                #with mp.Pool(processes=self.ncores) as p:
                #    results   = p.map(func=self.GetSmallestRMSD,iterable=args) 
                #best_rmsd = max(results)
                #best      = args[results.index(best_rmsd)]
            #else:
            #    best      = self.samplePool[0]
            #    best_rmsd = self.GetSmallestRMSD(best)
            #    for geomID in self.samplePool[1:]: 
            #        smallest_rmsd = self.GetSmallestRMSD(geomID)
            #        if smallest_rmsd > best_rmsd:
            #            best_rmsd = smallest_rmsd
            #            best      = geomID
            best      = self.samplePool[0]
            best_rmsd = self.GetSmallestRMSD(best)
            for geomID in self.samplePool[1:]: 
                smallest_rmsd = self.GetSmallestRMSD(geomID)
                if smallest_rmsd > best_rmsd:
                    best_rmsd = smallest_rmsd
                    best      = geomID
        else:
            #with mp.Pool(processes=self.ncores) as p:
            #    results   = p.map(func=self.GetSmallestRMSD2,iterable=args) 
            #best_rmsd = max(results)
            #best      = args[results.index(best_rmsd)]
            #best      = self.samplePool[0]
            #best_rmsd = self.GetSmallestRMSD2(best)
            #for geomID in self.samplePool[1:]: 
            #    smallest_rmsd = self.GetSmallestRMSD2(geomID)
            #    if smallest_rmsd > best_rmsd:
            #        best_rmsd = smallest_rmsd
            #        best      = geomID
            best,best_rmsd=self.GetNextHitStructure(self.samplePool)   
        if self.selectedGeoms == None:
            self.selectedGeoms = [best]
        else:
            self.selectedGeoms.append(best) 
        self.smallestRMSDs.append(best_rmsd)
        self.samplePool.remove(best)
        # Print selection results
        if len(self.selectedGeoms)==1:
            self.selectCurrentTime = time.time()
            print(f"POLUS: {'N':>7} {'Best':>10} {'Dis(A)':>10} {'Dur(s)':>14}")
            print(f"POLUS: {1:>7} {best:>10} {best_rmsd:>10.3f} {self.selectCurrentTime-self.selectStartTime:>14.4e}")
        else:
            if (len(self.selectedGeoms))%self.selectPrintPace==0:
                self.selectCurrentTime = time.time()
                print(f"POLUS: {len(self.selectedGeoms):>7} {best:>10} {best_rmsd:>10.3f} {self.selectCurrentTime-self.selectStartTime:>14.4e}")

    # GetCentroid method needed too
    def SetCentroid(self,poolIDs):
        if len(poolIDs)==1 and self.tmpcentroid == None:
            self.centroid    =  self.rotTraj[poolIDs[0]]
            self.tmpcentroid =  self.rotTraj[poolIDs[0]]
        else:
            tmp              = np.array(self.tmpcentroid)
            ngeoms           = len(poolIDs)
            expandedTmp      = float(ngeoms-1)*tmp
            newGeom          = np.array(self.rotTraj[poolIDs[-1]])
            self.centroid    = (expandedTmp+newGeom)/float(ngeoms) 
            self.tmpcentroid = self.centroid.copy()
            #ngeoms      = len(poolIDs)
            #centroid    = [[0.0,0.0,0.0] for i in range(self.natoms)]
            #geoms       = list()
            #for ID in poolIDs:
            #    geom=self.rotTraj[ID]
            #    for j in range(self.natoms):
            #        centroid[j] = centroid[j]+geom[j]
            #for k in range(self.natoms):
            #    centroid[k] = centroid[k]/float(ngeoms)
            #self.centroid    = centroid
            #self.tmpcentroid = centroid
        
        
    def GetSmallestRMSD2(self,Test_Geom_ID):   
        if self.selectedGeoms == None:
            minRMSD = ComputeRMSD(self.rotTraj[Test_Geom_ID],self.centroid,self.rotMethod,False,self.userWeights)
        else:
            minRMSD = ComputeRMSD(self.rotTraj[Test_Geom_ID],self.centroid,self.rotMethod,False,self.userWeights)
        return minRMSD

    def GetNextHitStructure(self,Test_Geom_IDs):   
        hit     = Test_Geom_IDs[0]
        dist    = ComputeRMSD(self.rotTraj[hit],self.centroid,self.rotMethod,False,self.userWeights)
        for gid in Test_Geom_IDs[1:]:
            dist_    = ComputeRMSD(self.rotTraj[hit],self.centroid,self.rotMethod,False,self.userWeights)
            if dist_ > dist:
                dist = dist_
                hit  = gid
        #cs      = math.ceil(len(Test_Geom_IDs)/self.ncores)
        #pool    = ThreadPool(self.ncores)
        #results = pool.map(func=self.GetSmallestRMSD2,iterable=Test_Geom_IDs,chunksize=cs)
        #dist    = max(results)
        #hit     = Test_Geom_IDs[results.index(dist)]
        return hit, dist

    def SelectAndWrite(self):
        self.UpdateSampleSize(self.writeFBSInputs)
        print(f"POLUS: Selection of diverse geometries in progress...")
        self.selectStartTime = time.time()
        if not self.autoStop:
            largestSampleSize = max(self.sampleSize)
            for i in range(largestSampleSize):
                if (i+1)<=self.ngeoms:
                    self.SelectGeoms()
                    if self.groupAverage:
                        self.SetCentroid(self.selectedGeoms)
                else:
                    break
            self.SetExternalSet()
            self.GenerateSample()
            self.WriteSampleXYZ() 
            self.WriteDiversityMetrics() 
            self.WriteCentroid()
            self.WriteIndexFile()
            self.largestSubSample = self.selectedGeoms.copy()
            sortedSampleSize  = sorted(self.sampleSize,reverse=True)
            copySelectedGeoms = self.selectedGeoms.copy()
            copySelectedXYZ   = copy.deepcopy(self.selectedXYZ)
            copySmallestRMSDs = self.smallestRMSDs.copy()
            geomKeys          = list(copySelectedXYZ.keys())
            if len(self.sampleSize)>1:
                for sampleSize in sortedSampleSize[1:]:
                    self.selectedGeoms = copySelectedGeoms[:sampleSize]
                    self.selectedXYZ   = {key:copySelectedXYZ[key] for key in geomKeys[:sampleSize]}
                    self.smallestRMSDs = copySmallestRMSDs[:sampleSize]
                    self.WriteSampleXYZ() 
                    self.WriteDiversityMetrics() 
                    self.WriteIndexFile()
        else:
            stop = False
            smallestThreshold = min(self.threshold)
            while not stop:
                self.SelectGeoms()
                if self.smallestRMSDs[-1]<=smallestThreshold or len(self.smallestRMSDs)==self.ngeoms:
                    stop = True
                else:
                    if self.groupAverage:
                        self.SetCentroid(self.selectedGeoms)

            self.SetExternalSet()
            self.GenerateSample()
            self.WriteSampleXYZ() 
            self.WriteDiversityMetrics() 
            self.WriteCentroid()
            self.WriteIndexFile()
            self.largestSubSample = self.selectedGeoms.copy()
            sortedThreshold   = sorted(self.threshold,reverse=False)
            copySmallestRMSDs = self.smallestRMSDs.copy()
            copySelectedGeoms = self.selectedGeoms.copy()
            copySelectedXYZ   = copy.deepcopy(self.selectedXYZ)
            copySmallestRMSDs = self.smallestRMSDs.copy()
            geomKeys          = list(copySelectedXYZ.keys())
            if len(self.threshold)>1:
                for threshold in sortedThreshold[1:]:
                    sampleSize = 0
                    for value in copySmallestRMSDs:
                        if value >= threshold:
                            sampleSize +=1
                        else:
                            break
                    self.selectedGeoms = copySelectedGeoms[:sampleSize]
                    self.selectedXYZ   = {key:copySelectedXYZ[key] for key in geomKeys[:sampleSize]}
                    self.smallestRMSDs = copySmallestRMSDs[:sampleSize]
                    self.WriteSampleXYZ() 
                    self.WriteDiversityMetrics() 
                    self.WriteIndexFile()
        # Print duration
        #PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def Execute(self,Ref_Geom=None):
        # Set Weights
        self.SetWeights(self.weightsVector)
        start_time = time.time()
        if not self.groupAverage:
            # Set RMSD matrix
            self.SetRMSDMatrix()
            end_time1 = time.time()
            print(f"POLUS: Time (s) for building RMSD matrix {end_time1 - start_time:10.6e}")
            # Select Geometries And Write Files
            if self.performSelection:
                self.SelectAndWrite()
                end_time2 = time.time()
                print(f"POLUS: Time (s) for selecting & writing geometries {end_time2 - end_time1:10.6e}")
                print(f"POLUS: Total duration (s) of the DAS procedure {end_time2 - start_time:10.6e}")
            print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")
        else:
            # Rotate Trajectory
            self.RotateTrajectory(self.refGeomFilename,self.rotateTraj,self.rotMethod)
            # Set sample pool
            if not isinstance(self.rotTraj,list):
                RaiseError(message=" Unable to execute diversity-based sampling")
            else:
                self.samplePool = [x for x in range(len(self.rotTraj))]
            # Compute Centroid
            self.ComputeCentroid()
            #RaiseError(message="GroupAverage method not yet implemented")
            # Select Geometries And Write Files
            start_time = time.time()
            self.SelectAndWrite()
            end_time1 = time.time()
            print(f"POLUS: Time (s) for selecting & writing geometries {end_time1 - start_time:10.6e}")
            print(f"POLUS: Total duration (s) of the DAS procedure {end_time1 - start_time:10.6e}")
            print(f"POLUS: Traced memory {tracemalloc.get_traced_memory()}")

