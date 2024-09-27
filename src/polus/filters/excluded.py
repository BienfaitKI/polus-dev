import os
import time
from polus.utils.logging import RaiseError, PrintInfo
from polus.utils.printing import PrintOnTerminal



class EXCLUDED():
    def __init__(self,systemName=None,fOrigStart=0,fFiltStart=0,inputDir=None,outputDir=None,filteredDir=None):
        self.inputDir              = inputDir
        self.outputDir             = outputDir
        self.filteredDir           = filteredDir
        self.systemName            = systemName
        self.fOrigStart            = fOrigStart
        self.fFiltStart            = fFiltStart
        self.atomOfChoice          = None
        self.inputFile             = None
        self.filteredFile          = None
        self.excludedIDs           = None
        self.excludedGeomsFilename = None
        self.supportedElements     = ["C","H","O","N","P","F","S","Cl"]
  
    def SetSystemName(self):
        if self.systemName == None:
            self.systemName = "MOL"
        else:
            if isinstance(self.systemName,str):
                self.systemName = self.systemName
            else:
                RaiseError(message="Program cannot specify system Name ")

    def SetDirs(self):
        if self.inputDir == None:
            self.inputDir = os.getcwd()
        else:
            self.inputDir = os.path.join(os.getcwd(),self.inputDir)

        if self.outputDir == None:
            self.outputDir = os.getcwd()
        else:
            self.outputDir = os.path,join(os.getcwd(),self.outputDir)
            if not os.path.isdir(self.outputDir):
                os.makedirs(self.outputDir)

        if self.filteredDir == None:
            self.filteredDir = os.getcwd()
        else:
            self.filteredDir = os.path.join(os.getcwd(),self.filteredDir)

   
    def SetAtomOfChoice(self):
        if self.inputDir == None:
            self.SetDirs()
        for file in os.listdir(self.inputDir):
            if file.endswith(".csv"):
                randomFile = file
                break
        for entry in file.split("_"):
            if entry[0] in self.supportedElements and len(entry)<5:
                self.atomOfChoice = entry
                break

        #print(f"POLUS: Atom of Choice  {self.atomOfChoice}")

    def GetFilesOfInterest(self):
        if self.inputDir == None or self.filteredDir == None:
            self.SetDirs()
        if self.atomOfChoice == None:
            self.SetAtomOfChoice()
        RawFiles   = os.listdir(self.inputDir)
        FiltFiles  = os.listdir(self.filteredDir)
        for file in RawFiles:
            if self.atomOfChoice in file.split("_"):
                 self.inputFile = os.path.join(self.inputDir,file)
                 break

        for file in FiltFiles:
            if self.atomOfChoice in file.split("_"):
                 self.filteredFile = os.path.join(self.filteredDir,file)
                 break


    def FindExcludedPointIDs(self):
        if self.inputFile == None or self.filteredFile == None:
            self.GetFilesOfInterest()
        with open(self.inputFile,"r") as file1:
            content1 = file1.readlines()[1:]
        with open(self.filteredFile,"r") as file2:
            content2 = file2.readlines()[1:]
        self.excludedIDs = list()
        testGeomIDs      = [x.split(",")[self.fFiltStart:5+self.fFiltStart] for x in content2]
        #msg_    = "Finding excluded geometries"
        #start   = time.time()
        #PrintInfo(message = msg_)
        #PrintOnTerminal(msg = msg_)
        for i in range(len(content1)):
            geomID = content1[i].split(",")[self.fOrigStart:5+self.fOrigStart]
            #print(geomID)
            #if geom in content2:
            if geomID not in testGeomIDs:
                self.excludedIDs.append(i)

        #PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

    def WriteExcludedGeomIDs(self,outputFilename=None):
        msg_    = "Writing excluded geometry indices"
        start   = time.time()
        PrintInfo(message = msg_)
        PrintOnTerminal(msg = msg_)
        if self.excludedIDs == None:
            self.FindExcludedPointIDs()
        if self.outputDir == None:
            self.SetDirs()
        if self.systemName == None:
            self.SetSystemName()
        if outputFilename == None or not isinstance(outputFilename,str):
            outputFilename = os.path.join(self.outputDir,self.systemName.upper()+"-EXCLUDED-INDEX.dat")
        else:
            outputFilename = os.path.join(self.outputDir,outputFilename)
        self.excludedGeomsFilename = outputFilename
        with open(outputFilename,"w") as myfile:
            for entry in self.excludedIDs:
                myfile.write(f"{entry:>7}\n")
        
        PrintOnTerminal(duration = time.time()-start,msgLength=len(msg_))

