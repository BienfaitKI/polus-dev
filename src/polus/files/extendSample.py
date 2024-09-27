import os
from polus.utils.logging import RaiseError




class DATA():
    def __init__(self,origDir=None,complDir=None,outputDir=None):
        self.origDir     = origDir
        self.complDir    = complDir
        self.outputDir   = outputDir
        self.origFiles   = None
        self.complFiles  = None
        self.outputFiles = None


    def SetDirs(self):
        Dirs = [self.origDir,self.complDir,self.outputDir]
        for Dir in Dirs:
            if Dir  == None:
                 Dir = os.getcwd()
            else:
                if isinstance(Dir,str):
                    Dir = Dir
                else: 
                    RaiseError(message="Invalid directory path "+Dir)

        if self.origDir == self.outputDir or self.complDir == self.outputDir:
            self.outputDir = os.path.join(os.getcwd(),"OUTPUT")

        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)

    def SetOrigFiles(self):
        if self.origDir == None:
            self.SetDirs()
        csvFiles = list()
        for file in os.listdir(self.origDir):
            if file.endswith(".csv"):
                csvFiles.append(os.path.join(self.origDir,file))
        csvFiles = list(set(csvFiles))
        self.origFiles = csvFiles
        for filePath in self.origFiles:
            if not os.path.isfile(filePath):
                RaiseError(message="Program cannot find file "+filePath)


    def SetComplFiles(self):
        if self.complDir == None:
            self.SetDirs()
        if self.origFiles == None:
            self.SetOrigFiles()
        self.complFiles = list()
        for file_ in self.origFiles:
            tail_ = os.path.split(file_)[1]
            self.complFiles.append(os.path.join(self.complDir,tail_))   
            if not os.path.isfile(self.complFiles[-1]):
                RaiseError(message="Program cannot find file "+self.complFiles[-1])


    def SetOutputFiles(self):
        if self.outputDir == None:
            self.SetDirs()
        if self.origFiles == None:
            self.SetOrigFiles()
        self.outputFiles = list()
        for file_ in self.origFiles:
            tail_ = os.path.split(file_)[1]
            self.outputFiles.append(os.path.join(self.outputDir,tail_))   


    def WriteExtendedFiles(self):
        if self.origFiles == None:
            self.SetOrigFiles()
        if self.complFiles == None:
            self.SetComplFiles()
        if self.outputFiles == None:
            self.SetOutputFiles()

        for i in range(len(self.origFiles)):
            origFile    = self.origFiles[i]
            complFile   = self.complFiles[i]
            outputFile  = self.outputFiles[i]
            with open(origFile,"r") as file1:
                content1 = file1.readlines()
            with open(complFile,"r") as file2:
                content2 = file2.readlines()
            newContent   = content1 + content2[1:]
            with open(outputFile,"w") as outFile:
                for line in newContent:
                    outFile.write(line)

    def Execute(self):
        self.SetDirs()
        self.SetOrigFiles()
        self.SetComplFiles()
        self.SetOutputFiles()
        self.WriteExtendedFiles()


class TRAINDATA(DATA):
    def __init__(self,origDir=None,complDir=None,outputDir=None):
        super().__init__(origDir,complDir,outputDir)
        self.origDir          = origDir
        self.complDir         = complDir
        self.outputDir        = outputDir
        self.origTrainFiles   = None
        self.complTrainFiles  = None
        self.outputTrainFiles = None

    def SetOrigtTrainFiles(self):
        if self.origDir == None:
            self.SetDirs()
        csvFiles = list()
        for file in os.listdir(self.origDir):
            if file.endswith(".csv") and "TRAINING_SET" in file:
                csvFiles.append(os.path.join(self.origDir,file))
        csvFiles = list(set(csvFiles))
        self.origTrainFiles = csvFiles
        for filePath in self.origFiles:
            if not os.path.isfile(filePath):
                RaiseError(message="Program cannot find file "+filePath)


    def SetComplTrainFiles(self):
        if self.complDir == None:
            self.SetDirs()
        if self.origTrainFiles == None:
            self.SetOrigTrainFiles()
        self.complTrainFiles = list()
        for file_ in self.origTrainFiles:
            atom      = os.path.split(file_)[1].split("_")[1]
            keyPhrase = atom.upper()+"_processed"
            for File_ in os.listdir(self.complDir):
                if keyPhrase in File_:
                    self.complTrainFiles.append(os.path.join(self.complDir,testFile))   
            if not os.path.isfile(self.complTrainFiles[-1]):
                RaiseError(message="Program cannot find file "+self.complFiles[-1])


    def SetOutputTrainFiles(self):
        if self.outputDir == None:
            self.SetDirs()
        if self.origTrainFiles == None:
            self.SetOrigTrainFiles()
        self.outputTrainFiles = list()
        for file_ in self.origTrainFiles:
            tail_ = os.path.split(file_)[1]
            self.outputTrainFiles.append(os.path.join(self.outputDir,tail_))   

    def MakeShallowCopies(self):
        self.origFiles   = self.origTrainFiles
        self.complFiles  = self.complTrainFiles
        self.outputFiles = self.outputTrainFiles

    def Execute(self):
        self.SetDirs()
        self.SetOrigTrainFiles()
        self.SetComplTrainFiles()
        self.SetOutputTrainFiles()
        self.MakeShallowCopies()
        self.WriteExtendedFiles()


