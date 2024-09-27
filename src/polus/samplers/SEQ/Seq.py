import os
import random
from  polus.utils.logging import RaiseWarning, RaiseError
from polus.files.outputs import WriteJobDetails
from polus.utils.userInputs import ReadAtomLabels, GetProps




class SeqSampler():
    def __init__(self,inputDir=None,skip=0,valTest=False,allProps=False,fromBottom=False,randomSelect=False,systemName=None,atoms=None,props=None,outputDir=None,trainSize=None,validSize=None,testSize=None):
        self.inputDir     = inputDir
        self.outputDir    = outputDir
        self.trainSize    = trainSize
        self.validSize    = validSize
        self.testSize     = testSize
        self.systemName   = systemName
        self.randomSelect = randomSelect
        self.fromBottom   = fromBottom
        self.atoms        = atoms
        self.props        = props
        self.allProps     = allProps
        self.skip         = skip
        self.valTest      = valTest
        self.inputFiles   = None
        if isinstance(self.trainSize,int):
            self.trainSize = [self.trainSize]
        if isinstance(self.validSize,int):
            self.validSize = [self.validSize]
        if isinstance(self.testSize,int):
            self.testSize = [self.testSize]

    def setListOfAtoms(self):
        if self.atoms      == None:
            self.atoms      = ReadAtomLabels(self.inputDir)
            if (self.atoms == None):
                RaiseError(message="POLUS: missing atoms")
        else:
            self.atoms      = self.atoms

    def setListOfProps(self):
        if self.props      == None:
            self.props      = GetProps(allProps=True,momentRank=4)
        else:
            self.props      = self.props

    def checkSystemName(self):
        if self.systemName == None:
            self.systemName = "MOL"
        else:
            if not isinstance(self.systemName,str):
                self.systemName = "MOL"

    def checkOutputDir(self):
        if self.outputDir == None:
            self.outputDir = os.getcwd()
        else:
            if not os.path.isdir(self.outputDir):
                os.mkdir(self.outputDir)

    def checkInputDir(self):
        if not os.path.isdir(self.inputDir):
            RaiseError(message="POLUS: Cannot find input directory")

    def setInputFilenames(self):
        if self.inputFiles == None:
            self.inputFiles = list()
            for file in os.listdir(self.inputDir):
                if file.endswith(".csv"):
                    self.inputFiles.append(os.path.join(self.inputDir,file))

    def GenerateJobDetailsFile(self):
        WriteJobDetails(InDir=self.inputDir,\
                        OutDir=self.outputDir, \
                        system=self.systemName, \
                        atoms=self.atoms, \
                        props=self.props)

    def generateFerebusInputs(self):
        if self.inputFiles == None:
            self.setInputFilenames()
        largestTrainSet = max(self.trainSize)
        largestValidSet = max(self.validSize)
        largestTestSet  = max(self.testSize)
        for file in self.inputFiles:
            print(f"POLUS: Processing {file}")
            with open(file,"r") as f:
                content = f.readlines()
            header      = content.copy()[0]
            body        = content.copy()[1+self.skip:]
            complem     = [i for i in range(len(body)-largestTrainSet,len(body))]
            for _ in range(100):
                random.shuffle(complem)
            validIds    = complem[:largestValidSet]
            for _ in range(100):
                random.shuffle(complem)
            testIds     = complem[:largestTestSet]
            for size1 in self.trainSize:
                for size2 in self.validSize:
                    for size3 in self.testSize:
                        # Set training data
                        trainData = body.copy()[:size1]
                        if not self.randomSelect:
                            if not self.fromBottom:
                                # Set validation data
                                validData = body.copy()[len(body)-size3-size2:len(body)-size3]
                                # Set test data
                                testData = body.copy()[len(body)-size3:]
                            else:
                                # Set validation data
                                validData = body.copy()[size1:size1+size2]
                                # Set test data
                                testData = body.copy()[size1+size2:size1+size2+size3]
                        else:
                            # Set validation data
                            validData = [body.copy()[x] for x in validIds[:size2]]
                            # Set test data
                            testData  = [body.copy()[x] for x in testIds[:size3]]
                        # Check valTest
                        if self.valTest:
                            for line in testData:
                                validData.append(line)
                            testData = validData.copy()
                        # Write files
                        subDir_        = os.path.join(self.outputDir,"SEQ-"+str(size1)+"-"+str(size2)+"-"+str(size3))
                        if not os.path.isdir(subDir_):
                            os.mkdir(subDir_)
                        # Find feature start
                        f_start       = None
                        c             = 0
                        for entry in header.split(","):
                            if entry[0]=="f" and str(entry[1]) in [str(x) for x in range(0,10)]:
                                f_start = c
                                break
                            c+=1 
                        f_end = f_start + 3*len(self.atoms) - 6
                        # Write files
                        for prop in self.props: 
                            subDir        = os.path.join(subDir_,prop)
                            if not os.path.isdir(subDir):
                                os.mkdir(subDir)
                            atom          = os.path.split(file)[1].split("_")[0]
                            trainFilename = os.path.join(subDir,self.systemName+"_"+atom+"_TRAINING_SET.csv")
                            validFilename = os.path.join(subDir,self.systemName+"_"+atom+"_INT_VALIDATION_SET.csv")
                            testFilename  = os.path.join(subDir,self.systemName+"_"+atom+"_EXT_VALIDATION_SET.csv")
                            newHeader     = list()
                            count         = 0
                            for entry in header.split(",")[f_start:]:
                                if entry[0]=="f" and str(entry[1]) in [str(x) for x in range(0,10)]:
                                    newHeader.append(entry.split("_")[0])
                                else:
                                    if not self.allProps:
                                        if prop in entry or entry==prop:
                                            newHeader.append(entry)
                                    else:
                                        if count>=2:
                                            newHeader.append(entry)
                                count = count + 1
                            newHeader = ",".join(newHeader).replace("\n"," ")
                            newHeader = newHeader + "\n"
                            tempHeader=header.replace("\n","")
                            prop_ID   = tempHeader.split(",").index(prop)
                            with open(trainFilename,"w") as f1:
                                f1.write(newHeader)
                                for i in range(len(trainData)):
                                    line_ = trainData[i].split(",")
                                    if not self.allProps:
                                        line  = list()
                                        for j in range(len(line_)):
                                            if j in range(f_start,f_end) or j==prop_ID:
                                                line.append(line_[j])
                                    else:
                                        line = [x for x in line_[2:]]
                                    line  = ",".join(line).replace("\n","")
                                    line  = line+"\n"
                                    f1.write(line)
                            with open(validFilename,"w") as f2:
                                f2.write(newHeader)
                                for i in range(len(validData)):
                                    line_ = validData[i].split(",")
                                    if not self.allProps:
                                        line  = list()
                                        for j in range(len(line_)):
                                            if j in range(f_start,f_end) or j==prop_ID:
                                                line.append(line_[j])
                                    else:
                                        line = [x for x in line_[2:]]
                                    line  = ",".join(line).replace("\n","")
                                    line  = line+"\n"
                                    f2.write(line)
                            with open(testFilename,"w") as f3:
                                f3.write(newHeader)
                                for i in range(len(testData)):
                                    line_ = testData[i].split(",")
                                    if not self.allProps:
                                        line  = list()
                                        for j in range(len(line_)):
                                            if j in range(f_start,f_end) or j==prop_ID:
                                                line.append(line_[j])
                                    else:
                                        line = [x for x in line_[2:]]
                                    line  = ",".join(line).replace("\n","")
                                    line  = line+"\n"
                                    f3.write(line)

    def Execute(self):
        self.checkSystemName()
        self.checkInputDir()
        self.checkOutputDir()
        self.setInputFilenames()
        self.setListOfAtoms()
        self.setListOfProps()
        self.GenerateJobDetailsFile()
        self.generateFerebusInputs()
