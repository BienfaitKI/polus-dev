import os
from polus.filters.recoveryError import recovEnergy, recovQ00



class IqaFilter():
    def __init__(self,atoms=None,threshold=None,systemName=None,workingDir=None,targetProp=None,geomIDs=None,inputDir=None,outputDir=None):
        self.atoms       = atoms
        self.systemName  = systemName
        self.workingDir  = workingDir
        self.targetProp  = targetProp
        self.inputDir    = inputDir
        self.outputDir   = outputDir
        self.geomIDs     = geomIDs
        self.dual        = False
        self.threshold   = threshold
        self.filteredDir = None

    def Execute(self):
       job = recovEnergy(atoms=self.atoms,
                         systemName=self.systemName,
                         workingDir=self.workingDir,
                         geomIDs=self.geomIDs,
                         inputDir=self.inputDir,
                         outputDir=self.outputDir, 
                         dualFilter=self.dual,
                         threshold=self.threshold)
       job.Execute()
       if (self.filteredDir==None):
           self.filteredDir = job.filteredDir

class Q00Filter():
    def __init__(self,atoms=None,threshold=None,systemName=None,workingDir=None,targetProp=None,geomIDs=None,inputDir=None,outputDir=None):
        self.atoms       = atoms
        self.systemName  = systemName
        self.workingDir  = workingDir
        self.targetProp  = targetProp
        self.inputDir    = inputDir
        self.outputDir   = outputDir
        self.geomIDs     = geomIDs
        self.threshold   = threshold
        self.filteredDir = None

    def Execute(self):
       job = recovQ00(atoms=self.atoms,
                      systemName=self.systemName,
                      workingDir=self.workingDir,
                      geomIDs=self.geomIDs,
                      inputDir=self.inputDir,
                      outputDir=self.outputDir, 
                      threshold=self.threshold)
       job.Execute()
       if (self.filteredDir==None):
           self.filteredDir = job.filteredDir


class DualFilter():
    def __init__(self,atoms=None,thresholdIqa=None,thresholdQ00=None,systemName=None,workingDir=None,targetProp=None,geomIDs=None,inputDir=None,outputDir=None):
        self.atoms        = atoms
        self.systemName   = systemName
        self.workingDir   = workingDir
        self.targetProp   = targetProp
        self.inputDir     = inputDir
        self.outputDir    = outputDir
        self.geomIDs      = geomIDs
        self.thresholdIqa = thresholdIqa
        self.thresholdQ00 = thresholdQ00
        self.dual         = True
        self.filteredDir  = None

    def Execute(self):
       jobq00 = recovQ00(atoms=self.atoms,
                         systemName=self.systemName,
                         workingDir=self.workingDir,
                         geomIDs=self.geomIDs,
                         inputDir=self.inputDir,
                         outputDir=self.outputDir, 
                         threshold=self.thresholdQ00)
       jobq00.Execute()
       "+++ Update dual filteredDir +++"
       if (self.filteredDir==None):
           self.filteredDir = jobq00.filteredDir
           self.inputDir    = self.filteredDir

       jobiqa = recovEnergy(atoms=self.atoms,
                         systemName=self.systemName,
                         workingDir=self.workingDir,
                         geomIDs=self.geomIDs,
                         inputDir=self.inputDir,
                         outputDir=self.outputDir, 
                         dualFilter=self.dual,
                         threshold=self.thresholdIqa)
       jobiqa.Execute()
       "+++ Update dual filteredDir +++"
       if (self.filteredDir==None or self.filteredDir==jobq00.filteredDir):
           self.filteredDir = jobiqa.filteredDir


