import os
import mdtraj as md
#from openbabel import openbabel
from sklearn.decomposition import PCA
from polus.utils.logging import RaiseError



class PCA():
    def __init__(self,filename):
        self.filename        = filename
        self.filenamePDB     = None
        self.trajectory      = None
        self.rotated         = False
        self.atomPairIDs     = None
        self.distanceVectors = None

    def LoadTrajectory(self):
        if not os.path.isfile(self.filename):
            RaiseError(message=" Program cannot find file "+self.filename)
        else:
            if self.trajectory  == None:
                trajObject       = md.formats.XYZTrajectoryFile(file=self.filename,mode="r")
                self.trajectory  = trajObject.read()
                trajObject.close()
            
