import os
import sys
from polus.utils.logging import RaiseError, RaiseWarning
from rdkit import Chem
from rdkit.Chem import AllChem,rdchem,rdDistGeom
from rdkit.Chem.rdmolfiles import MolFromXYZFile, MolFromSmiles
from rdkit.Chem.rdMolAlign import AlignMolConformers


class CONFINDER():
    def __init__(self,xyz=None,SMILES=None,alignConfs=False,forceField="UFF",embedMethod="DG",useETA=True,printETA=True,outputDir=None, forceTol=0.001,numConfs = 20, maxSiters = 2000, maxOiters = 500):
        self.inputFile = xyz
        self.topology  = SMILES
        self.embedMeth = embedMethod
        self.useETA    = useETA
        self.printETA  = printETA
        self.outputDir = outputDir
        self.forceTol  = forceTol
        self.numConfs  = numConfs
        self.maxSiters = maxSiters
        self.maxOiters = maxOiters
        self.forceField= forceField
        self.alignConfs= alignConfs
        self.mol       = None
        self.confs     = None
        self.atoms     = None
        self.symbols   = None
        if self.inputFile == None and self.topology==None:
            RaiseError(message="Confinder exited without performing the requested task ")
        if self.outputDir == None:
            self.outputDir = os.apth.join(os.getcwd(),"CONFINDER-OUTPUT")
            if not os.path.isdir(self.outputDir):
                os.mkdir(self.outputDir)
        if isinstance(self.outputDir,str) and not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

    def getMolObject(self):
        if isinstance(self.topology,str):
            self.mol = MolFromSmiles(self.topology)
            self.mol = Chem.AddHs(self.mol)
        if self.mol == None and isinstance(self.inputFile,str):
            try:
                #self.mol = MolFromXYZFile(self.inputFile)
                self.mol = Chem.MolFromXYZBlock(self.inputFile)
            except BaseException as e:
                print(e)
                sys.exit()
        self.atoms   = self.mol.GetAtoms()
        self.symbols = [a.GetSymbol() for a in self.atoms]
        print(f"POLUS: MOL object created successfully")
    
    def generateConformers(self):
        if self.mol == None:
            self.getMolObject()
        if not self.mol == None:
            print(f"POLUS: Embedding method  {self.embedMeth}")
            if self.embedMeth == "DG":
                self.confs=rdDistGeom.EmbedMultipleConfs(mol=self.mol, 
                                                    numConfs=self.numConfs,
                                                    maxAttempts=self.maxSiters, 
                                                    printExpTorsionAngles=self.printETA,
                                                    useExpTorsionAnglePrefs=self.useETA,    
                                                    forceTol=self.forceTol,
                                                    numThreads= os.cpu_count(),
                                                    enforceChirality=False,
                                                    useRandomCoords=True,
                                                    useBasicKnowledge=True)
                self.numConfs = len(self.confs)
                # Check conformer embedding
                if(len(self.confs)>=1):
                    print(f"POLUS: Number of conformers generated {len(self.confs)}")
                    for i in range(len(self.confs)):
                        # Check force field
                        if self.forceField == "UFF":
                            opt  = AllChem.UFFOptimizeMolecule(self.mol,confId=self.confs[i],maxIters=self.maxOiters)
                        else:
                            RaiseError(message="POLUS: Force field unavailable")
                        # Check optimization success
                        if opt == 0:
                            conf_i  = self.mol.GetConformer(i)
                            #print(f"POLUS: Conformer {i} {conf_i}")
                            
                        else:
                            RaiseWarning(message="POLUS: Optimization of conf {i} failed. Try increasing maxOiters")
                else:
                    RaiseError(message="POLUS: Confinder exited abnormally. Embedding failed. Try with higher maxSiter")
            else:
                RaiseError(message="POLUS: Confinder exited abnormally. Only DG method is available for conf. embedding")

        if self.confs == None:
            RaiseError(message="Confinder exited suddenly as it cannot find the requested conformations")

    def writeConformerGeometries(self):
        # Check confs
        if self.confs == None:
            self.generateConformers()
        # ALign conformers of mol object
        if self.alignConfs:
            AlignMolConformers(self.mol)

        # Write conformers
        all_confs_file = os.path.join(self.outputDir,"CONFORMERS.xyz")
        conf_file      = open(all_confs_file,"w")
        for i in range(self.numConfs):
            file_= os.path.join(self.outputDir,"CONF-"+str(i+1)+".xyz")
            conf = self.mol.GetConformer(i)
            with open(file_,"w") as f:
                n = len(self.symbols)
                f.write(str(n)+"\n\n")
                conf_file.write(str(n)+"\n\n")
                for j in range(n):
                    pos  = conf.GetAtomPosition(j)
                    f.write(f"{self.symbols[j]:<5} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f} \n")
                    conf_file.write(f"{self.symbols[j]:<5} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f} \n")
        conf_file.close()

    def analyseConformers(self):
        print("POLUS: Analysis function not yet implemented")
        
    def Execute(self):
        self.getMolObject()
        self.generateConformers()
        self.writeConformerGeometries()
        
