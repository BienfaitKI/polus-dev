import os
from polus.utils.logging import RaiseError

def ReadXYZFile(Filename,readLabels=False):
    if not os.path.isfile(Filename):
        RaiseError(message=f"File {Filename} not found")
    else:
        XYZ = dict()
        with open(Filename,"r") as f:
            content = f.readlines()
        natoms = 0
        for line in content[2:]:
            if len(line.split())>=4:
                natoms+=1
            else:
                break
        ngeoms = int(len(content)/(natoms+2))
        for i in range(ngeoms):
            geom_ = content[(natoms+2)*i:(natoms+2)*(i+1)]
            geom  = list()
            countAtoms = 0
            for atom in geom_[2:]:
                if not readLabels:
                    geom.append([eval(x) for x in atom.split()[1:4]])
                else:
                    geom.append([x for x in atom.split()[:4]])
                    geom[countAtoms][1:4] = [eval(x) for x in geom[countAtoms][1:4]]
                countAtoms+=1
            XYZ[i] = geom
    return XYZ

