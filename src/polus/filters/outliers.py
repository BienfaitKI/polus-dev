import os
import numpy as np
from scipy.stats import iqr
from polus.utils.read_module import readfile
from polus.utils.logging import RaiseError




class Odd():
    def __init__(self,inputDir=None,outputDir=None, prop="iqa", method="Z-score"):
        self.inputDir    = inputDir
        self.outputDir   = outputDir
        self.prop        = prop
        self.method      = method
        self.inputFiles  = None
        self.outputFiles = None
        self.outliers    = None
        self.cleanIntErr = None
        self.oddIntErr   = None
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)

    def setFilenames(self):
        if os.path.isdir(self.inputDir):
            self.inputFiles  = list()
            self.outputFiles = list()
            items = os.listdir(self.inputDir)
            for item in items:
                if item.endswith(".csv"):
                    self.inputFiles.append(os.path.join(self.inputDir,item))
                    self.outputFiles.append(os.path.join(self.outputDir,item))
        else:
            RaiseError(message="POLUS: Input directory not found ")

    def identifyOutliers(self):
        self.outliers = list()
        if self.inputFiles == None:
            self.setFilenames()
        # read file content / prop vector
        self.cleanIntErr = list()
        self.oddIntErr   = list()
        for file in self.inputFiles:
            print(f"POLUS: Processing {file}")
            data  = np.array(readfile(filename=file,prop=self.prop)[2])
            data2 = np.array(readfile(filename=file,prop="integration_error")[2])
            # find temp outliers
            if self.method == "Z-score" or self.method == "ZS":
                upper = data.mean() + 3.0*data.std()
                lower = data.mean() - 3.0*data.std()
            if self.method == "eZ-score" or self.method == "eZS":
                upper = data.mean() + 3.5*data.std()
                lower = data.mean() - 3.5*data.std()
            if self.method == "extrZ-score" or self.method == "extrZS":
                upper = data.mean() + 4.0*data.std()
                lower = data.mean() - 4.0*data.std()
            if self.method == "pZ-score" or self.method == "pZS":
                upper = data.mean() + 5.0*data.std()
                lower = data.mean() - 5.0*data.std()
            if self.method == "odd" or self.method == "ODD":
                upper = data.mean() + 9.0*data.std()
                lower = data.mean() - 9.0*data.std()
            if self.method == "Z-score-median" or self.method == "ZSM":
                _,m,_ = np.percentile(data,[25,50,75])
                upper = m + 3.0*data.std()
                lower = m - 3.0*data.std()
            if self.method == "eZ-score-median" or self.method == "eZSM":
                _,m,_ = np.percentile(data,[25,50,75])
                upper = m + 3.5*data.std()
                lower = m - 3.5*data.std()
            if self.method == "extrZ-score-median" or self.method == "extrZSM":
                _,m,_ = np.percentile(data,[25,50,75])
                upper = m + 4.0*data.std()
                lower = m - 4.0*data.std()
            if self.method == "pZ-score-median" or self.method == "pZSM":
                _,m,_ = np.percentile(data,[25,50,75])
                upper = m + 5.0*data.std()
                lower = m - 5.0*data.std()
            if self.method == "odd-median" or self.method == "ODDM":
                _,m,_ = np.percentile(data,[25,50,75])
                upper = m + 9.0*data.std()
                lower = m - 9.0*data.std()
            if self.method == "IQR":
                iqrange = iqr (data)
                q1,m,q3 = np.percentile(data,[25,50,75])
                upper   = q3 + 1.5*iqrange
                lower   = q1 - 1.5*iqrange
            temp  = list()
            print(f"POLUS: {'upper-bound':>12}     {'lower-bound':>12}    {'abs_int_err':>12}    {'decision':>10}    {'geom-id':>10} ")
            print(f"POLUS: {upper:>12.6e}    {lower:>12.6e}    {'NA':>12}    {'NA':>10}    {'NA':>10}")
            for i in range(data.shape[0]):
                int_err = abs(data2[i])
                if data[i]>upper or data[i]<lower:
                    temp.append(i)
                    self.oddIntErr.append(int_err)
                    print(f"POLUS: {data[i]:>12.6e}    {data[i]:>12.6e}    {int_err:>12.6e}    {'ODD':>10}    {i:>10}")
                else:
                    self.cleanIntErr.append(int_err)

            # update self.outliers
            for odd in temp:
                if odd not in self.outliers:
                    self.outliers.append(odd)
        
        return self.outliers

    def writeCleanedFiles(self):
        if self.outliers == None:
            self.identifyOutliers()
        for i in range(len(self.inputFiles)):
            in_ = self.inputFiles[i]
            out = self.outputFiles[i]
            with open(in_,"r") as f:
                content = f.readlines()
            h,b  = content[0], content[1:]
            with open(out,"w") as f_:
                f_.write(h)
                for j in range(len(b)):
                    if j not in self.outliers:
                        f_.write(b[j])

    def writeOddGeometries(self):
        if self.outliers == None:
            self.identifyOutliers()
        file_ = os.path.join(os.getcwd(),"outliers.txt")
        with open(file_,"w") as f:
            for i  in range(len(self.outliers)):
                f.write(str(self.outliers[i])+"\n")
        n = len(self.outliers)
        print(f"POLUS: Summary")
        print(f"POLUS: Number of outliers {n}")
        print(f"POLUS: Min. Integration error (pass)  {min(self.cleanIntErr):>12.8e}  {2625.5*min(self.cleanIntErr):>12.8e}")
        print(f"POLUS: Max. Integration error (pass)  {max(self.cleanIntErr):>12.8e}  {2625.5*max(self.cleanIntErr):>12.8e}")
        if len(self.oddIntErr)>=1:
            print(f"POLUS: Min. Integration error (odd)   {min(self.oddIntErr):>12.8e}  {2625.5*min(self.oddIntErr):>12.8e}")
            print(f"POLUS: Max. Integration error (odd)   {max(self.oddIntErr):>12.8e}  {2625.5*max(self.oddIntErr):>12.8e}")

    def Execute(self):
        self.setFilenames()
        self.identifyOutliers()
        self.writeCleanedFiles()
        self.writeOddGeometries()
