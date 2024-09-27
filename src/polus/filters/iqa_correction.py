import os
from polus.utils.read_module import readfile
from polus.filters.recovery_error import recovEnergy
from polus.utils.logging import RaiseError
from polus.utils.userInputs import ReadAtomLabels
import sys
import shutil

class iqa_correct():
    def __init__(self, system_name, molEnergyFile=None,inputDir=None,allProps=True,outputDir=None,working_directory=None, geom_IDs=None):
        self.outputDir = outputDir
        self.input_dir = inputDir
        self.allProps  = allProps
        self.system_name = system_name
        if working_directory==None:
            self.wdir = os.getcwd()
        else:
            self.wdir = working_directory
        if self.input_dir == None:
            self.input_dir = os.path.join(os.getcwd(),"input_files") 
        self.geom_IDs = geom_IDs
        self.ngeoms = None
        self.correction_terms = None
        self.corrected_iqa_energies = None
        self.corrected_molecular_iqa_energies = None
        self.integration_errors = None
        self.raw_integration_errors = None
        self.sum_integration_errors = None
        self.wfn_energies = None
        self.atomic_iqa_energies = None
        self.molecular_iqa_energies = None
        self.list_files = None
        self.list_atoms = ReadAtomLabels(self.input_dir)
        self.natoms = len(self.list_atoms)
        self.molEnergyFile = molEnergyFile

    def generate_list_files(self):
        system_name = self.system_name
        list_input_dir = os.listdir(self.input_dir)
        if self.list_files==None:
            self.list_files = list()
            for atom in self.list_atoms:
                found = False
                for file_ in list_input_dir:
                    if atom in file_.split("_"):
                        self.list_files.append(os.path.join(self.input_dir,file_))
                        found = True
                if not found:
                    RaiseError(message=f"POLUS: Input file for atom {atom} not found")

    def get_integration_errors(self):
        self.generate_list_files()
        input_data = readfile(self.list_files[0],prop="integration_error")
        file_body, prop_index = input_data[1], input_data[4]
        if self.geom_IDs==None:
            list_geom_IDs = [i for i in range(len(file_body))]
        else:
            list_geom_IDs = self.geom_IDs
        self.ngeoms = len(list_geom_IDs)

        if self.integration_errors==None:
            self.integration_errors ={}
            self.raw_integration_errors ={}
            for i in range(self.natoms):
                atom = self.list_atoms[i]
                self.integration_errors[atom] = []
                self.raw_integration_errors[atom] = []
                data = readfile(self.list_files[i],prop="integration_error")
                file_content, int_error_idx = data[1], data[4]
                for j in range(len(file_content)):
                    if j in list_geom_IDs:
                        line = file_content[j]
                        self.integration_errors[atom].append(abs(eval(line.split(",")[int_error_idx])))
                        self.raw_integration_errors[atom].append(eval(line.split(",")[int_error_idx]))

        if self.sum_integration_errors==None:
            self.sum_integration_errors = []
            for i in range(self.ngeoms):
                int_err_sum = 0.0
                for j in range(self.natoms):
                    atom = self.list_atoms[j]
                    int_err_sum += self.integration_errors[atom][i]
                self.sum_integration_errors.append(abs(int_err_sum))

    def get_energies(self):
        if self.list_files==None:
            self.generate_list_files()
        input_data = readfile(self.list_files[0],prop="integration_error")
        file_body, prop_index = input_data[1], input_data[4]
        if self.geom_IDs==None:
            self.geom_IDs = [i for i in range(len(file_body))]

        if self.molEnergyFile!=None and os.path.isfile(self.molEnergyFile):
            with open(self.molEnergyFile,"r") as f:
                content = f.readlines()
            if len(content)!=len(self.geom_IDs):
                sys.exit("POLUS: Molecular energies not matching atomic contributions! Job exited.")
            self.wfn_energies = list()
            for line in content:
                energy = eval(line.split()[0].replace("\n",""))
                self.wfn_energies.append(energy)
            task = recovEnergy(self.list_atoms, self.system_name, self.list_files,input_directory=self.input_dir,working_directory=self.wdir, geom_IDs=self.geom_IDs)
            task.list_wfn_energies = self.wfn_energies
            task.get_molecular_iqa_energies()
            task.get_recovery_energies()
            self.molecular_iqa_energies = task.molecular_IQA_energies
            self.atomic_iqa_energies = task.atomic_IQA_energies
            self.recovery_energies = task.recov_energies
        else:
            task = recovEnergy(self.list_atoms, self.system_name, self.list_files,input_directory=self.input_dir,working_directory=self.wdir, geom_IDs=self.geom_IDs)
            task.get_list_wfn_energies()
            task.get_molecular_iqa_energies()
            task.get_recovery_energies()
            self.wfn_energies = task.list_wfn_energies
            self.molecular_iqa_energies = task.molecular_IQA_energies
            self.atomic_iqa_energies = task.atomic_IQA_energies
            self.recovery_energies = task.recov_energies
        
    def get_correction_terms(self):
        # Get energies
        self.get_energies()
        # Get integration errors
        self.get_integration_errors()

        if self.correction_terms==None:
            self.correction_terms = {}
            for i in range(self.natoms):
                atom = self.list_atoms[i]
                self.correction_terms[atom] = []
                for j in range(self.ngeoms):
                    corr = self.recovery_energies[j]*self.atomic_iqa_energies[atom][j]/self.molecular_iqa_energies[j] +\
                           self.recovery_energies[j]*(self.integration_errors[atom][j]-self.sum_integration_errors[j]/float(self.natoms))
                    #corr = self.recovery_energies[j]*(self.integration_errors[atom][j]/self.sum_integration_errors[j])
                    #corr = self.recovery_energies[j]*self.atomic_iqa_energies[atom][j]/self.molecular_iqa_energies[j]
                    self.correction_terms[atom].append(corr)


    def get_corrected_atomic_iqa_energies(self):
        # Get atomic iqa energies
        if self.atomic_iqa_energies==None:
            self.get_energies()
        # Get correction terms
        self.get_correction_terms()

        if self.corrected_iqa_energies == None:
            self.corrected_iqa_energies = {}
            for i in range(self.natoms):
                atom = self.list_atoms[i]
                self.corrected_iqa_energies[atom] = []
                for j in range(self.ngeoms):
                    corr_iqa = self.atomic_iqa_energies[atom][j]+self.correction_terms[atom][j]
                    self.corrected_iqa_energies[atom].append(corr_iqa)

    def get_corrected_molecular_iqa_energies(self):
        pass

    def write_raw_and_corrected_atomic_iqa_energies(self,output_dir=None):
        if output_dir!=None:
            outdir = output_dir
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
        else:
            outdir = os.getcwd()

        # Get data
        self.get_corrected_atomic_iqa_energies()
        
        if  os.path.isdir("corr_ref_data"):
            shutil.rmtree("corr_ref_data")

        for i in range(self.natoms):
            atom = self.list_atoms[i]
            out  = os.path.join(outdir,"raw_corr_iqa")
            if not os.path.isdir(out):
                os.mkdir(out)
            filename = os.path.join(out,self.system_name.upper()+"_"+atom.upper()+".dat")
            with open(filename,"w") as myfile:
                myfile.write("#FIELDS GEOM_ID, RAW_IQA_ENERGY, CORR_IQA_ENERGY, ATOMIC_REC_ENERGY, INTEGRATION_ERROR\n")
                for j in range(self.ngeoms):
                    raw_iqa_energy = self.atomic_iqa_energies[atom][j]
                    corr_iqa_energy = self.corrected_iqa_energies[atom][j]
                    atomic_rec_energy = corr_iqa_energy - raw_iqa_energy
                    integration_error = self.raw_integration_errors[atom][j]
                    if j!=self.ngeoms-1:
                        myfile.write(f"{j:>6}    {raw_iqa_energy:>12.8e}    {corr_iqa_energy:>12.8e}    {atomic_rec_energy:>12.8e}    {integration_error:>12.8e}\n")
                    else:
                        myfile.write(f"{j:>6}    {raw_iqa_energy:>12.8e}    {corr_iqa_energy:>12.8e}    {atomic_rec_energy:>12.8e}    {integration_error:>12.8e}")

    def write_corrected_reference_data(self):
        if self.outputDir!=None:
            outdir = os.path.join(self.outputDir,"corr_ref_data")
        else:
            outdir = os.path.join(os.getcwd(),"corr_ref_data")
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        write_all_prop = self.allProps

        # Get data
        self.get_corrected_atomic_iqa_energies()
        self.generate_list_files()
        
        #if  os.path.isdir("corr_ref_data"):
        #    shutil.rmtree("corr_ref_data")

        #print(self.corrected_iqa_energies)
        # Get header and features
        for i in range(self.natoms):
            atom = self.list_atoms[i]
            input_file  = readfile(self.list_files[i],prop="iqa")
            header            = input_file[0]
            file_body         = input_file[1]
            feats_indices     = input_file[3] 
            prop_index        = input_file[4]
            all_prop_ind      = input_file[5]

            #Setting the header of the output file
            csv_header        = [header.split(',')[x] for x in feats_indices]
            if (write_all_prop):
                for index in all_prop_ind:
                    csv_header.append(header.split(',')[index]) 
            else:
                all_prop_ind = [prop_index]
                csv_header.append(header.split(',')[prop_index])
            csv_header = ",".join(csv_header)
            # Get dataset
            dataset = []
            k = 0

            for j in range(self.ngeoms):
                if j in self.geom_IDs:
                    dataset.append([])
                    line = file_body[j].split(",")
                    for idx in feats_indices:
                        dataset[j].append(line[idx])
                    for idx in all_prop_ind:
                        if idx !=prop_index:
                            dataset[k].append(line[idx])
                        else:
                            iqa_energy = str(self.corrected_iqa_energies[atom][k])
                            dataset[k].append(iqa_energy)
                    if prop_index == all_prop_ind[-1]:
                        dataset[k][-1] = dataset[k][-1]+"\n"
                
                    k +=1
  
            #filename = os.path.join(outdir,"corr_ref_data/"+self.system_name.upper()+"_"+atom.upper()+".csv")
            filename = os.path.join(outdir,os.path.split(self.list_files[i])[1])
            with open(filename,"w") as myfile:
                if write_all_prop:
                    myfile.write(csv_header)
                    for i in range(len(dataset)):
                        geom = ",".join(dataset[i])
                        myfile.write(geom)
                else:
                    myfile.write(csv_header)
                    myfile.write('\n')
                    for i in range(len(dataset)):
                        geom = ",".join(dataset[i])
                        myfile.write(geom)


    def get_corrected_molecular_iqa_energies(self):
        pass

    def write_corrected_iqa_energies(self):
        pass

    def Execute(self):
        pass

        
        


        
        
        
