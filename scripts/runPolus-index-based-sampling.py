from polus.samplers.SEQ.Seq import SeqSampler
from polus.filters.RecoveryManager import IqaFilter, Q00Filter, DualFilter
from polus.filters.excluded import EXCLUDED
from polus.samplers.INDEX.indexSampling import SELECT
from polus.filters.outliers import Odd
from polus.trajectories.diversity import Sampler
from polus.filters.iqa_correction import iqa_correct
import os

TRAIN = [100,250,500,750,1000,1250,1500,1750,2000]
OUT   = "datasets-composite-5000-pZS-iqa-corr-q00-valTest-index-based"

# Outlier removal
outlier_job = Odd(inputDir="input_files",outputDir="OUTLIER-CHECK",prop="iqa", method = "pZS")
outlier_job.Execute()

# IQA correction
iqa_corr_job = iqa_correct(inputDir="OUTLIER-CHECK",allProps=True,system_name="ETL",outputDir=None,working_directory=None,geom_IDs=None)
iqa_corr_job.write_raw_and_corrected_atomic_iqa_energies()
iqa_corr_job.write_corrected_reference_data()

# Recovery test filter
q00_job = Q00Filter(threshold=0.001,systemName="ETL",inputDir="corr_ref_data")
q00_job.Execute()

# Find excluded geometries
exclude_job = EXCLUDED(systemName="ETL",fOrigStart=2,fFiltStart=0,inputDir="input_files",filteredDir=q00_job.filteredDir)
exclude_job.WriteExcludedGeomIDs()

# Diversity-based sampling
div_sampling =Sampler(systemName="ETL",chunkSize=1000,writeFerebusInputs=True,weightsVector="HL1:1",rotateTraj=True,rotMethod="KU",refGeom=None,parallel=True,autoStop=False,outputDir="OUTPUT-ETL",filename="ETL-SAMPLE-5000.xyz",sampleSize=TRAIN.copy())
div_sampling.Execute()


# Index-based selection
if not os.path.isdir(OUT):
    os.mkdir(OUT)
for train_set_size in TRAIN:
    outdir = os.path.join(OUT,"INDEX-SAMPLING-"+str(train_set_size))
    job4 = SELECT(allProp=False,props=["iqa"],valTest=True,externalSet=div_sampling.externalSet,considerExtSet=True,indexFile=div_sampling.divIndexFiles[train_set_size],trainSize=train_set_size,validSize=1000,testSize=1000,systemName="ETL",excludedIndexFile=exclude_job.excludedGeomsFilename,inputDir="input_files",outputDir=outdir,randomSeed=20)
    job4.Execute()
