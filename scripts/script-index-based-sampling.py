from polus.filters.RecoveryManager import DualFilter
from polus.samplers.INDEX.indexSampling import SELECT
from polus.trajectories.diversity import Sampler
from polus.filters.excluded import EXCLUDED

trainSetSizes  = [100,250,500,750,1000,1500,2000,2500,3000,4000,5000,6000]
validSetSizes_ = [100,250,500,1000]
trainSetSizes_ = trainSetSizes.copy()

# Create DualFilter object
job1 = DualFilter(inputDir="input_files",thresholdIqa=1.0,thresholdQ00=0.001,systemName="MAL")
# Perform Filtering
job1.Execute()
# Create Exclusion object
job2 = EXCLUDED(systemName=job1.systemName,inputDir="input_files",filteredDir=job1.filteredDir)
# Perform Exclusion
job2.WriteExcludedGeomIDs()
# Create diversity sampling object
job3 =Sampler(systemName=job1.systemName,writeFerebusInputs=True,weightsVector="HL1:1",rotateTraj=True,rotMethod="KU",refGeom=None,seedGeom="MAL-OPT.xyz",parallel=True,autoStop=False,outputDir="OUTPUT-MAL",filename="MALONDIALDEHYDE-TRAJ.xyz",sampleSize=trainSetSizes)
# Perform diversity based selection
job3.Execute()
# Create index based sampler
for train_set_size in trainSetSizes_:
    for valid_set in validSetSizes_:
        job4 = SELECT(allProp=False,props=["iqa"],externalSet=job3.externalSet,considerExtSet=True,indexFile=job3.divIndexFiles[train_set_size],trainSize=train_set_size,validSize=valid_set,testSize=1000,systemName=job1.systemName,excludedIndexFile=job2.excludedGeomsFilename,inputDir="input_files",outputDir="INDEX-SAMPLING-"+str(train_set_size),randomSeed=20)
    #Perform index based selection
        job4.Execute()
