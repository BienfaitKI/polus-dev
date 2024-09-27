from polus.samplers.RS.RSmanager import RSampler
from polus.filters.RecoveryManager import IqaFilter, Q00Filter, DualFilter
from polus.filters.outliers import Odd


job = DualFilter(thresholdIqa=1.0,thresholdQ00=0.001,systemName="MAL")
job.Execute()
job2 = Odd(inputDir=job.filteredDir,outputDir="OUTLIER-CHECK",prop="iqa", method = "pZS")
job2.Execute()
job3 = RSampler(systemName="MAL",props=["iqa"], inputDir=job2.outputDir,trainSize=[100,250,500,750,1000],validSize=[100,250,500], testSize=[500])
job3.Execute()
