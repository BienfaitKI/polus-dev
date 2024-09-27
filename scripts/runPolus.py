from polus.samplers.SEQ.Seq import SeqSampler
from polus.filters.RecoveryManager import IqaFilter, Q00Filter, DualFilter
from polus.filters.outliers import Odd


job = DualFilter(thresholdIqa=1.0,thresholdQ00=0.001,systemName="PZD")
job.Execute()
job2 = Odd(inputDir=job.filteredDir,outputDir="OUTLIER-CHECK",prop="iqa", method = "Z-score")
job2.Execute()
job3=SeqSampler(inputDir=job2.outputDir,randomSelect=False,fromBottom=False,props=["iqa"],systemName="PZD",atoms=None,outputDir="datasets-composite-2000",trainSize=[100,250,500,750,1000],validSize=[100,250,500],testSize=[500])
job3.Execute()
