from polus.trajectories.commons import File
from polus.trajectories.diversity import Sampler
import numpy as np


job = Sampler(systemName="UREA",rotateTraj=True,rotMethod="KU",weightsVector="HL1:2",parallel=True,autoStop=False,threshold=[0.05,0.10,0.025,0.010],outputDir="OUTPUT-UREA",filename="urea.xyz",sampleSize=[500])
job.Execute()
