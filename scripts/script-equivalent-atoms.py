from polus.distributions.analysis import ANALYSER


job = ANALYSER(inputDir="input_files",smoothDistr=False,atoms=None,prop="q00",systemName="BENZENE",outputDir="DISTROS",pValue=0.05,compareMethod="CVM",eqAtomsSM="RS")
#job.Execute()
job.CreateCompositeFiles("BENZENE-ATOM-TYPES.dat")
