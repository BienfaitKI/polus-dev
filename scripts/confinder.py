from polus.conformers.search import CONFINDER


job = CONFINDER(SMILES="FC2CCCC(C1CCOC1)C2",alignConfs=False,numConfs=100,outputDir="CONFINDER-OUTPUT")
job.Execute()
