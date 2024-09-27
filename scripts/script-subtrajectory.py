from polus.trajectories.commons import File


Traj = "PYRAZINAMIDE-TRAJ.xyz"
fileObject = File(filename=Traj)
fileObject.ExtractSubTrajectory(atoms=[1,4],outputFilename="PYRAZINAMIDE-TRAJ-O1-N4.xyz")
