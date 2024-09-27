#!/bin/bash 
#SBATCH -p multicore 
#SBATCH -n 2 
#SBATCH --job-name=submit-polus 


env=~/.venv/polus/
if [ -d $env]; then
    # Activate polus environment
    source $env/bin/activate

    # Run script
    python3 runPolus.py
else
    echo " Cannot find polus environment"
fi

