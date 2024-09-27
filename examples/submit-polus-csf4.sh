#!/bin/bash 
#SBATCH -p multicore 
#SBATCH -n 2 
#SBATCH --job-name=submit-polus 


env=~/.venv/polus_csf4/
if [ -d $env]; then
    # Activate polus environment
    source $env/bin/activate

    # Run script
    python3 script.py
else
    echo " Cannot find polus environment"
fi

