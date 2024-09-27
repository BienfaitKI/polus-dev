#!/bin/bash --login
#$ -cwd 
#$ -pe smp.pe 16
#$ -l mem512

source ~/.venv/polus/bin/activate

export OMP_NUM_THREADS=$SNLOTS

python script-index-based-sampling.py

