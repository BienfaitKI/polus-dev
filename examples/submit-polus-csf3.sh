#!/bin/bash --login
#$ -cwd 
#$ -pe smp.pe 8

source ~/.venv/polus_test/bin/activate

export OMP_NUM_THREADS=$SNLOTS

python script.py

