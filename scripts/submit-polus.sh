#!/bin/bash --login
#$ -cwd 
#$ -pe smp.pe 16
#$ -l mem512


source ~/.venv/polus/bin/activate

python script.py


