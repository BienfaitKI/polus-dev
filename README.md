# Overview

POLUS is a data processing package integrated in the workflow of the machine-learned force field FFLUX. The package mainly works on
output from molecular dynamics simulators and electronic structure calculators. It handles a list of input files, including .csv, .txt,
and .xyz files. 

POLUS performs various tasks, such as: 
- Extraction of representative subsets of geometries based on structural diversity.
- Correction of atomic IQA energies following a penalised stakeholder master equation.
- Identification of spatially equivalent atoms based on the distribution of local properties along a dynamic path.
- Stratification of datasets based on the distribution of both atomic and molecular properties.

# Installation

POLUS can be easily installed by following these steps:

1. Download or clone a copy of the copy. For the moment, we have only
made available the development version called "polus-dev".

2. Move inside the downloaded/cloned folder: ```cd polus-dev```

3. ```pip install .```

We recommend that the package be installed in a dedicated virtual environment.
Assuming the name of the environment is "polus", it is easy to create one following
one of the two solutions:

1.  ```python -m venv ~/.venv/polus```

2. ```source ~/.venv/polus/bin/activate```

OR

1. ```conda create --name polus```

2. ```conda activate polus```



# Examples

The package comes with an example folder containing input files and a python script that
tests different functionalities of the code. If successfully executed, the previous script 
generates several folders and files.


