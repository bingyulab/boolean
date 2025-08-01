# Boolean

This repository contains tools for analyzing and optimizing Boolean networks using various methods, including CellNOpt, MEIGO, and Caspo.


## Load HPC Modules
To run the code on a High-Performance Computing (HPC) cluster, you need to load the necessary modules. Here are the commands to do so:

```bash
salloc -N 1 -n 1 --exclusive
salloc -N 1 -n 1 --exclusive -C skylake
module load lang/R/4.4.1-gfbf-2023b
module load lang/Python/3.11.5-GCCcore-13.2.0
module load lang/SciPy-bundle/2023.11-gfbf-2023b
source ${HOME}/.virtualenvs/boolean/bin/activate
module load lang/Python/3.11.5-GCCcore-13.2.0
module load math/MATLAB/2024a-r6
```
Since the architectures of the clusters are different, and the partitions in each cluster have different architectures, you may need to adjust the module versions accordingly.

```bash
module load lang/SciPy-bundle/2023.11-gfbf-2023b
# python -m pip install --user virtualenv # For local installation
virtualenv --system-site-packages ${HOME}/.virtualenvs/boolean
```
The --system-site-packages option inherits all packages already installed in your system site-packages directory. To use the environment, just load the module and the environment:
```bash
module load lang/SciPy-bundle/2023.11-gfbf-2023b
source ${HOME}/.virtualenvs/boolean/bin/activate
deactivate
```

## Installation
Please follow the instructions in the [installation.md](doc/installation.md) file for installation.
```bash
pip install -r requirements.txt
Rscript Install/install.R
```

## Data Preparation
```bash
casq -s -c -n data/file.xml data/apoptosis.xml
mv data/apoptosisbnet data/apoptosis.bnet
mv data/apoptosissif data/apoptosis.sif
mv data/apoptosiscsv data/apoptosis.csv

python sbmlqual_to_sif_simple.py data/apoptosis.xml data/apoptosis.sif
```
More info on SBMLqual to SIF conversion: [SBMLqual to SIF Conversion Guide](doc/SBMLqual_to_SIF_Guide.md).

