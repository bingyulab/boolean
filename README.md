# Boolean

A Boolean model implementation.

```
salloc -N 1 -n 1 --exclusive
module load lang/Python/3.11.5-GCCcore-13.2.0
module load lang/R/4.4.1-gfbf-2023b
module load math/MATLAB/2024a-r6
```

## Installation

CaSQ is available on PyPI and can be installed using pip:

```bash
pip install casq
pip install rpy2  # R interface for Python (recommended for CellNOpt)
# pip install cellnopt.wrapper # cellnopt (alternative, but less reliable)
pip install caspo joblib networkx numpy pandas pydot pyparsing scikit-learn scipy seaborn clingo #caspo
```

### R Packages Installation (for CellNOpt)
```bash
Rscript Install/install_cellnopt.R
```

```bash
https://github.com/sysbiolux/optPBN.git
https://github.com/gingproc-IIM-CSIC/MEIGO64.git
```

```bash
unzip optPBN_stand_alone_v2.2.3
cd optPBN/optPBN_stand_alone_v2.2.3/
matlab -nodisplay -nosplash -r "run('$HOME/boolean/Install/install_optpbn.m'); exit"
```

```bash
cd ~/MEIGO64/MEIGO
matlab -nodisplay -nosplash -r "\
  addpath('$HOME/MEIGO64/MEIGO');               % add MEIGO64 folder
  addpath(genpath('$HOME/MEIGO64/MEIGO'));
  install_MEIGO;                   % run the installer script
  savepath('~/matlab_startup/pathdef.m');                        % save the updated pathdef
  exit" 
```

```bash
mkdir -p ~/matlab_startup

cd ~/matlab_startup && matlab -nodisplay -nosplash
run('~/matlab_startup/startup.m')
```
### Local Installation

#### Linux
```bash
# Install dependencies
apt-get install -y libtinfo5
apt-get install graphviz libgraphviz-dev pkg-config

# Install Python packages
pip install git+https://github.com/DEAP/deap@master
pip install git+https://github.com/hklarner/pyboolnet@master
pip install pydot igraph cairosvg pygraphviz
```

#### macOS
```bash
# Install dependencies
brew install libtinfo5
brew install graphviz libgraphviz-dev pkg-config

# Install Python packages
pip install git+https://github.com/DEAP/deap@master
pip install git+https://github.com/hklarner/pyboolnet@master
pip install pydot igraph cairosvg pygraphviz
```

### Alternative Usage

For a simpler setup, you can open the project in Google Colab and run all cells.


## Usage
### Data Preparation
```bash
casq -s -c -n data/file.xml data/apoptosis.xml
mv data/apoptosisbnet data/apoptosis.bnet
mv data/apoptosissif data/apoptosis.sif
mv data/apoptosiscsv data/apoptosis.csv

python sbmlqual_to_sif_simple.py data/apoptosis.xml data/apoptosis.sif
```

More info on SBMLqual to SIF conversion: [SBMLqual to SIF Conversion Guide](SBMLqual_to_SIF_Guide.md).

