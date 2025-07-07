# Boolean

A Boolean model implementation.

## Installation

CaSQ is available on PyPI and can be installed using pip:

```bash
pip install casq
pip install cellnopt.wrapper # cellnopt
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
casq data/file.xml data/apoptosis.xml
```

