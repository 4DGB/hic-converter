# hic-converter
A collection of scripts used in preprocessing and setup of Hi-C data for the 4DGB project.

# Installation
### Download and set up
    ## Clone this repo
    git clone git@github.com:4DGB/hic-converter.git

    ## Change directory to scripts
    cd /hic-converter/scripts

    ## Make scripts executable
    chmod +x h5.to.hic.sh
    chmod +x h5.to.short.py

### Generate a conda environment with an installation of HiCExplorer
    conda create -n hicexplorerenv hicexplorer -c bioconda -c conda-forge

# Example Usage
### Converting contact matrix from hicexplorer as an .h5 file to a contact matrix with .hic
    ## Acitvate hic-explorer environment
    conda activate hicexplorerenv

    ## Change directory
    cd /hic-converter/scripts

    ## Envoke h5 to hic conversion
    ./h5.to.hic.sh -m ../data/ENCLB571GEP.chr22.200kb.h5 -g ../sizes/chr22.size.bed -o ../data/ -V true -R false

# Dependencies
### Juicer
[Juicer tools jar file](https://github.com/aidenlab/juicer/wiki/Download) 

Scripts here were developed using version 1.22.01 ([also stored here on this repo](https://github.com/4DGB/hic-converter/tree/main/jar)).

### Python
[Anaconda](https://www.anaconda.com/products/individual) 

A conda environment is generated to install [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html)
