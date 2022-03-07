# hic-converter
A collection of scripts used in preprocessing and setup of Hi-C data for the 4DGB project.

# Example Usage
### Converting contact matrix from hicexplorer as an .h5 file to a contact matrix with .hic
    ## Acitvate hic-explorer environment
    conda activate hicexplorerenv

    ## Change directory
    cd ./hic-converter/scripts

    ## Envoke h5 to hic conversion for chromosome 22
    ./h5.to.hic.sh -m ../data/h5/ENCLB571GEP.chr22.200kb.h5 -g ../sizes/chr22.size.bed -o ../data/hic/

# Installation
### Download and set up
    ## Clone this repo
    git clone git@github.com:4DGB/hic-converter.git

    ## Change directory to scripts
    cd ./hic-converter/scripts

    ## Make scripts executable
    chmod +x h5.to.hic.sh
    chmod +x h5.to.short.py

### Generate a conda environment with an installation of HiCExplorer
    conda create -n hicexplorerenv hicexplorer -c bioconda -c conda-forge

# Dependencies
### Juicer
[Juicer tools jar file](https://github.com/aidenlab/juicer/wiki/Download) 

Scripts here were developed using version 1.22.01 ([also stored here on this repo](https://github.com/4DGB/hic-converter/tree/main/jar)).

### Python
[Anaconda](https://www.anaconda.com/products/individual) 

A conda environment is generated to install [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html)

# Data used
Data used here was taken from the ENCODE project which has produced paired-end sequencing libraries on A549 cells, a cancer lung cell line. Specifically, the paired-end sequencing library of the first isogenic replicate ([ENCLB571GEP](https://www.encodeproject.org/experiments/ENCSR662QKG/)) with paired files [ENCFF039FYU](https://www.encodeproject.org/files/ENCFF039FYU/) and [ENCFF479RSE](https://www.encodeproject.org/files/ENCFF479RSE/) were used to generate hic files via HiCExplorer. 