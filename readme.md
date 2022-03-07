# hic-converter
A collection of scripts used in preprocessing and setup of Hi-C data for the 4DGB project.

# Dependencies
    Juicer tools
    Anaconda (conda)

# Installation
### Generate a conda environment with an installation of hicexplorer
`conda create -n hicexplorerenv hicexplorer -c bioconda -c conda-forge`

# Example Usage
### Activate python environment with installation of hicexplorer
conda activate hicexplorerenv

### Converting contact matrix from hicexplorer as an .h5 file to a contact matrix with .hic
./h5.to.hic.sh -m ../data/ENCLB571GEP.chr22.200kb.h5 -g ../sizes/chr22.size.bed -o ../data/ -V true -R false