# hic-converter
A collection of tools, scripts, and commands used in preprocessing and setup of Hi-C data for the 4DGB project.

# Example usage
### Converting contact matrix from HiCExplorer as an .h5 file to a contact matrix with .hic
    ## Activate hic-explorer environment
    conda activate hicexplorerenv

    ## Change directory to tools
    cd ./hic-converter/tools

    ## Convert the h5 file to hic formated file for chromosome 22
    ./h5.to.hic.sh -m ../data/h5/ENCLB571GEP.chr22.200kb.00.h5 -g ../sizes/chr22.size.bed -o ../data/hic/

### Generating a .hic file from .summary.txt.gz file
    ## Change directory to tools
    cd ./hic-converter/tools

    ## Convert a summary.txt.gz file to an .hic file for a single chromosome
    ./summary.to.chrom.hic.py -i ../data/summary/GSM2667262_WT1.HiC.rep1.mus.chr13.summary.txt.gz -g mm9 -c chr13 -O ../data/hic/GSM2667262_WT1.HiC.rep1.mus.chr13.hic

# Installation
### Download and set up
    ## Clone this repo
    git clone git@github.com:4DGB/hic-converter.git

    ## Change directory to tools
    cd ./hic-converter/tools

    ## Make scripts executable
    chmod +x *.sh
    chmod +x *.py

### Generate a conda environment with an installation of HiCExplorer
    conda create -n hicexplorerenv hicexplorer -c bioconda -c conda-forge

#### Other needed python libraries
    os, argparse, numpy, pandas, gzip, subprocess

# Dependencies
### Juicer
[Juicer tools jar file](https://github.com/aidenlab/juicer/wiki/Download) 

Scripts here were developed using version 1.22.01 ([also stored here on this repo](https://github.com/4DGB/hic-converter/tree/main/tools)).

### Python
[Anaconda](https://www.anaconda.com/products/individual) 

A conda environment is generated to install [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html)

# Data used
The data used here was generated by the ENCODE project which has produced paired-end, Hi-C sequencing libraries on A549 cells, a cancer lung cell line. Specifically, the library of the first isogenic replicate ([ENCLB571GEP](https://www.encodeproject.org/experiments/ENCSR662QKG/)) with paired files [ENCFF039FYU](https://www.encodeproject.org/files/ENCFF039FYU/) and [ENCFF479RSE](https://www.encodeproject.org/files/ENCFF479RSE/) for the initial time point and the first isogenic replicate ([ENCLB870JCZ](https://www.encodeproject.org/experiments/ENCSR499RVD/)) with paired files [ENCFF668EDF](https://www.encodeproject.org/files/ENCFF668EDF/) and [ENCFF398SQH](https://www.encodeproject.org/files/ENCFF398SQH/) for the twelve hour time point were used to generate hic files via HiCExplorer.

# Useful commands for 4DGB project
### Example using juicer tools for Hi-C data preprocessing
    ## change directory and 
    ## make alias callable in this instance
    cd ./hic-converter
    shopt -s expand_aliases

    ## Set the juicer alias
    alias juicer='java -Xms512m -Xmx2048m -jar ./tools/juicer_tools_1.22.01.jar'

    ## Splitting off single chromosome (chromosome 22) data from large Hi-C file
    ## that is KR balanced and resolved at 200 kb
    juicer pre in.short out.chr22.200kb.hic ./sizes/chr22.size.bed -r 200000 -k KR

# Example Hi-C creation and formatting
### Pipeline(s) outline and file outputs
R1.fastq.gz, R2.fastq.gz ->  [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/)  -> output.h5 -> (.h5 to .hic, see above)

or
 
R1.fastq.gz, R2.fastq.gz -> [juicer pipeline](https://github.com/aidenlab/juicer/wiki) -> merged_nodups.txt.gz -> juicer tools (see below)

### Formatting with Juicer Pre
### Using juicer pre to extract contacts along chromosome 1, with 200 kb resolution, and Knight Runiz correction.
juicer pre -k KR -r 200000 merged_nodups.txt.gz output.chr1.hic chr1.size.bed

### Example data within chr1.size.bed file used above
chr1	248956422
