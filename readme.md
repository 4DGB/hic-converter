# hic-converter
A collection of tools, commands, and data used in preprocessing and setup of Hi-C data for the 4D-Genome Browser (4DGB) Project.

## Example usage
### Converting contact matrix from HiCExplorer as an .h5 file to a contact matrix with .hic
    
```
## Activate hic-explorer environment
conda activate hicexplorerenv

## Change directory to tools
cd ./hic-converter/tools

## Convert the .h5 file to .hic formated file for Mus musculus chromosome 13
./h5.to.hic.sh -m ../data/h5/SRR1956527_chr13.h5 -g ../data/sizes/mm10.chr13.size.bed -o ../data/hic/SRR1956527_chr13.200kb.hic
```

### Generating a .hic file from .summary.txt.gz file

```
## Change directory to tools
cd ./hic-converter/tools

## Convert a summary.txt.gz file to an .hic file for a single chromosome
./summary.to.chrom.hic.py -i ../data/summary/GSM2667262_WT1.HiC.rep1.mus.chr13.summary.txt.gz -g mm9 -c chr13 -O ../data/hic/GSM2667262_WT1.HiC.rep1.mus.chr13.200kb.hic
```

### Convert juicer merged_nodups (long format) file for chromosome 22 to .hic

```
## Activate hic-explorer environment
conda activate hicexplorerenv

## Change directory to tools
cd ./hic-converter/tools

## Envoke our long to chrom hic function
./long.to.chrom.hic.py -i ../data/long/merged_nodups.chr22.subsampled.txt.gz -g ../data/sizes/GRCh38.chr22.size.bed -c chr22 -O ../data/hic/chr22.10kb.hic -R 10000
````

## Installation
### Download and set up

```
## Clone this repo
git clone git@github.com:4DGB/hic-converter.git

## Change directory to *tools*
cd ./hic-converter/tools

## Make scripts within the *tools* directory executable
chmod +x *.sh
chmod +x *.py
```

### Generate a conda environment with an installation of HiCExplorer
```
conda create -n hicexplorerenv hicexplorer -c bioconda -c conda-forge
```

#### Other needed python libraries
    os, argparse, numpy, pandas, gzip, subprocess

## Dependencies
### Juicer
[Juicer tools jar file](https://github.com/aidenlab/juicer/wiki/Download) 

Scripts here were developed using version 1.22.01 ([also stored here on this repo](https://github.com/4DGB/hic-converter/tree/main/tools)).

### Python
[Anaconda](https://www.anaconda.com/products/individual) 

A conda environment is generated to install [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html)

## Useful commands for 4DGB project
### Example using juicer tools for Hi-C data preprocessing

```
## change directory and make alias callable in this instance
cd ./hic-converter
shopt -s expand_aliases

## Set the juicer alias
alias juicer='java -Xms512m -Xmx2048m -jar ./tools/juicer_tools_1.22.01.jar'

## Splitting off single chromosome (chromosome 22) data from large Hi-C file that is KR balanced and resolved at 200 kb
juicer pre in.short out.chr22.200kb.hic ./sizes/chr22.size.bed -r 200000 -k KR
```

## Example Hi-C creation and formatting
### Pipeline(s) outline and file outputs
R1.fastq.gz, R2.fastq.gz ->  [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/)  -> output.h5 -> (.h5 to .hic, see above)

or
 
R1.fastq.gz, R2.fastq.gz -> [juicer pipeline](https://github.com/aidenlab/juicer/wiki) -> merged_nodups.txt.gz -> juicer tools (see below)

### Formatting with Juicer Pre
### Using juicer pre to extract contacts along chromosome 1, with 200 kb resolution, and Knight Runiz correction.
    juicer pre -k KR -r 200000 merged_nodups.txt.gz output.chr1.hic chr1.size.bed

### Example data within chr1.size.bed file used above
    chr1	248956422

## Work Cited

- [Marks *et al.*](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0698-x) "Dynamics of gene silencing during X inactivation using allele-specific RNA-seq." *Genome biology* 16.1 (2015): 1-20.

- [Durand *et al.*](https://www.sciencedirect.com/science/article/pii/S2405471216302198) "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." *Cell systems* 3.1 (2016): 95-98.

- [Wang *et al.*](https://www.sciencedirect.com/science/article/pii/S0092867418305841) "SMCHD1 merges chromosome compartments and assists formation of super-structures on the inactive X." *Cell* 174.2 (2018): 406-421.

- [Wolff *et al.*](https://hicexplorer.readthedocs.io/en/latest/index.html) "Galaxy HiCExplorer 3: a web server for reproducible Hi-C, capture Hi-C and single-cell Hi-C data analysis, quality control and visualization." *Nucleic acids research* 48.W1 (2020): W177-W184.

- [Lappala *et al.*](https://www.pnas.org/doi/abs/10.1073/pnas.2107092118) "Four-dimensional chromosome reconstruction elucidates the spatiotemporal reorganization of the mammalian X chromosome." Proceedings of the National Academy of Sciences 118.42 (2021): e2107092118.