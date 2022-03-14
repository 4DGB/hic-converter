#!/bin/zsh

## List of commands used to convert h5 files to hic file format for 
## chromosomes 6 and 22 at hours 00 and 12 with 100 nM dexamethasone.
## by: Cullen Roth

## Sourced user's zsh file
source ~/.zshrc

## Activate hic-environment
conda activate hicexplorerenv

## For hic data (in h5 format) taken from hour 00 and 12 of 
## dexamethason exsposure and for chromosomes 6, convert to hic file. 
for sample in ENCLB571GEP.chr6.200kb.00.h5 ENCLB870JCZ.chr6.200kb.12.h5
do
    ./h5.to.hic.sh -m ../data/h5/${sample} -g ../sizes/chr6.size.bed -o ../data/hic/
done

## For hic data (in h5 format) taken from hour 00 and 12 of 
## dexamethason exsposure and for chromosomes 22, convert to hic file. 
for sample in ENCLB571GEP.chr22.200kb.00.h5 ENCLB870JCZ.chr22.200kb.12.h5
do
    ./h5.to.hic.sh -m ../data/h5/${sample} -g ../sizes/chr22.size.bed -o ../data/hic/
done