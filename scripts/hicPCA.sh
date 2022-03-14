#!/bin/zsh

## List of commands from hic-explorer 
## used to generate PC1 and PC2 bed files
## for chromosomes 6 and 22 at hours 00 and 12 with 100 nM dexamethasone.
## by: Cullen Roth

## Sourced user's zsh file
source ~/.zshrc

## Activate hicenvironment env
conda activate hicexplorerenv

## For hic data (in h5 format) taken from hour 00 and 12 of dexamethason exsposure 
## and for chromosomes 6 and 22 call pc1 and pc2 vi hic-explorer
for sample in ENCLB571GEP.chr22.200kb.00.h5 ENCLB870JCZ.chr22.200kb.12.h5 ENCLB571GEP.chr6.200kb.00.h5 ENCLB870JCZ.chr6.200kb.12.h5
do
    hicPCA -m ../data/h5/${sample} -o ../data/bed/${sample}.pc1.bed ../data/bed/${sample}.pc2.bed -f bedgraph -we 1 2
done