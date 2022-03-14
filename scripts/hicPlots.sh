#!/bin/zsh

## List of hicexplorer commands used to plot hic-matrices
## by: Cullen Roth

## Sourced user's zsh file
source ~/.zshrc

## Activate hic-environment
conda activate hicexplorerenv

## For hic data (in h5 format) taken from hour 00 and 12 of 
## dexamethason exsposure and for chromosomes 6 and 22 plot hic matrix 
for sample in ENCLB571GEP.chr22.200kb.00.h5 ENCLB870JCZ.chr22.200kb.12.h5 ENCLB571GEP.chr6.200kb.00.h5 ENCLB870JCZ.chr6.200kb.12.h5
do
    hicPlotMatrix --matrix ../data/h5/${sample} --outFileName ../figures/h5/${sample}.png \
                  --rotationX 90 --log1p --dpi 150 --colorMap coolwarm \
                  --increaseFigureHeight 1.25 --increaseFigureWidth 1.25 --fontsize 10
done