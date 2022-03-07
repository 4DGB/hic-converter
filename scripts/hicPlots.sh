##
## List of hicexplorer commands used to plot hic-matrices
## by - Cullen Roth
## 
## NOTE: be sure to activate python environment with 
##       an installation of hicExplorer. for e.g.
## 
conda activate hicexplorerenv
## 
##
## For chromosome 6
hicPlotMatrix --matrix ../data/h5/ENCLB571GEP.chr6.200kb.h5 --outFileName ../data/figures/ENCLB571GEP.chr6.png \
    --title "A549 Chromosome 6 Hi-C Contacts" --rotationX 90 --log1p --dpi 300 \
    --colorMap coolwarm --increaseFigureHeight 1.25 --increaseFigureWidth 1.25 --fontsize 6 
## 
## 
## For chromosome 22
hicPlotMatrix --matrix ../data/h5/ENCLB571GEP.chr22.200kb.h5 --outFileName ../data/figures/ENCLB571GEP.chr22.png \
    --title "A549 Chromosome 22 Hi-C Contacts" --rotationX 90 --log1p --dpi 300 \
    --colorMap coolwarm --increaseFigureHeight 1.25 --increaseFigureWidth 1.25 --fontsize 6 
##
## 
## fin