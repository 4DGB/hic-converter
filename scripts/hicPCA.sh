##
## List of commands from hic-explorer 
## used to generate PC1 and PC2 bed files
## for chromosomes 6 and 22.
##
## Activate hicenvironment env
conda activate hicexplorerenv
##
## For chromosome 22 call pc1 and pc2
hicPCA -m ../data/h5/ENCLB571GEP.chr22.200kb.h5 -o ../data/bed/ENCLB571GEP.chr22.200kb.pc1.bed ./data/bed/ENCLB571GEP.chr22.200kb.pc2.bed -f bedgraph -we 1 2
##
## For chromosome 6 call pc1 and pc2
hicPCA -m ../data/h5/ENCLB571GEP.chr6.200kb.h5 -o ../data/bed/ENCLB571GEP.chr6.200kb.pc1.bed ./data/bed/ENCLB571GEP.chr6.200kb.pc2.bed -f bedgraph -we 1 2
## 