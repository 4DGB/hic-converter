#!/bin/bash
## Setup 
## ------------------------------------------------------------------------------------------ ##
## NOTE: Be sure to add the following to your .bashrc. 
## Also, replace the path after -jar with the path to juicer tools jar. 
## 
## alias juicer='java -Xms512m -Xmx2048m -jar /Users/croth/Downloads/juicer_tools_1.22.01.jar'
## 
## Example usage:
## ------------------------------------------------------------------------------------------ ##
## 
## ## Source user's bashrc explicitly
## source ~/.bashrc
## 
## ## Load dependencies
## module load Anaconda3/2020.07
## 
## ## Activate hicexplorer python environment
## conda activate /opt/apps/anaconda/3.2021.05/envs/hicexplorer
##
## ## Call the conversion script
## h5.to.hic.sh -m /panfs/biopan04/4DGENOMESEQ/HIC/VERO/2401_019/2401.019.h5 \
##	-b 1000,10000,100000,1000000,10000000,100000000 \
##	-g /panfs/biopan04/4DGENOMESEQ/REFERENCES/VERO/VERO.genome.size.txt \
##	-o /panfs/biopan04/4DGENOMESEQ/HIC/VERO/2401_019/ -V true -R false
##
## ------------------------------------------------------------------------------------------ ##
## Set default variables and help message
verbose=false
remove=true
current=`pwd`
outdir=${current}/
jarpath='/home/croth/juicer_tools_1.22.01.jar'
helpmessage="\nh5.to.hic.sh [options] -m [h5 matrix] -b [binsize] -g [genome size file]

Inputs (required):
-m The HiC matrix in h5 format to be converted. 
-b The binsizes used in the conversion of h5 to hic. 
-g A genome (or chromosome) size file. This is a tsv file with two columns:
    1) the chromosome(s) or contig(s) name(s)
    2) the length of the chromosome(s) or contig(s)

Options:
-o Output directory to place new HiC matrix in hic format (default: $outdir)
-V Bollean flag to run script in verbose mode, good for debugging (default: $verbose). 
-R Boolean flag to remove tempory files (.toremove.) after completing analysis (default: $remove). 
-J Path to juicer tools jar file (default: $jarpath).

Dependencies include: juicer tools, python3, and HiCExplorer.\n\n"

## Gather Variables
## ------------------------------------------------------------------------------------------ ##
while getopts "m:b:g:o:V:R:J:h" opt; do
    case $opt in
        m) h5=$OPTARG;;
        b) binsize=$OPTARG;;
        g) genomesizes=$OPTARG;;
        o) outdir=$OPTARG;;
        V) verbose=$OPTARG;;
        R) remove=$OPTARG;;
        J) jarpath=$OPTARG;;
        h) printf "$helpmessage"
            exit;;
        ?) printf "$helpmessage"
            exit;;
    esac
done

## Check if any arguments are unknown or missing and print help
if [ $OPTIND -eq 1 ]; then
    printf "$helpmessage"
    exit
fi

## Convert to a hic matrix
## Gather basename
math5=`basename $h5`

## Make alias' callable 
shopt -s expand_aliases

## Add juicer path
alias juicer='java -Xms512m -Xmx2048m -jar $jarpath'

## Print juicer path if verbose
if $verbose ; then
    echo 
    echo "WARNING: Printing juicer pre help message."
	echo `juicer pre --help`
    echo 
fi

## ------------------------------------------------------------------------------------------ ##
## Convert to gene interactions file from h5
echo "    h5 --> tsv  "
echo "Converting h5 matrix to ginteractions (tsv) file via hicConvertFormat."

hicConvertFormat -m $h5 -o ${outdir}${math5}.toremove.ginteractions --inputFormat h5 --outputFormat ginteractions

## Convert tsv to shorted short file
echo " "
echo "   tsv --> short"
echo "Converting ginteractions file to sorted, short file."
h5.to.short.py ${outdir}${math5}.toremove.ginteractions.tsv $genomesizes ${outdir}${math5}.toremove.short $verbose

## Call juicer
echo " "
echo " short --> hic  "
echo "Activating juicer pre and converting short file to hic."

## Gather the number of bins within the variable binsize
bincount=$(($(echo $binsize | grep -o "," | wc -l)+1))

## Test if we have only one binzie
if [[ $bincount -eq 1 ]] ; then

    ## If so make a directory for that bin size
    mkdir -p ${outdir}BP_${binsize}/

    ## Call juicer pre for each chromosome for just this one bin size
    ## iterate thru chromosomes: in future versions we will make this a variable
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
    do
        echo "Converting short file for:" $chr
        juicer pre -r $binsize -c $chr ${outdir}${math5}.toremove.short ${outdir}BP_${binsize}/${math5}.${chr}.hic $genomesizes
        echo "Finished:" $chr
    done

    ## Call juicer for the genome for just this binsize
    juicer pre -r $binsize ${outdir}${math5}.toremove.short ${outdir}BP_${binsize}/${math5}.hic $genomesizes
else
    ## Call juicer for the genome for all the binsizes in binsize
    juicer pre -r $binsize ${outdir}${math5}.toremove.short ${outdir}${math5}.hic $genomesizes

fi 
## ------------------------------------------------------------------------------------------ ##
## Remove the temporary files
if $remove ; then

## If true print so
    if $verbose ; then
        echo " "
        echo "Removing temporary files ..."
    fi 

    ## remove files
    rm ${outdir}${math5}.toremove*
fi 

## Fin
echo
echo "Done :-)"