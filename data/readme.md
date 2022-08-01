# Data used

## Activate and Inactive Chromosome X (.summary.txt.gz files)
Data in .summary.txt.gz format were downloaded from the Gene Expression Omnibus. The first wild-type (WT) replicate from the record GSE99991 with accession record GSM2667262 was used in analysis and development of the *.summary.txt.gz* to *.hic* protocol. Here we store the two haplotype paired files with data on chromosome 13: GSM2667262_WT1.HiC.rep1.cas.chr13.summary.txt.gz and GSM2667262_WT1.HiC.rep1.mus.chr13.summary.txt.gz. The *cas* and *mus* haplotype names correspond to the activate (Xa) and inactive (Xi) chromosome X states. The mm9, Mus musculus reference genome was used in analysis.

## Mouse Embryonic Stem Cells Marks et al. 

### Downloading, aligning, and constructing .h5 data
The following code and instructions were modified from a [HicExplorer tutorial]( https://hicexplorer.readthedocs.io/en/latest/content/mES-HiC_analysis.html) on aligning Hi-C data from female mouse, embryonic stem cells. It was adapted to only analyze data for chromosome 13. Additionally only data from the first replicate of [Marks et al. 2015]( https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0698-x) is used in analysis. See the [main page]( https://github.com/4DGB/hic-converter) of this repository for instruction on setting up a python computing environment.

```
conda activate hicexplorerenv
```

### Make directories for analysis
    mkdir -p genome_mm10 fastq bam bwa GRCh38

### Gather mm10 reference, use wget and tar to unpack the data
    wget http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz -O genome_mm10/chromFa.tar.gz
    tar -xvzf genome_mm10/chromFa.tar.gz

### Index chromosome 13 via bwa
    bwa index -p bwa/chr13_index genome_mm10/chr13.fa

### Find restrction sites across the genome 
Here the SauIII restriction site was used in experiments; it targets and leaves a GATC sites. 
Use the findsites function from HicExplorer.

    hicFindRestSite --fasta ./genome_mm10/chr13.fa -o chr13_SauIII_cut_sites.bed --searchPattern GATC

### Download data from Marks et al. for alignment
    wget SRR1956527_1.fastq.gz -O ./fastq/SRR1956527_1.fastq.gz
    wget SRR1956527_2.fastq.gz -O ./fastq/SRR1956527_2.fastq.gz

### Align data to chromosome 13 via bwa mem
Here we align data seperately for each fastq file. Note we do not sort or filter the output bam files.

    bwa mem -A 1 -B 4 -E 50 -L 0 -t 14 bwa/chr13_index ./fastq/SRR1956527_1.fastq.gz | samtools view -Shb - > ./bam/SRR1956527_1_chr13.bam
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 14 bwa/chr13_index ./fastq/SRR1956527_2.fastq.gz | samtools view -Shb - > ./bam/SRR1956527_2_chr13.bam

### Construct Hi-C contact matrix in .h5 format
Use the build matrix function from HicExplorer.

    hicBuildMatrix --samFiles ./bam/SRR1956527_1_chr13.bam ./bam/SRR1956527_2_chr13.bam \
    --restrictionSequence GATC --danglingSequence GATC --restrictionCutFile chr13_SauIII_cut_sites.bed \
    --threads 5 --inputBufferSize 400000 \
    --outFileName SRR1956527_chr13.h5 --QCfolder ./SRR1956527_chr13_QC 

## A549 Cells From the ENCODE Project
The data listed here was generated by the ENCODE project which has produced paired-end, Hi-C sequencing libraries on A549 cells, a cancer lung cell line. Specifically, the library of the first isogenic replicate ([ENCLB571GEP](https://www.encodeproject.org/experiments/ENCSR662QKG/)) with paired files [ENCFF039FYU](https://www.encodeproject.org/files/ENCFF039FYU/) and [ENCFF479RSE](https://www.encodeproject.org/files/ENCFF479RSE/) for the initial time point and the first isogenic replicate ([ENCLB870JCZ](https://www.encodeproject.org/experiments/ENCSR499RVD/)) with paired files [ENCFF668EDF](https://www.encodeproject.org/files/ENCFF668EDF/) and [ENCFF398SQH](https://www.encodeproject.org/files/ENCFF398SQH/) for the twelve hour time point were used to generate hic files via HiCExplorer. The output .h5 files are not stored here as they are quite large (> 50mb). As an example, we have included the instructions for the construciton of the .h5 file of the twelve hour timepoint below.

### Download the two paired fastq files for the first replicate
    wget https://www.encodeproject.org/files/ENCFF668EDF/@@download/ENCFF668EDF.fastq.gz -O ./fastq/ENCFF668EDF_R1_001.fastq.gz
    wget https://www.encodeproject.org/files/ENCFF398SQH/@@download/ENCFF398SQH.fastq.gz -O ./fastq/ENCFF398SQH_R2_001.fastq.gz

### Download the human reference genome from the ENCODE project 
    wget https://www.encodeproject.org/files/ -O ./GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

### Find restriction sites for MobI in the GRCh38 genome
    hicFindRestSite --fasta ./GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -o GRCh38_MobI_cut_sites.bed --searchPattern GATC

### Index the human genome (this can take awhile but we only need to do it once!)
    bwa index -p bwa/GRCh38_index ./GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta 

### Align data to the human reference genome 
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 14 bwa/GRCh38_index ./fastq/ENCFF668EDF_R1_001.fastq.gz | samtools view -SHb - > ./bam/ENCFF668EDF_R1_001.bam
    bwa mem -A 1 -B 4 -E 50 -L 0 -t 14 bwa/GRCh38_index ./fastq/ENCFF398SQH_R2_001.fastq.gz | samtools view -SHb - > ./bam/ENCFF398SQH_R2_001.bam

### Build .h5 matrix for chromosome 22
    hicBuildMatrix --samFiles ./bam/ENCFF668EDF_R1_001.bam ./bam/ENCFF398SQH_R2_001.bam \
    --region chr22 --threads 4 --inputBufferSize 200000 \
    --restrictionSequence GATC --danglingSequence GATC --restrictionCutFile GRCh38_MboI_cut_sites.bed \
    --QCfolder ENCLB870JCZ.chr22.12 --outFileName ENCLB870JCZ.chr22.12.h5