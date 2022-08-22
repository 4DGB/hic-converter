# The hic-converter tools
Description of executalbe tools.

## Conversion of ginteractions file to short file
"gin.to.short.py" converts an input ginteractions file a short file.

    ./gin.to.short.py -i ./path/to/ginteractions.tsv -g ./path/to/chrom.sizes.tsv -o ./output/path/out.short

## Converting an .h5 file to an .hic
"h5.to.hic.sh" converts input .h5 files to a juicer compatable .hic file.

    conda activate hicexplorerenv

    ./h5.to.hic.sh -m ./path/to/in.h5 -g ./path/to/chrom.size.bed -o ./path/to/out.hic

It requires a python environment with HiCExplorer installed (see main page for details). This script utlizes conversion functions from HiCExplorer, custom python scripting, and juicer tools to convert an .h5 file to a .hic file. 

1) in.h5 -> [hicConvertFromat](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html) -> out.ginteractions.tsv

2) in.ginteractions.tsv -> gin.to.short.py -> out.short.gz

3) in.short.gz -> [juicer pre](https://github.com/aidenlab/juicer/wiki/Pre) -> out.hic
       
## Converting from .summary.txt.gz to .hic file

"summary.to.chrom.hic.py" converts input .summary.txt.gz to an output .hic file. 

    ./summary.to.chrom.hic.py -i ./path/to/the.summary.gz.txt -g Genome ID -c chromosome name

As input, the script expects a g-zipped summary.txt file (as seen and stored on the Gene Expression Omnibus), a genome.id (for example mm9, mm10, hg18 or hg19), and a chromosome name. The current version of this script only converts data for a single chromosome. This script was built using a python environment (v3.8.8) with the following python libraries:
    * argparse
    * os
    * pandas
    * numpy
    * gzip
    * subprocess

## Formationg a .hic file from merged_nodups.txt file

"long.to.chrom.hic.py" generates an .hic file from the long format file (i.e. a merged_nodupts.txt) [from the juicer pipeline output](https://github-wiki-see.page/m/aidenlab/juicer/wiki/Pre).

    ./long.to.chrom.hic.py -i ./path/to/merged_nodups.txt.gz -g Genome ID -c chromosome name

## Preparing a .hic file for the 4D-Genome Browser

"hic.prep.sh" brings in data from a multi-resolution hic file and converts it to a .hic file compatable with the 4D-Genome Browser (i.e single resolution with KR matrix correction). 

    conda activate hicexplorerenv

    ./hic.prep.sh -m ./path/to/in.hic -g ./path/to/chrom.size.bed -o ./path/to/out.hic

This function utlizes conversion functions from HiCExplorer, custom python scripting, and juicer tools to format an input .hic file to a .hic file analysis in the 4D-Genome Browser. 

1) in.hic -> [hicConvertFromat](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html) -> out.cool

2) in.cool -> hicConvertFromat -> out.ginteractions.tsv

3) in.ginteractions.tsv -> gin.to.short.py -> out.short.gz

4) in.short.gz -> [juicer pre](https://github.com/aidenlab/juicer/wiki/Pre) -> out.hic