# The hic-converter tools
Description of executalbe tools.

## Converting to .h5 to .hic
"./h5.to.hic.sh" converts input .h5 files to a juicer compatable .hic file.

    conda activate hicexplorerenv
    ./h5.to.hic -m in.h5 -g chr.size.bed

It requires a python environment with HiCExplorer installed. See main page for details. This script utlizes conversion functions from HiCExplorer, custom python scripting, and juicer tools to convert an .h5 file to a .hic file. 

1) in.h5 -> [hicConvertFromat](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html) -> out.ginteractions.tsv

2) in.ginteratcions.tsv -> h5.to.short.py -> out.short.gz

3) in.short.gz -> [juicer pre](https://github.com/aidenlab/juicer/wiki/Pre) -> out.hic
       
## Converting from .summary.txt.gz to .hic

"./summary.to.chrom.hic.py" converts input .summary.txt.gz to an output .hic file. 

    ./summary.to.chrom.hic.py -i input.summary.txt.gz -g genome.id -c chrom.name

As input, the script expects a g-zipped summary.txt file (as seen and stored on the Gene Expression Omnibus), a genome.id (for example mm9, mm10, hg18 or hg19), and a chromosome name. The current version of this script only converts data for a single chromosome. This script was built using a python environment (v3.8.8) with the following python libraries:
    * argparse
    * os
    * pandas
    * numpy
    * gzip
    * subprocess