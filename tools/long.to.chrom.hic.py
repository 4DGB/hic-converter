#!/usr/bin/env python
# coding: utf-8

## merged_nodups.txt.gz file to .hic protocol (for a single chromosome).
## Developed by: Cullen Roth, PhD, B-GEN Group.
## For any issues please email: croth@lanl.gov

# In[0]:
## Set default variables, the level of verbosity, the hic tolerance, out put dir, the path to juicer tools jar
tolerance, jarpath, quality = int(0.5*(10**6)), './juicer_tools_1.22.01.jar', 10

## Set the resolutions, the genomeid genome name, and the correction method
resolutions, correction, outpath = '50000,100000,150000,200000,250000,300000,350000,400000', 'KR', './'

## Import argparse and set parser
import argparse

## Make the parse
parser = argparse.ArgumentParser()

## remove actions groups
parser._action_groups.pop()

## Set new actions groups
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

## Set required variables, 
required.add_argument("-i", type=str, required=True, help="Path to input merged_nodups.txt.gz file.\nThis script assumes the file is g-zipped.\nMust provide a relative path.", metavar='./path/to/merged_nodups.txt.gz')
required.add_argument("-g", type=str, required=True, help="Genome ID (for example: hg18, hg19, hg38, mm9, or mm10) or path to chrom.size bed file (chrom\tlength) used in analysis.")
required.add_argument("-c", type=str, required=True, help="Name of the chromosome to parse Hi-C contacts on.")

## Set option variables
optional.add_argument("-R", type=str, required=False, help="Comma-seperated list of resolutions (in bp) for Hi-C map construction (default: %s)."%resolutions, default=resolutions)
optional.add_argument("-K", type=str, required=False, help="Calculate specific matrix normalization and correction.\nList of built-in normalizations: VC, VC_SQRT, KR, and SCALE (default: %s)."%correction, default=correction)
optional.add_argument("-O", type=str, required=False, help="Output path and name of Hi-C file generated from this script.\nDefault is the same as the input path, replacing the .summary.txt.gz extension with .hic",default=None)
optional.add_argument("-J", type=str, required=False, help="Path to juicer tools jar file (default: %s)."%jarpath, default=jarpath)
optional.add_argument("-T", type=int, required=False, help="Lower threshold of Hi-C contacts to flag a potential error (default: %s)."%tolerance, default=tolerance)
optional.add_argument("-Q", type=int, required=False, help="Mapping quality score for filtering Hi-C contacts (default: %s)."%quality, default=quality)
optional.add_argument("-N", type=str, required=False, help="The name used to rename the chromosome within the short file, for example: chr_6a to chr6.\nWARNING: The new name must be in the genome ID or chrom.size file (see -g).", default=None)

## Set optional, boolean variables
optional.add_argument("-V", help="Flag to run in verbose mode.\nDefault behavior is false.", action='store_true')
optional.add_argument("-S", help="Flag to keep intermitant .short file generated during conversion.\nDefault behavior is to remove the .short file", action='store_false')
optional.add_argument("-C", help="Flag to save out the number of Hi-C contacts to file.\nDefault behavior is to skip this step.", action='store_true')

## Set the paresed values
values = parser.parse_args()

# In[1]:
## Set requried variables
inputpath, genomeid, coi = values.i, values.g, values.c

## Set optional variables
resolutions, correction, outpath, jarpath, tolerance, quality, newname = values.R, values.K, values.O, values.J, values.T, values.Q, values.N

## Set boolean variables
verbose, toclean, tocount = values.V, values.S, values.C

# In[2]:
## Check input arugments, bring in os
import os

## Set the list of available referances
referances = ['hg18','hg19','hg38','dMel','mm9','mm10','anasPlat1','bTaurus3','canFam3','equCab2','galGal4','Pf3D8','sacCer3','sCerS288c','susScr3','TAIR10']

## Set the referror message
referror = 'Must be one of: %s.\nAlternatively, this can be the path of the chrom.sizes file that lists on each line the name and size of the chromosomes.'%', '.join(referances)

## Check the genome id is in the list of referances
assert (genomeid in referances) or os.path.exists(genomeid), "ERROR: The genome %s is not available!\n%s"%(genomeid,referror)

## Set the list of matrix corrections
corrections = ['VC', 'VC_SQRT', 'KR', 'SCALE']

## Check that the specified correction method is within this list
assert (correction in corrections), "ERROR: Unknown specified matrix correction.\nSelect from: %s!"%', '.join(corrections)

# In[3]:
## Check that the input path existis
assert os.path.exists(inputpath), "ERROR: The input file (%s) does not exist!"%(inputpath)

## Check that the input file is gzipped
assert (inputpath[-3:] == '.gz'), "ERROR: The input file (%s) is not g-zipped!\nPlease gzip and try again!"%inputpath

## Set the out name and path of the hic file
if outpath is None: ## If the default was given
    hic_path = inputpath.split('.txt')[0]+'.%s.hic'%coi

else: ## otherwise set the path
    hic_path = outpath

## Gather the parent path of the hic exsits
parent_path = '/'.join(hic_path.split('/')[:-1])

## Assert the parent path is real! 
assert os.path.exists(parent_path), "ERROR: The output parent directory (%s) does not exist!"%parent_path

## Check the extension of the hic path is .hic
assert (hic_path.split('/')[-1][-4:] == '.hic'), "ERROR: The specified extension of the output file name (%s) is not .hic!"%hic_path.split('/')[-1]

## Set the out name and path of the short file
short_path = hic_path[:-4] + '.short.txt'

# In[4]:
## Bring in other needed mods
import numpy as np, gzip

## Load in the data, nitilize lists
chrlines = []

## Set the fields
fields = np.arange(8)

## Open the gzipped file for reading
with gzip.open(inputpath,'rt') as f:
    
    ## Iterate thru the lines in file
    for line in f:

        ## If the chr1 is in the line text, then
        if coi in line:
        
            ## split the line and make into an array
            split = np.array(line.split(' '))

            ## Check the chromosome fields for chr1
            if (split[1] == coi) and (split[5] == coi) and (int(split[8]) >= quality) and (int(split[11]) >= quality):

                ## Join the line
                newline = ' '.join(split[fields])+'\n'

                ## Correct name if needed
                if newname is not None:

                    ## Make the new line
                    newline = newname.join(newline.split(coi))
            
                ## If both chr1 and chr2 are the chr of interest, then append fields of interest to list
                chrlines.append(newline)

        else: ## else pass this line
            pass
    ## Close the file
    f.close()

## Reassign the chromosome name
if newname is not None:
    coi = newname

# In[5]:
## How many contacts are on this chromosome?
ncontacts = len(chrlines)

## Print to screen
if verbose:
    print('Total number of Hi-C contacts detected: %s'%ncontacts)

    ## Check if the number of contacts is low
    if (ncontacts < tolerance):
        print('WARNING: The number of Hi-C contacts (%s) along chromosome %s is lower than set tolerance (%s).'%(ncontacts,coi,tolerance))

## Assert we have data
assert (ncontacts>0), "ERROR: No Hi-C contacts were detected along chromosome %s"%coi

## Save out the hi-c counts
if tocount:

    ## Make a path to the count txt file
    count_path = hic_path+'.counts.txt'

    ## Tell the user we are printing to file
    if verbose:
        print("WARNING: Printing counts to file: %s"%count_path)
        
    ## Open the count file path
    with open(count_path,'w') as f:

        ## Write the counts
        f.writelines('%s\t%s\n'%(coi,ncontacts))
    f.close()

## Open a file
with open(short_path,'w') as f:

    ## Write to file
    f.writelines(chrlines)

## Close the file
f.close()

## Check that the short path exists
assert os.path.exists(short_path), "ERROR: Unable to locate .short file: %s"%short_path

# In[6]:
## Load in subprocess mod
import subprocess 

## Format juicer command , set the jarpath and java call
juicer_pre = f'java -Xms512m -Xmx2048m -jar {jarpath} pre {short_path} {hic_path} {genomeid} -c {coi} -r {resolutions} -k {correction}'

## If true, print the command to screen
if verbose:
    print(juicer_pre+'\n')

## Submit juicer command to bash
subcheck = subprocess.call(juicer_pre, shell=True)

## Check the call
assert (subcheck == 0), "ERROR: Juicer pre command did not execute properly!\n%s"%juicer_pre

## Remove short file
if toclean:
    if verbose: 
        print('WARNNING: Removing .short file.')
    os.remove(short_path)

## if to clean is false and we are in verbose mode
elif not toclean and verbose:
    print("WARNING: Retaining .short file on path and with name: %s"%short_path)

else:
    pass

## Finish
print('Finished :-)')