#!/usr/bin/env python
# coding: utf-8
## gin.to.short.py
## 
## Converts ginteractions input to short format to help with hicexplorer to juicer compatability.
##
## Developed by:
##
## Cullen Roth, Ph.D.
##
## Postdoctoral Research Associate
## Genomics and Bioanalytics (B-GEN)
## Los Alamos National Laboratory
## Los Alamos, NM 87545
##
## For help with any issues please email: croth@lanl.gov

# In[0]:

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
required.add_argument("-i", type=str, required=True, help="Path to input .summary.txt.gz file for analysis.\nThis script assumes the file is g-zipped.\nMust provide a relative path.", metavar='./path/to/ginteractions.tsv')
required.add_argument("-g", type=str, required=True, help="Path of the chrom.sizes.tsv file that lists the name and size of the chromosomes in order.")
required.add_argument("-o", type=str, required=True, help="Output path and name of the short file generated from this script.")

## Set optional, boolean variables
optional.add_argument("-V", help="Flag to run in verbose mode.\nDefault behavior is false.", action='store_true')

## Set the paresed values
values = parser.parse_args()

## Set input values
ginpath, sizepath, savepath, verbose = values.i, values.g, values.o, values.V

# In[1]:

## Bring in needed mods
import pandas as pd, numpy as np

## Set other needed variables and functions
## Set column names of ginteractions file and other vars
ginnames, outsep, tsep = ['chrom1','left1','right1','chrom2','left2','right2','hiccount'], ' ', '\t'

## write ftn for loading in ginteractions file
def loadgin(path):
    return pd.read_csv(path, sep=tsep, names = ginnames)

## Write ftn for expanding the ginteractions file
def expandgin(df):

    ## Initilize list
    newrows=[]
    
    ## Iterate thru the unique counts
    for k,temp in df.groupby('hiccount'):

        ## Initiliase iterator and iterate
        n = 0
        for i in np.arange(1,k+1):

            ## Append the new row
            newrows.append(temp)

            ## Add to our iterator
            n = n + 1
        
        ## Check our work
        assert (n == k), "ERROR: There was an error during expansion!"
    
    ## Concat rows
    newrows = pd.concat(newrows,ignore_index=True).reset_index(drop=True)
    
    ## Check our work
    assert (newrows.shape[0] == df.hiccount.sum()), "ERROR: Missing hi-c counts!"
    
    ## Return the new rows
    return newrows

## Write a ftn for converting ginteractions to short
def gintoshort(df):
    
    ## Expand the dataframe
    exdf = expandgin(df)
    
    ## Add a strand and fragment columns
    exdf['str1'], exdf['str2'], exdf['frag1'], exdf['frag2'] = 0, 1, 0, 1
    
    ## Return expanded df
    return exdf[['str1','chrom1','left1','frag1','str2','chrom2','left2','frag2']]

## Write ftn for checking index
def checkix(x,y):
    x = np.array(sorted(x.index.values))
    y = np.array(sorted(y.index.values))
    assert (np.sum(x-y) == 0), "ERROR: The indices of the dataframes to not match!"

# In[ ]:

## Load in genomesize and contact data
## Log if verbose
if verbose:
    print("Loading genome size and contact (h5) files.")

## Load genome size file
genomesize = pd.read_csv(sizepath,sep=tsep,names=['Chrom','Size'])

## Make a list of chromosomes
chrlist = genomesize.Chrom.tolist()

# In[ ]:

## Load in and set columns
temp = loadgin(ginpath)

## Take the total contact counts
contacts = temp.hiccount.sum()

## Print size of temp file
if verbose:
    print('Detected %s hi-c contacts.'%contacts)

## Make sure we have contacts! 
assert (contacts > 0), "ERROR: No hi-c contacts detected!"   

# In[ ]:

## Subset data for data in genomesizes file
temp = temp[(temp['chrom1'].isin(chrlist)) & (temp['chrom2'].isin(chrlist))].reset_index(drop=True)

## Count these hic records
chrlist_count = temp.hiccount.sum()

## Number of contacts dropped
ndrop = contacts - chrlist_count

## calculate total number of conatacts dropped
nperc = np.round(100*ndrop/contacts,3)

## Print the number of dropped contacts
if verbose:
    print("WARNING: Removed %s ( %s"%(ndrop,nperc) + " % ) contacts from unlisted chromosomes."  )

# In[ ]:

## Check that we have contacts for all chromosomes in chrlist
## Gather chromosomes still in the filtered h5
tempchrlist = list(np.unique(np.concatenate([temp['chrom1'].unique(),temp['chrom2'].unique()])))

## Gather the names of the missing chromosomes
missing = [c for c in chrlist if c not in tempchrlist]

## If any chromosomes are missing
if len(missing) > 0:
    print("WARNING: No contacts were detected for chromosomes:")
    print("\n".join(missing))

# In[ ]:

## Split by contact type
## Expand temp and make into a ginteractions file
temp = gintoshort(temp)

## Count the expanding, each row represents a hi-c contact
short_counts = temp.shape[0]

## Check work
assert (short_counts == chrlist_count), "ERROR: The number of hi-c counts was corrupted after expansion!"

## Log if verbose
if verbose:
    print("Splitting inter- & intra-chromosomal contacts.")
    
## Gather the between chrom contacts
inter = temp[(temp.chrom1!=temp.chrom2)]

## Check the shape and number of inter-chromosome contacts
if verbose and (inter.shape[0] == 0):
    print("WARNING: Zero inter-chromosomal contacts detected.")
else:
    print("Number of between chromosome contacts: %s"%inter.shape[0])
    
## Gather the within chromosome contacts
intra = temp[(temp.chrom1==temp.chrom2)]

## Check that there are contacts within chromosomes
assert (intra.shape[0] != 0), "ERROR: Zero intra-chromosomal contacts detected."

## Check the shape and number of intra-chromosome contacts
if verbose:
    print("Number of within chromosome contacts: %s"%intra.shape[0])
    
## What is the ratio of intra vs inter
if verbose and (intra.shape[0] > 0):
    
    ## Calculate ratio
    interintra = np.round(100*inter.shape[0]/intra.shape[0],3)
    
    ## Print to screen
    print('Ratio of inter- to intra-chromosome contacts: %s %s'%(interintra,'%'))

## Correct the order of the intra chromosomal contacts, remove temp it is no longer needed
del temp

## Log if verbose
if verbose:
    print("Sorting intra-chromosomal contacts.")

## Sort the within chromcontacts by chromosome and left read postition
intrac = pd.concat([intra[(intra.chrom1==c)].sort_values('left1') for c in chrlist])

## Delete the old intra
del intra

# In[ ]:

## Split inter chromosomal (between chromosomes) contacts into left and right pairs
## Log status
if verbose and (inter.shape[0]>0):
    print("Gathering pairs of inter-chromosomal contacts.")

## Gather left 
left = inter[['str1','chrom1','left1','frag1']]

## Gather right
righ = inter[['str2','chrom2','left2','frag2']]

## Set the new gin names
newginnames = ['str1','chrom1','left1','frag1','str2','chrom2','left2','frag2']

## Take the correction index
tocorrect_ix = inter.index.values

# In[ ]:

## Reorder pairs of inter chromosomal contacts
if verbose and (inter.shape[0]>0):
    print("Reordering inter-chromosomal contacts by chromosome.")

## Initilize inter list
inter = []

## Iteratively correct the inter chromosome names
for i in tocorrect_ix:
    
    ## Gather chromosome names from the left and right pair of the inter chromosome contacts
    c1, c2 = left.loc[i,'chrom1'], righ.loc[i,'chrom2']
    
    ## Gather chromosome index of the left and right read of those in contact
    c1ix, c2ix = genomesize[(genomesize.Chrom==c1)].index.min(), genomesize[(genomesize.Chrom==c2)].index.min()

    ## If the "Left" chromosome is the first in order make in this order
    if (c1ix < c2ix):
        newline = left.loc[i].tolist() + righ.loc[i].tolist()
        
    ## Else if "right" chromosome is the first in order make in this order
    else:
        newline = righ.loc[i].tolist() + left.loc[i].tolist()
        
        ## assert that the chromosomes may not have the same index
        assert (c1ix != c2ix), "ERROR: The chromosomes are not inter-chromosomal contacts! "
        
    ## append to inter list 
    inter.append(newline)
    
## Make list into dataframe
inter = pd.DataFrame(inter,columns=newginnames,index=tocorrect_ix)

## Check that we have the same size dataframe
assert (inter.shape[0] == left.shape[0]), "ERROR: We lost contacts between chromosomes"

# In[ ]:

## Sort inter pairs by chromosome positon
if verbose and (inter.shape[0]>0):
    print("Sorting inter-chromosomal contacts by chromosome.")

## Gather list of chromosomes with trans contacts and initlizse corrected inter chrom contact list
interchrs, interc = [c for c in chrlist if c in inter['chrom1'].tolist()], []

## Iterate thru the interacting chromosomes
for c in interchrs:
    
    ## Slice the single chromosome
    temp = inter[(inter.chrom1==c)]
    
    ## Gather the inter chromosomes
    interchrom = genomesize[(genomesize.Chrom.isin(temp['chrom2'].unique()))].Chrom.tolist()
    
    ## Sort the right side of the interchromosomes
    tempc = pd.concat([temp[(temp['chrom2']==ic)].sort_values(['left1','left2']) for ic in interchrom])

    ## append to the corrected between chromosome contact list
    interc.append(tempc)
    
## Concatonate into a dataframe
if (len(interc)>0):
    interc = pd.concat(interc)
    
    ## Check our work
    assert (inter.shape[0] == interc.shape[0])
    
    ## Check the index
    checkix(inter,interc)
    
    ## Delete memory hogs
    del tempc

else:
    ## Set interc to the empty dataframe made above
    interc = inter
    
    ## Check work
    assert (interc.shape[0] == 0), "ERROR: Inter-chromosomal contacts are detected where none were expected!"

# In[ ]:

## Combine both sorted inter and intra by sorted chromosome in chrlist
if verbose:
    print("Blocking contacts of %s chromosome(s)."%len(chrlist))

## Initilize list and set counter
short, ci = [], 0

## Iterate thru each chromosome
for c in chrlist:
    
    ## Slice intra (within) and inter (between) contacts
    temp1, temp2 = intrac[(intrac['chrom1']==c)], interc[(interc['chrom2']==c)]
    
    ## Print a warning if both intra and inter chrom contacts are zero!
    if (temp1.shape[0]==0) and (temp2.shape[0]==0):
        print('WARNING: No contacts found for %s'%c)
        continue
    
    ## If there are no between chrom contacts
    if (temp2.shape[0]==0):
        
        ## Set new temp to just the within chrom contacts
        temp = temp1
        
    ## Other wise concatinate them
    else:
        temp = pd.concat([temp1,temp2])
        
    ## append to list
    short.append(temp)
 
## make into a dataframe
short = pd.concat(short)

## Check the final shape
assert (short.shape[0] == short_counts), "ERROR: There are missing valid hi-c contacts!"

## check our assignment
assert (short.dropna().shape[0] == short.shape[0]), "ERROR: There is missing data in the hi-c dataframe!"

# In[ ]:

## Generate a short file
if verbose:
    print("Generating hi-c short file: %s"%savepath)

## Gather columns to be converted to integers
toint = [c for c in short.columns if c not in ['chrom1','chrom2']]

## Convert these columns to integers
for c in toint:
    short[c] = short[c].apply(int)    
    
## SAve out dataframe
short.to_csv(savepath,sep=outsep,header=False,index=False)

## Print finish
print("Finished :D")