#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## Converts h5 input to short format
## By: Cullen Roth
## Bring in system mod
import sys


# In[ ]:


## Set user defined variables
## Check we have three inputs!
assert (len(sys.argv) >= 4), "ERROR: This script must include:\n(1) The full path to a ginteractions (tsv) file (which is assumed to be an h5 matrix converted via HicExplorer).\n(2) A genome size (tsv) file with chromosome and size columns.\n(3) A valid output path to save the hic short file."

## Gather data inputs
datapath = str(sys.argv[1])
sizepath = str(sys.argv[2])
savepath = str(sys.argv[3])

## Set verbosity if passed
if (len(sys.argv) == 5):
    if str(sys.argv[4]) == 'true':
        verbose = True
    else:
        verbose = False
else:
    verbose = False


# ## Set user defined variables
# ## Set input path
# datapath = '/Users/croth/HIC/MRC5/2401.006.h5.toremove.ginteractions.tsv'
# 
# ## Set output path
# savepath = '/Users/croth/HIC/MRC5/2401.006.h5.toremove.short'
#     
# ## Set path to size file
# sizepath = '/Users/croth/REFERENCES/ENCODE/genome.size.txt'
# #sizepath = '/Users/croth/REFERENCES/ENCODE/test1.size.txt'
# #sizepath = '/Users/croth/REFERENCES/ENCODE/test2.size.txt'
# 
# ## Set verbose 
# verbose = False

# In[ ]:


## Set other needed variables
## Set verbosity
#verbose = True

## Set input sep
mysep = '\t'

## Set output output sep
outsep = ' '

## Set column names
colname =  ['Chrom1','Left1','Right1','Chrom2','Left2','Right2','Quality']


# In[ ]:


## Bring in needed mods
import pandas as pd, numpy as np

## Write a ftn to check index between two dataframes
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
genomesize = pd.read_csv(sizepath,sep=mysep,names=['Chrom','Size'])

## Make a list of chromosomes
chrlist = genomesize.Chrom.tolist()


# In[ ]:


## Load in and set columns
temp = pd.read_csv(datapath,sep=mysep,header=None,names=colname)

## Take total contact counts
contacts = temp.shape[0]

## Print size of temp file
if verbose:
    print('Detected %s HiC contacts.'%contacts)
    
if (contacts == 0):
    print('ERROR: No HiC contacts detected!')
    sys.exit(1)


# In[ ]:


## Subset data for data in genomesizes file
temp = temp[(temp[colname[0]].isin(chrlist)) & (temp[colname[3]].isin(chrlist))].reset_index(drop=True)

## Gather the new index after dropping samples
theindex = temp.index.values

## Number of contacts dropped
ndrop = contacts - temp.shape[0]

## calculate total number of conatacts dropped
nperc = np.round(100*ndrop/contacts,3)

## Print the number of dropped contacts
if verbose:
    print("WARNING: Removed %s ( %s"%(ndrop,nperc) + " % ) contacts from unlisted chromosomes."  )


# In[ ]:


## Check that we have contacts for all chromosomes in chrlist
## Gather chromosomes still in the filtered h5
tempchrlist = list(np.unique(np.concatenate([temp[colname[0]].unique(),temp[colname[3]].unique()])))

## Gather the names of the missing chromosomes
missing = [c for c in chrlist if c not in tempchrlist]

## If any chromosomes are missing
if len(missing) > 0:
    print("WARNING: No contacts were detected for chromosomes:")
    print("\n".join(missing))


# In[ ]:


## Split by contact type
## Log if verbose
if verbose:
    print("Splitting inter- & intra-chromosomal contacts.")
    
## Gather the between chrom contacts
inter = temp[(temp.Chrom1!=temp.Chrom2)]

## Check the shape and number of inter-chromosome contacts
if verbose and (inter.shape[0] == 0):
    print("WARNING: Zero inter-chromosomal contacts detected.")
else:
    print("Number of between chromosome contacts: %s"%inter.shape[0])
    
## Gather the within chromosome contacts
intra = temp[(temp.Chrom1==temp.Chrom2)]

## Check the shape and number of intra-chromosome contacts
if verbose and (intra.shape[0] == 0):
    print("ERROR: Zero intra-chromosomal contacts detected.")
    sys.exit(1)
else:
    print("Number of within chromosome contacts: %s"%intra.shape[0])
    
## What is the ratio of intra vs inter
if verbose and (intra.shape[0] > 0):
    
    ## Calculate ratio
    interintra = np.round(100*inter.shape[0]/intra.shape[0],3)
    
    ## Print to screen
    print('Ratio of inter- to intra-chromosome contacts: %s %s'%(interintra,'%'))


# In[ ]:


## Correct intra chromosomal contacts
## Remove temp
del temp

## Log if verbose
if verbose:
    print("Sorting intra-chromosomal contacts.")

## Sort the within chromcontacts by chromosome and left read postition
intrac = pd.concat([intra[(intra.Chrom1==c)].sort_values('Left1') for c in chrlist])

## Delete the old intra
del intra


# In[ ]:


## Split inter chromosome contacts into left and right pairs
## Log status
if verbose and (inter.shape[0]>0):
    print("Gathering pairs of inter-chromosomal contacts.")

## Gather left 
left = inter[inter.columns[:3]]

## Check work
assert (left.shape[1] == 3), "ERROR: Missing columns of left pairs.\nThere should be three and there are %s"%left.shape[1]

## Gather right
righ = inter[inter.columns[3:-1]]

## Check work
assert (righ.shape[1] == 3), "ERROR: Missing columns of right pairs.\nThere should be three and there are %s"%righ.shape[1]

## Take the correction index
tocorrect = inter.index.values

## Take the quality of between chromosome contacts
interquality = inter[colname[-1]]


# In[ ]:


## Reorder pairs of inter chromosomal contacts
if verbose and (inter.shape[0]>0):
    print("Reordering inter-chromosomal contacts by chromosome.")

## Initilize inter list
inter = []

## Iteratively correct the inter chromosome names
for i in tocorrect:
    
    ## Gather chromosome names from 
    ## The left pair and .. 
    c1 = left.loc[i,colname[0]]
    
    ## the right pair of the inter chromosome contact
    c2 = righ.loc[i,colname[3]]
    
    ## Gather chromosome index of the left read and ..
    c1ix = genomesize[(genomesize.Chrom==c1)].index.min()
    
    ## the right read of the pair in contact
    c2ix = genomesize[(genomesize.Chrom==c2)].index.min()
    
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
inter = pd.DataFrame(inter,columns=colname[:-1],index=tocorrect)

## Check that we have the same size dataframe
assert (inter.shape[0] == left.shape[0])


# In[ ]:


## Sort inter pairs by chromosome positon
if verbose and (inter.shape[0]>0):
    print("Sorting inter-chromosomal contacts by chromosome.")

## Initilize corrected inter (between) chrom contact list
interc = []

## Gather list of chromosomes with trans contacts
interchrs = [c for c in chrlist if c in inter[colname[0]].tolist()]

for c in interchrs:
    
    ## Slice the single chromosome
    temp = inter[(inter.Chrom1==c)]
    
    ## Gather the inter chromosomes
    interchrom = genomesize[(genomesize.Chrom.isin(temp[colname[3]].unique()))].Chrom.tolist()
    
    ## Sort the right side of the interchromosomes
    tempc = pd.concat([temp[(temp[colname[3]]==ic)].sort_values([colname[1],colname[4]]) for ic in interchrom])

    ## append to the corrected between chromosome contact list
    interc.append(tempc)
    
## concatonate into a dataframe
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
    assert (interc.shape[0] == 0)


# In[ ]:


## Combine both sorted inter and intra by sorted chromosome in chrlist
if verbose:
    print("Blocking contacts of %s chromosome(s)."%len(chrlist))

## Initilize list
hic = []

## Set counter
ci = 0

## Iterate thru each chromosome
for c in chrlist:
    
    ## Slice intra (within)
    temp1 = intrac[(intrac[colname[0]]==c)]
    
    ## Slice inter (between)
    temp2 = interc[(interc[colname[0]]==c)]
    
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
    hic.append(temp)
    
    ## Count
    ci += 1
    
## Check our count
assert ci == len(chrlist)
    
## make into a dataframe
hic = pd.concat(hic)

## Check the final shape
assert (hic.shape[0] == len(theindex)), "ERROR: There are missing valid HIC contacts!"

## Check inter chrom contacts last column
checkix(hic[(hic[colname[-1]].isna())],interquality)

## Reassign last column to inter chrom contacts
hic.loc[interquality.index,colname[-1]] = interquality.values

## check our assignment
assert (hic.dropna().shape[0] == hic.shape[0]), "ERROR: There is missing data in the HIC dataframe!"

## Check final index
checkix(hic,pd.DataFrame(index=theindex))


# In[ ]:


## Generate a short file
if verbose:
    print("Generating hic short file: %s"%savepath)

## gather colunm names to be held over
convertix = np.array([0,1,3,4,6])

## Make new column names
newcols = ['buffer1'] + hic.columns[:2].tolist() + ['buffer2','buffer3'] + hic.columns[3:5].tolist() + ['buffer4'] + hic.columns[-1:].tolist()

## Check that their are nine of these
assert len(newcols) == 9, "ERROR: The short file columns were not generated correctly."

## Initilize short dataframe
short = pd.DataFrame(columns=newcols,index=hic.index)

## For each old column name
for c in colname:
    
    ## If its in the new short dataframe assigne it
    if c in newcols:
        short[c] = hic[c]
        
    else:
        pass
    
## Assign zeros to buffer columns 1,2, and 3
short[['buffer1','buffer2','buffer3']] = 0

## and a one to buffer column 4
short[['buffer4']] = 1

## Convert all the columns except those with the chromosome name to integers
## Gather columns to be converted
toint = [c for c in short.columns if c not in [colname[0],colname[3]]]

## Convert to integers
for c in toint:
    short[c] = short[c].apply(int)    
    
## Check that we didn't lose any records
checkix(short,hic)

## SAve out dataframe
short.to_csv(savepath,sep=outsep,header=False,index=False)

## Print finish
if verbose:
    print("Finished :D")

