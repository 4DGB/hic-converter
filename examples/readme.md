# Example Plotting in Python

The jupyter notebook [Chromosome_3D_plotting](Chromosome_3D_plotting.ipynb) contains exmaples of calling python functions from our custom python library, [chromoviz](chromoviz.py). 
These functions rely on [pyvista](https://pyvista.org) to generate static images (.png) of 3D chormosome structures. 
The chromoviz library contains several helper functions for loading in 3D structure files, fitting 3D splines to points, and aligning chromosome models in 3D space. 

## Create jupyter environment
With the environment.yml file in this folder make the python environment for plotting 3D structures from the 4DGBWorkflow. 
```
conda env create -f environment.yml -n chromo-env
```

## Activate jupyter lab 
After generating this environment, call jupyter lab
```
jupyter lab
```

## Data used
We host here the [relative locations](CHM13v2.0.centromeres.bed) of centromeres in the T2T human genome reference. 
These are used to replot and correct the 3D structures around centromeres, which tend to be unmappable, forming long unbound loops/chains in the 3D models. 
The actual Hi-C data was taken from [Venu et al.](https://link.springer.com/article/10.1186/s13072-025-00631-4), and mapped to the T2T human refernece genome. 
These data, for chromosome X were processed (in .hic format) with [4DGBWorkflow](https://github.com/4DGB/4DGBWorkflow) and the 3D structures (in .csv format) from LAMMPS are stored [here](LAMMPS). 