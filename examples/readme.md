# Example Plotting in Python

The jupyter notebook [Chromosome_3D_plotting].ipynb contains exmaples of calling python functions from the python library chromoviz.py. 
These functions rely on pyvista to generate static images of 3D chormosome structures. 

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