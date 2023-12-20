# full_morpho
## Description
This directory contains the source codes and exemplar data for full morphology analyses.

## Usage
### Estimation inter-neuronal and inter-regional diversity [`./analyses/sdmatrix.py`, `./analyses/sdfeatures.py`]
This Python script calculates the pairwise correlations between all neurons based on the their L-Measure features (`../data/features.txt`). Additionally, the inter-regional diversity will be estimated by the average correlation among neuron pairs from the corresponding regions.


### Visualization of exemplar morphologies [`analyses/viz_swc_in_vtk_ccfbrain.py`]
Use

	python viz_swc_in_vtk_ccfbrain.py --v root.vtk --folder '../data/axon_sort'

### Generation of Figure 4 [`Fig4_full_morpho_final.ipynb`]
One can easily conduct the analyses following the ipynb file, and the related exemplar input files, intermediate files, and output files can be found in `data_notebook` folder of this project.
