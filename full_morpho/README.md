# full_morpho
## Description
This directory contains the source codes and exemplar data for full morphology analyses.

## Usage
### Estimation inter-neuronal and inter-regional diversity 
- Calculates the pairwise correlations between all neurons based on the their L-Measure features (`../data/features.txt`) using the scripts `./analyses/sdmatrix.py`, `./analyses/sdfeatures.py`. Additionally, the inter-regional diversity will be estimated by the average correlation among neuron pairs from the corresponding regions.


### Visualization of exemplar morphologies [`analyses/viz_swc_in_vtk_ccfbrain.py`]
Use

	python viz_swc_in_vtk_ccfbrain.py --v root.vtk --folder '../data/axon_sort'

### Whole-brain spatial-tuned clustering of full morphologies (Figure 4)
One can easily conduct the analyses following the ipynb file (`Fig4_full_morpho_final.ipynb`), and the related exemplar input files, intermediate files, and output files can be found in `../data/Fig4` folder of this project.
