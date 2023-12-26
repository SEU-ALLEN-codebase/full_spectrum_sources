# Bouton
This section includes scripts of two main components: the first involves the bouton morphology analysis codes, specifically those generating Figure 7 of the paper (`Fig7_bouton_analysis.ipynb`). The second set of scripts is dedicated to generating bouton features for cross-scale analyses. Specifically:

## Source codes
1. **Cross-Scale Bouton Feature Generation** 
    - `analyses/bouton_analyzer.py`: A Python script designed for extracting bouton (synaptic terminal) features tailored for cross-scale analyses. Exemplar input file and output files can be found at `./data/bouton_v20231211_swc` and `./data/bouton_features.csv`
    - `analyses/sdfeature.py` and `analyses/sdmatrix.py`: Scripts responsible for aggregating regional-level and neuronal-level bouton features, as well as creating similarity matrices between different neurons and neuron types. 

2. **Bouton Morphology Analysis** 
    The Jupyter Notebook `Fig7_bouton_analysis.ipynb` contains codes for analyzing and visualizatin of bouton preference on various morphological scales (Figure 7 of the paper). Examplar and intermediate files are stored at `data/Fig7`. Users can execute the entire analysis step-by-step by utilizing the notebook script.

## Morphological data
The morphologies containing boutons/varicosities (in format of *SWC* file), are archived on both GitHub and Google Drive. The archive includes morphologies in their original image space (`RAW.zip`) and in isotropic 1um CCFv3 space (`CCFv3.zip`). Boutons are identified as nodes with a *type* value of 5 in the second column. The morphologies on GitHub are released in conjunction with this project and can be found in the `./R1876_dataset` directory. These datasets are sparsely sampled using 20 voxels or 10um for the original space version and the registered version, respectively.  For access to the un-sampled version on Google Drive, please visit this link: https://drive.google.com/drive/folders/1jrLMmEUjvl-_XR6c5rDnhcm5DNhLJpSq.

