## Bouton
This section includes scripts of two main components: the first involves the bouton morphology analysis codes, specifically those generating Figure 7 of the paper (`Fig7_bouton_analysis.ipynb`). The second set of scripts is dedicated to generating bouton features for cross-scale analyses. Specifically:

1. Bouton Morphology Analysis: `Fig7_bouton_analysis.ipynb`: A Jupyter Notebook containing codes for analyzing bouton morphology, resulting in Figure 7 of the paper. 
    - Some of examplar and intermediate files are stored at `data/Fig7`

2. Cross-Scale Bouton Feature Generation: 
    - `analyses/bouton_analyzer.py`: A Python script designed for extracting bouton (synaptic terminal) features tailored for cross-scale analyses. 
    - `analyses/sdfeature.py` and `analyses/sdmatrix.py`: Files responsible for aggregating regional-level and neuronal-level bouton features, as well as creating similarity matrices between different neurons and neuron types. 
    - `analyses/convert_eswc_to_swc.py`: A script attempting to convert extended swc files (*ESWC*) that contain additional information, such as bouton types, to standardized *SWC* files.
