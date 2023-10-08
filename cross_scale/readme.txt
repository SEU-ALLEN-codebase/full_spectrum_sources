# Feature generation of morphometry at various scales


# Analyses
1. Relationship between neuron distances (soma-to-soma, axon-to-axon) and neuron similarities
Two types of neuron distances: soma-to-soma and axon-to-axon distances are anlyzed, where the axon-to-axon "distance" is estimated by the cosine distance between their axonal projections. The neuron similarity of two neurons is the pearson correlation between their multiscale features.
    - script: "analyses/distance_vs_correlation.py"
    - input: 
        - soma positions for all features, "../../common_lib/misc/soma_pos_1891_v20230110.csv"
        - correlations between neurons, "../multi-scale/corr_neuronLevel_sdmatrix_heatmap_stype_all.csv"
        - projection file: "../../common_lib/41586_2021_3941_MOESM4_ESM_proj.csv"
        - meta informations for the neurons: "../../common_lib/41586_2021_3941_MOESM4_ESM.csv"
    
2. 
    


