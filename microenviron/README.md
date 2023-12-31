# microenviron
## Usage
### Ⅰ. Local Morphology Generation (Optional)
The first step in constructing the microenvironment involves the generation of morphologies, which can be dendrites, complete morphologies, or any sub-neuronal arbors derived from full morphologies. In this project, we have utilized automatically traced local morphologies to represent local image blocks. As such, the initial phase necessitates the acquisition of an adequate number of morphologies from neuronal images. Below, we provide a comprehensive set of exemplary procedures for tracing neurons within extensive neuronal blocks.
#### Automatic Tracing Codes Location
The automatic tracing codes can be found at `generation/42k_script`. Prior to utilizing these codes, please ensure that the following prerequisites are met:

- Cropped images with dimensions of 512x512x256 voxels (in `xyz` order, corresponds to ~236x236x512um^3).

#### Procedure
1. **Removal of Dense Images with Multiple Somas**
   
   - Location: `generate/42k_script/rm_multi`
   - This is an optional step and is only necessary if there are numerous image blocks with multiple somas.
   - To perform this step, you should provide at least one text file containing the coordinates of all somas. Examples of such files include "/home/vkzohj/data/enhanced_soma_images/all.txt" or the pickled file "coords.pickle".
   - You can estimate the grid-based soma counts matrix using the script "statis.py".
   - Finally, using these files, you can filter out all dense image blocks and save the remaining image blocks into a text file "passed.txt" using a similar approach as demonstrated in the script "filter.py".

2. **Enhance the Image Using "imPreProcess" Plugin**
   
   - Location: `generation/42k_script/enh/run.sh`
   - The enhancement method "imPreProcess" is bound with Vaa3D.
   - This enhancement process is detailed in the paper "Image enhancement to leverage the 3D morphological reconstruction of single-cell neurons" (Guo et al., Bioinformatics, 2022, [doi: 10.1093/bioinformatics/btab638](https://doi.org/10.1093/bioinformatics/btab638)).

3. **Tracing with APP2**
   
   - Location: `generation/42k_script/app2`
   - APP2 is a popular automatic tracing algorithm (Xiao et al., bioinformatics, 2013, [doi: 10.1093/bioinformatics/btt170](https://doi.org/10.1093/bioinformatics/btt170)), and it is available in Vaa3D.
   - Exemplary batch running codes can be found in the specified directory.

4. **Tracing with neuTube**
   
   - Location: `generation/42k_script/ntb`
   - Similar to step 3 (APP2 tracing), but using another popular auto-tracing algorithm neuTube (Feng et al., eNeuro, [doi: 10.1523/ENEURO.0049-14.2014](https://doi.org/10.1523/ENEURO.0049-14.2014)).
   - Batch running codes are located in the specified directory.

5. **Pruning the APP2 Reconstructions**
   
   - Location: `generation/42k_script/seg_prune`
   - Refer to the script "run.sh" for guidance on pruning APP2 reconstructions.

6. **Consensus Calculation**
   
   - Location: `generation/42k_script/consensus`
   - Estimation of the consensus morphology between pruned APP2 reconstruction and neuTube reconstruction.
   - Relevant codes can be found in the specified directory.

### Ⅱ. Microenvironment Construction
A microenvironment is a spatially-tuned ensemble of neighboring neurons. In this project, we illustrate the process of generating microenvironments using 5 nearby neurons, each represented by their morphologies within 512x512x256 image blocks. Here are the steps for generating microenvironments:

1. **Calculate L-Measure Features for All Morphologies**

   - We utilize the `global_neuron_feature` plugin integrated with Vaa3D to calculate L-Measure features for each neuron. You may opt to use other morphological features instead of L-Measure features. Alternatively, you can estimate L-Measure features using the L-Measure server (http://cng.gmu.edu:8080/Lm), which implements some features differently. This step results in obtaining 22-dimensional features for each neuron.
   - Script: `calc_global_features.py`
   - Location: `generation`
   - Input: swc_dir - Directory containing all registered morphologies. One exemplar reconstructed local morphology (`Img_X_6116.92_Y_10268_Z_4371.18.v3dpbd_stps.swc`) can be found under the directory `data`, and the whole set of reconstructions are archived under the `42k_dataset` folder
   - Output: L-Measure feature file for all neurons in each brain. An exemplar output file is: `lm_features/lm_d22_150k_0602_201605.csv`

2. **Aggregate Additional Semantic Information**

   - This step involves gathering additional semantic information, such as brain regions for all neurons and their soma locations. This information will be merged with all features to create a unified feature file containing features for all neurons.
   - Script: `preprocessing.py`
   - Location: `generation`
   - Input:
     - swc_dir - Directory containing all morphologies
     - feature_dir - Directory containing all feature files, with each brain having a feature file (containing multiple neurons)
   - Output: Merged feature file (e.g., `lm_features_d22_all.csv`). An exemplar file containing 10 neurons can be found at `data/lm_features_d22_all_part.csv`

3. **Construction of Microenvironments**

   - As outlined in our manuscript, we first identify the top 5 neurons within a given radius. This radius is estimated as the 50th percentile of all distances between neurons and their top 5 nearest neighbors (refer to the "estimate_radius" function). A microenvironment is then defined as the collection of the target neuron and the selected top 5 (at most) neurons. The feature of a microenvironment is a distance-weighted summary of all neurons within that microenvironment. Here we used an heuristic 
   - Script: `micro_env_features.py`
   - Location: `generation`
   - Input: The merged feature file (e.g., `./data/lm_features_d22_all.csv`)
   - Output: Microenvironment feature file (e.g., `./data/micro_env_features_nodes300-1500_withoutNorm.csv`)

### Ⅲ. Microenvironment Analysis
1. **Select the Most Discriminating Features Using mRMR**

   - Based on the microenvironment features and their soma locations (soma-types), we identify the most discriminating features to facilitate visualization. We recommend using the top 3 features for ease of visualization in colorized 2D images.
   - Script: `mRMR_fsel.py`
   - Location: `analyses`
   - Input: Microenvironment feature file (e.g., `micro_env_features_nodes300-1500_withoutNorm.csv`)
   - Output: Features sorted by the mRMR score

2. **Plot Microenvironment Feature Landscapes Across the Whole Mouse Brain**

   - Script: `generate_me_map.py`
   - Location: `analysis`
   - Generation of feature maps: `generate_me_maps`
   - Calculate feature transition along different paths:
     - For paths 1 to 3: Run the function `sectional_dsmatrix`, and then `plot_me_dsmatrix`
     - For path 4 (within CP region): Use `feature_evolution_CP_radial`
   - Input: Microenvironment feature file (e.g., `micro_env_features_nodes300-1500_withoutNorm.csv`)

3. **Evaluation of reconstructed neurons**

    - Firstly, crop the manual annotations (gold standards), using the script: `analyses/gold_standards/crop.py`
    - Secondly, find out the correspondence of neurons from the full morphologies and local morphologies based on the brain and soma positions. This can be done using the script `analyses/gold_standards/utils/match_gs_recon.py`. Parameters are:
        - "gs_dir": The directory containing the morphologies (in format of SWC) with coordinates in original image space
        - "crop_dir": directory containing cropped SWC in 25um isotropic coordinate space
        - "recon_dir": directory containing reconstructed SWCs
        - "outfile": output file
    - Then, calculate the morphological features of the matched cropped gold standards. This can be done by running the script: `analyses/gold_standards/calc_glob_features.py`
    - Bench-testing the reconstructed neurons using distance metrics (e.g. spatial distance, substantial spatial distance, and percentage of different structures), or topological morphology descriptor (TMD) by running the script: `analyses/gold_standards/benchmark.py`
    - Finally, plot the morphological features distributions and TMD distributions, by running the script: `analyses/gold_standards/plot.py`

4. **Estimate the microenvironment feature for each manual annotation**
    The microenvironment features of a manually annotated single neuron is estimated by finding the microenvironment with the most similar L-measure features from the SEU-D15K. This is done by running the script `analyses/predict_types.py`



