# arbors
## Description

The automatic arbor detection algorithm operates under the assumption that a morphology can be represented as a graph, wherein each node corresponds to a vertex, and the parent-child compartment relations are treated as edges. It leverages spectral clustering techniques to partition all nodes effectively. 

Within our current framework, we treat each apical dendrite and basal dendrite as distinct arbors, and the axons are segmented into one or several arbors using the autoarbor script (`autoarbor_v1_yf.py`). While users can utilize original version of the registered version of SWCs (SEU-A1891) as input for the arborization script, we highly recommend employing down-sampled SWCs to enhance calculation speed. In our own practice, we down-sampled each SWC, achieving a morphology with approximately 80-micron spaced intervals, utilizing the `resample_swc` plugin and sorting the nodes using `sort_neuron_swc` plugin available in Vaa3D. 

For testing purposes, we have thoughtfully included several down-sampled examples within the 'data_examples' folder. These examples are readily accessible, and anyone is welcome to explore and experiment with them.
## Usage
### Axonal arbors generation
The generation of axonal arbors can be divided into 3 steps:
- Run the initial round of automatic arbor detection using the `run_arbor_r2.py` script with the specified parameters, you can run the following command in bash: 
        `python run_arbor_r2.py --r 1 --L 1 --H 4`
    This command will execute the script with the following parameter settings:
        --r 1: This sets the value of r to 1.
        --L 1: This sets the lower bound of the range for selecting neurons to 1.
        --H 4: This sets the upper bound of the range for selecting neurons to 4.
    This command will generate the arbors in the format of SWC files, and their statistics in a separate `txt` file.

- Calculate the median number of arbors for each projection-differentiated soma-types. Based on the statistical `txt` files generated in the previous step and the neuronal categorizing file: "../../commlib/41586_2021_3941_MOESM4_ESM.csv", you can calculate the median number of arbors for each neuronal subtypes. An example of such file can be found in "log/axonal_arbor_params_round2.csv"
- Run the second round of automatic arbor detection using the median arbor numbers with the following command in your Bash terminal:
        `python run_arbor_r2.py --r 2 --L 0 --H 0`
    You can create a new folder, for example, "../data/axon_arbors_round2", and move all generated files into the folder. 
### Arbor feature extraction
To extract features using the arbor_analyzer_r2.py script located in the "analyses" directory, you can follow these instructions:

To Get Features for All Soma-Types:

Execute the following command to extract features for all soma-types:

    `python arbor_analyzer_r2.py`
This command will run the feature extraction script and generate feature files for all soma-types.

To Get Features for Projection-Differentiated Soma-Types:

Set the `soma_type_merge` argument to False in the script. You may need to edit the script to change this parameter.

Re-run the script with the updated parameter:

    `python arbor_analyzer_r2.py`
This command will run the feature extraction script with the soma_type_merge set to False, resulting in feature files specifically for projection-differentiated soma-types.

After running these commands, you should have two sets of feature files for each neurite type (axonal, basal, or apical):

features_r2_somaTypes_{neurite_type}.pkl: Features for all soma-types.
features_r2_projAndSomaTypes_{neurite_type}.pkl: Features for projection-differentiated soma-types.
Make sure you are in the correct working directory where the script arbor_analyzer_r2.py is located and where the data files are present.
### Arbor plotting
1. Generate the feature dotted heatmap for arbor features: Turn on the specific plotting functions in the script `arbor_analyzer_r2.py`. 
2. Other panels can be found in the R scripts and ipython notebook files in the `analyses` directory