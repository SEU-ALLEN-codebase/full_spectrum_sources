# motif
## Description
This directory contains the source codes for extracting the primary axonal tract of each neuron, as well as evaluation of their diversity and stereotypy across the whole-brain. 

## Usage
### Generation [`./generation`]
1. **Extraction of the primary axonal tracts from digitalized neuron morphologies in the format of SWC.** To do so, run the script `extract_main_path.py`, after correctly configuring the parameters. Parameters are:
    - swc_dir: directory containing all SWC files. The SWC files must be registered to the standardized CCFv3 atlas space, so that they be comparably analyzed in the subsequent steps. An example of neuron file can be found in the `../data` directory.
    - out_dir: output directory for saving the extracted axonal projection path files, in the format of SWC. 
    - ctype_file: a CSV file containing the manually curated neuron types, projection subtypes, and cortical laminations. 
    - scale_factor: a scaling factor for the coordinates of input morphologies, in order to make sure all coordinates are in the isotropic 1um CCFv3 space.

### Analyses [`./analyses`]
1. **3D visualization of the primary axonal tracts.** Users can plot the primary axonal tracts of each given neuron types or types in the 3D CCFv3 space, with the somas represented by red circles, termini of tracts by black triangles, and the tracts connecting them are randomly colored. It can be done by turning on the if-clause in line 110. Customized parameters are:
    - mpath_dir: the directory containing all axonal tracts. The same as the "out_dir" in the `Generation` section.
    - show_pts: number of points to show for each primary axonal tract. 
    - show_instances: maximal number of neurons to be plotted.
    - figname: output figname
    - scale: scale factor (divisor) to transforming the coordinates to millimeter scale.

2. **Along-path radius evolution of the axonal tracts.** Turn on the if-clause in line 237, and run the script. The parameters are similar to these in 3D visualization part. 

3. **Sub-dividing neurons into subtypes based on their termini locations.** Turn on the if-clause in line 474, and customize the neuron type (`ctx_type`) parameter. 

4. **Estimation of features for motif.**
    - We should firstly generate the cached file (`main_tract.pkl`, in picke format) by executing the script `generate_pkl.py`. 
    - Secondly, we can run the script `sdfeatures.py` to calculate the morphological features for the primary axonal tract (motif) of each neuron.
    An example of file file can be found at the directory `../data`.

5. **Calculate the diversity and stereotypy of neurons from various brain areas.**
    Based on the cached file `main_tract.pkl`, you can quantify the diversity and stereotypy of neurons from various brain areas, with the script `sdmatrix.py`. It currently support the quantification between different soma-types (s-types), projection-base subtypes, and lamination-based subtypes.
