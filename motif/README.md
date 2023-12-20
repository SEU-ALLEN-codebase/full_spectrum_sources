# motif
## Description
This directory contains the source codes for extracting the primary axonal tract of each neuron, as well as evaluation of their diversity and stereotypy across the whole-brain. 

## Usage
### Generation
1. **Extraction of the primary axonal tracts from digitalized neuron morphologies in the format of SWC.** To do so, run the script `extract_main_path.py`, after correctly configuring the parameters. Parameters are:
    - swc_dir: directory containing all SWC files. The SWC files must be registered to the standardized CCFv3 atlas space, so that they be comparably analyzed in the subsequent steps. 
    - out_dir: output directory for saving the extracted axonal projection path files, in the format of SWC. 
    - ctype_file: a CSV file containing the manually curated neuron types, projection subtypes, and cortical laminations. 
    - scale_factor: a scaling factor for the coordinates of input morphologies, in order to make sure all coordinates are in the isotropic 1um CCFv3 space.

### Analysis
1. **3D visualization of the primary axonal tracts.** Users can plot the primary axonal tracts of each given neuron types or types in the 3D CCFv3 space, with the somas represented by red circles, termini of tracts by black triangles, and the tracts connecting them are randomly colored. It can be done by turning on the if-clause in line 110. Customized parameters are:
    - mpath_dir: the directory containing all axonal tracts. The same as the "out_dir" in the `Generation` section.
    - show_pts: number of points to show for each primary axonal tract. 
    - show_instances: maximal number of neurons to be plotted.
    - figname: output figname
    - scale: scale factor (divisor) to transforming the coordinates to millimeter scale.

2. Radius estimation along the axonal tracts
Just turn on the if-clause at line 236. 

3. Sub-dividing neurons into subtypes based on their termini locations
Just turn on the if-clause at line 473
