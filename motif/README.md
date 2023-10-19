# motif
## Description
## Usage
### Generation
1. Just run the script extract_main_path.py, after correctly set the parameters. The parameters are:
    - swc_dir: directory containing all SWC files. The SWC files must be registered to the standardized CCFv3 atlas space, so that they be comparably analyzed in the subsequent steps
    - out_dir: output directory for the extracted axonal projection path files, in the format of SWC
    - min_files: number of minimal files for each cell type. Deprecated
    - scale_factor: the scale factor to multiple for the coordinates of input SWC nodes, in order to make sure all coordinates are in the 1um CCFv3 space.

### Analysis
1. Plotting the axonal tracts for different subtypes
Just turn on the if-clause at line 109 of plot_main_tracts_v2.py, and turn off other if-clause. The plotting involves the function "plot_main_tracts", which will plot primary tracts in 3D the CCFv3 space, with the somas represented by red circles, termini of tracts by black triangles, and the tracts connecting them are randomly colored.

2. Radius estimation along the axonal tracts
Just turn on the if-clause at line 236. 

3. Sub-dividing neurons into subtypes based on their termini locations
Just turn on the if-clause at line 473