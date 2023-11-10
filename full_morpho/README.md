# full_morpho
## Description
### analyses/sdmatrix.py
This Python script is focused on loading and analyzing data related to neuronal cell types. Here's a breakdown of the code:

1. The script begins by importing necessary modules:
   - `os`: Provides a portable way to interact with the operating system.
   - `sys`: Provides access to some variables used or maintained by the interpreter and to functions that interact with the interpreter.
   - `itertools`: Provides functions to create iterators for efficient looping.
   - `numpy` (`np` alias): Provides support for arrays and matrices, as well as a large number of mathematical functions.
   - `pandas` (`pd` alias): Offers data structures for efficiently storing and manipulating data.
   
2. Paths to specific directories are appended to the system path.

3. Custom modules from `common_lib` and `microenviron/generation` are imported. These modules likely contain utility functions, configurations, and data structures used elsewhere in the code.

4. Several functions are defined to load and process data:

   - `load_data`: Loads and merges data from feature and cell type files, and filters out cell types with a number of neurons below a specified threshold.

   - `load_data_with_cortical_layer`: Extends `load_data` by also considering cortical layers. It creates a new column `cstype` combining soma region and cortical layer information.

   - `load_data_with_ptype`: Extends `load_data` by also considering specific target cell types. It creates a new column `ptype` combining soma region and target cell type information.

   - Each of these functions returns a Pandas DataFrame containing the selected data.

5. Functions for plotting are defined:

   - `plot_all`: Loads data and calculates correlations between features, then generates a structural diversity matrix based on soma region annotations.

   - `plot_struct`: Generates a structural diversity matrix based on soma region annotations. Additional options are available for customizing the plot.

   - `plot_struct_with_cortical_layer`: Similar to `plot_struct`, but also considers cortical layers.

   - `plot_struct_with_ptype`: Similar to `plot_struct`, but also considers specific target cell types.

6. The `if __name__ == '__main__':` block is used to execute the script when run directly. It specifies paths to the feature and cell type files, and an output directory for saving the generated plots. The script then proceeds to generate structural diversity matrices for different categories of cell types and save them as images.

Overall, this script is a tool for loading and analyzing data related to neuronal cell types, and for generating visualizations based on this data. It includes functionality for customizing the analysis based on soma region annotations, cortical layers, and specific target cell types.
### analyses/viz_swc_in_vtk_ccfbrain.py
This Python script involve the visualization of 3D neuronal structures using the VTK (Visualization Toolkit) library. Here's a breakdown of the code:

1. The script starts by importing necessary modules:
   - `pandas` (`pd` alias): Offers data structures for efficiently storing and manipulating data.
   - `numpy` (`np` alias): Provides support for arrays and matrices, as well as a large number of mathematical functions.
   - `os`: Provides a portable way to interact with the operating system.
   - `math`: Provides mathematical functions.
   - `argparse`: Allows the script to accept command-line arguments.
   - Various modules from VTK for 3D visualization.

2. The function `WriteImage` is defined, which takes a file name, render window, and a boolean `rgba` indicating whether to use an RGBA buffer. It writes the render window view to an image file with various supported formats (e.g., BMP, JPEG, PNM, PNG, PostScript, TIFF).

3. The script defines a command-line argument parser using `argparse`. It expects three optional arguments: `--v` (for a vtk file), `--s` (for an swc file), and `--folder` (for a folder containing swc files).

4. A VTK `vtkPolyDataReader` is initialized with the provided `vtk` file (`args.v`) and the data is stored in `polydata`.

5. VTK components (mapper, actor, renderer, render window, render window interactor) are set up for 3D visualization. A `vtkActor` is created with properties like color and opacity.

6. The `colordict` dictionary is defined, which maps numerical values to color names.

7. The script then enters a loop to process each SWC file found in the specified folder:
   - It reads the SWC file, extracts the relevant points, and creates lines connecting them to represent neuronal structures.
   - An `vtkActor` is created for each set of lines and added to the renderer.
   - The color of the lines is determined based on the `colordict`.

8. The size of the render window is set, and the window is rendered.

9. The function `WriteImage` is called to save the rendered view as an image file named 'render.png'.

10. The render window's interactor is started, allowing the user to interact with the 3D visualization.

Overall, this script reads in neuronal structure data from SWC files, visualizes them in 3D using VTK, and saves the rendered view as an image file.

## Usage
### analyses/sdmatrix.py
Use

	python sdmatrix.py
in shell or other IDEs to execute this file.

Users can also

Modify the variable *feat_file* to define the input feature file (e.g. feature.txt).

Modify the variable *celltype_file* to define the input celltype file (e.g. 41586_2021_3941_MOESM4_ESM.csv).

Modify the variable *outdir* in function to define the output directory.

The input feature data should be *.txt* format.

The input feature data should be *.csv* format.
### analyses/viz_swc_in_vtk_ccfbrain.py
Use

	python viz_swc_in_vtk_ccfbrain.py --v root.vtk --folder '../data/axon_sort'
in shell or other IDEs to execute this file.

The --folder is to define the output directory.