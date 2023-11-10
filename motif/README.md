# motif
## Description
### generation/extract_main_path.py
This Python script is a data processing and analysis tool for handling neuron morphologies, particularly SWC (Standardized SWC) files. Here's a breakdown of the code:

1. The script begins by importing necessary libraries and modules, including `os`, `sys`, `glob`, `collections`, `numpy`, `pandas`, `matplotlib`, and custom-defined modules related to SWC file handling.

2. It defines two classes:
   - `MainPathValidator`: This class is intended for analyzing path lengths and creating histograms of path lengths for different classes of neuron morphologies.
   - `MainPathProjection`: This class is used for processing SWC files to extract the main axonal tract, pruning short segments, and saving the result as a new SWC file.

3. The `MainPathValidator` class has a method called `analyze_path_lengths` for reading a text file containing path length data, creating histograms, and saving the histograms as image files.

4. The `MainPathProjection` class has methods for processing SWC files:
   - `__init__`: Initializes the class with an SWC file, optional scaling, and a flip option to standardize the orientation.
   - `calc_seg_lengths`: Calculates segment lengths before pruning.
   - `extract_main_tract`: Extracts the main axonal tract while pruning short segments. It also saves the pruned SWC file.

5. A function called `wrapper` is defined, which wraps the process of creating a `MainPathProjection` object and extracting the main axonal tract for a given SWC file.

6. The main part of the script is wrapped in the `if __name__ == '__main__':` block, which ensures that the following code is executed when the script is run directly:
   - It defines input and output directories, a cell type file, and other parameters for processing SWC files.
   - It loads cell type information from the provided Excel file.
   - It creates a list of arguments for processing multiple SWC files in parallel using the `multiprocessing` library.
   - A pool of worker processes is created to process the SWC files concurrently.

Overall, the script is used to process a collection of SWC files representing neuron morphologies. It analyzes and extracts main axonal tracts from these morphologies and saves the pruned results in new SWC files, categorizing them by cell type. The code also provides the functionality to create histograms of path lengths for different neuron classes.
### plot_main_tracts_v2.py
This Python code is a collection of functions and scripts related to the analysis and visualization of neuronal data. Here's a brief description of what the code does:

1. It starts by importing various libraries including `sys`, `os`, `glob`, `cv2`, `numpy`, `random`, and several others. These libraries provide functionality for tasks like file system operations, image processing, numerical computations, data visualization, and more.

2. The code defines several functions:

   - `estimate_radius2d(pts, method=0)`: Estimates the radius of 2D points using different methods (mean distance, median distance, etc.).
   
   - `convert_to_proj_name(class_name)`: Takes a class name and converts it into a projection name.
   
   - `truncated_cone(ax, p0, p1, R0, R1, color, alpha=0.5)`: Plots a truncated cone in a 3D plot.
   
   - `calc_best_viewpoint(pts)`: Calculates the best viewpoint for a set of points based on Principal Component Analysis (PCA).

   - `plot_main_tracts(class_name, key='tract.swc')`: Plots primary axonal tracts for a given class name.

   - `calc_tract_radii1(class_name)`: Calculates the radii of tracts for a given class using one method.

   - `calc_tract_radii2(class_name)`: Calculates the radii of tracts for a given class using another method.

   - `plot_tracts_radii(ctypes, figname, radius_type=2)`: Plots the radii of tracts for specified classes.

   - `plot_tracts_radii_all(type_dict, figname, radius_type=2)`: Plots the radii of tracts for multiple classes.

   - `plot_CTX_ET_neurons(class_name)`: Plots neurons in 3D space and clusters them based on their terminal points.

3. The code also includes a section (controlled by the `if` statements) where different operations are performed based on the chosen condition:

   - first `if`: Plots main tracts for various classes.

   - second `if`: Estimates radii along main tracts and creates visualizations.

   - third `if`: Clusters neurons and divides them into sub-types.

Overall, this code is a tool for analyzing and visualizing neuronal data, particularly related to the estimation of radii along tracts and the clustering of neurons.
## Usage
### Generation
1. Just run the script "extract_main_path.py", after correctly set the parameters. The parameters are:
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