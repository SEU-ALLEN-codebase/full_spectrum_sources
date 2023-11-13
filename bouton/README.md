# bouton
## Description
The Python script "analyses/bouton_analyzer.py" is a program for analyzing features of boutons (synaptic terminals) in neuron morphologies. Here's a breakdown of the code:

1. The script begins by importing necessary modules:
   - `sys`: Provides access to some variables used or maintained by the interpreter and to functions that interact with the interpreter.
   - `os`: Provides a portable way to interact with the operating system.
   - `glob`: Allows for file pattern matching.
   - `numpy` (`np` alias): Provides support for arrays and matrices, as well as a large number of mathematical functions.
   - `pandas` (`pd` alias): Offers data structures for efficiently storing and manipulating data.
   - `seaborn` and other modules are imported for potential visualization purposes.

2. There are custom modules being imported, presumably containing functions (`load_image`, `parse_swc`, `scale_swc`, etc.) that are used later in the code.

3. The recursion limit for Python is set to a high value with `sys.setrecursionlimit(100000)`.

4. A class `BoutonFeatures` is defined. It appears to be designed for analyzing bouton features from a provided SWC (neuron morphology) file.

   - The `__init__` method loads a SWC file, creates a morphology object, and extracts information about axons, boutons, and their IDs.
   
   - There is a static method `load_swc` that loads and scales a SWC file.

   - Methods like `get_num_bouton`, `get_teb_ratio`, `bouton_density_by_axon`, `geodesic_distances`, `bouton_intervals`, `get_number_of_segments`, and `get_targets` calculate various bouton-related features.

   - The `calc_overall_features` method aggregates and prints several metrics based on bouton features.

5. The `single_processor` function takes a SWC file, processes its bouton features, and writes the results to an output file.

6. The `calc_overall_features_all` function processes a directory of SWC files for bouton analysis using multiple processes for parallelization.

7. The `if __name__ == '__main__':` block is used to execute the script when run directly, and it specifies a directory containing SWC files (`bouton_dir`) for analysis.

Overall, this script is a tool for automatically analyzing and extracting features from bouton structures in neuron morphologies.
## Usage
### analyses/bouton_analyzer.py
Use

	python bouton_analyzer.py

in shell or other IDEs to execute this file.

Users can also

Modify the variable *bouton_dir* to define the input directory.

Modify the variable *outfile* in function *calc\_overall\_features\_all* to define the output directory.

The input data should be *.swc* format.

