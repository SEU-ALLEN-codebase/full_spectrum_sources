# bouton
## Description
The Python script "analyses/bouton_analyzer.py" is a program for analyzing features of boutons (synaptic terminals) in neuron morphologies. Here's a breakdown of the code:

The files "analyses/sdfeature.py" and "analyses/sdmatrix.py" are Python scripts for estimation of regional-level and neuronal-level bouton features for cross-scale feature map generation

The script "analyses/convert_eswc_to_swc.py" tries to convert the extended swc file (*ESWC*) that containing more informations, such as the bouton types to standardized *SWC* file.

## Usage
### analyses/bouton_analyzer.py
Use

	python bouton_analyzer.py

in shell or other IDEs to execute this file.

Users can also

Modify the variable *bouton_dir* to define the input directory.

Modify the variable *outfile* in function *calc\_overall\_features\_all* to define the output directory.

The input data should be *.swc* format.

