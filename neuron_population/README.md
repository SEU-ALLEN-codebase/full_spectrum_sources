# neuron_population
## Description
### generation/calc_brain_statis.py
This Python script is a collection of functions and a main block of code for processing brain images. Here's a high-level description of what the code does:

1. It imports various libraries and modules for image processing and data manipulation.

2. It defines several functions for tasks such as reading image files, performing image processing operations (e.g., thresholding, filtering), and computing statistics on brain images.

3. The script defines a function `get_zeng_threshs` which reads a CSV file containing brain ID, current threshold, and supposed threshold, and returns a dictionary mapping brain IDs to threshold values.

4. It defines a function `get_filesize` that computes the file sizes of TeraFly block images, potentially saving 2D MIP images for debugging purposes.

5. It defines a function `get_block_counts` which calculates block counts based on provided dimensions and parameters.

6. The function `ada_thresholding` applies adaptive thresholding to an image.

7. The class `CalcBrainStatis` appears to be a brain image statistics calculator. It has methods for initializing the object, setting various parameters, and performing statistical calculations on brain images.

8. The script defines a function `brain_statis_wrapper` which serves as a wrapper for processing brain images. It loads a mask image, sets various parameters, and then performs statistical calculations on the brain image.

9. The main block of code at the bottom of the script sets up parameters and uses multiprocessing to process brain images in parallel. It specifies source directories and IDs, sets thresholds, and calls the `brain_statis_wrapper` function for each brain image.

Please note that the code relies on external modules and libraries, so it requires the appropriate environment to run successfully. Additionally, some functionalities may be specific to neuroscience or image processing applications.
## Usage
### Generation
In directory *neuron_population/generation*, use

	python calc_brain_statis.py
in shell or other IDEs to execute the generation code.

Users should

Modify the variable *tera\_downsize\_file* to define the input tera-downsize file.

Modify the variable *mask\_file\_dir* to define the input mask directory.

Modify the variable *out\_dir* to define the output directory.

The input tera-downsize file should be *.csv* format.

The input mask file should be *.v3draw* format.

Users should also

Define the variable *match\_str* and *tera\_path* to specify the directory and filename-format of the terafly data.
### Analyses

1. Download and install ParaView based on python3.9 from [https://www.paraview.org/download/](https://www.paraview.org/download/).
2. Add "ParaView-5.11.0-Windows-Python3.9-msvc2017-AMD64\bin\Lib\site-packages" to your *PYTHONPATH* or add the following code in paraview\_obj.py and paraview\_vtk.py.

		import sys
		sys.path.append(r'D:\ParaView\ParaView-5.11.0-Windows-Python3.9-msvc2017-AMD64\bin\Lib\site-packages') # use your install directory
3. Execute paraview\_obj.py and paraview\_vtk.py in shell or other IDEs and input the input data's directory.

		python paraview_obj.py
		D:/temp_need/paper/classic_region_soma

The input data in paraview\_obj.py should be *.obj* format.

The input data in paraview\_vtk.py should be *.vtk* format.