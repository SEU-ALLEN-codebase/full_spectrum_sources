# neuron_population
## Description
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