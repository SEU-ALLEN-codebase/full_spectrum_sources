# full spectrum sources
Source codes for the full-spectrum neuronal morphometry analysis project. This project is led by Prof. Hanchuan Peng at SEU, and it is a collaborative work across several groups worldwide, with major contributions from Yufeng Liu, Shengdian Jiang, Yingxin Li, Sujun Zhao, and source code refactoring conducted by Yufeng Liu and Qiaobo Gong. 

The project offers a comprehensive framework for extracting and analyzing neuronal morphometry across various scales, ranging from centimeters to micrometers. It also provides tools for cross-scale comparative and integrative analytical paradigm. The primary modules within the full spectrum sources repository encompass **motif**, **arbor**, **full morphology**, **microenviron**, **bouton**, and **neuron population**. Each module includes detailed documentation with step-by-step guidance. These resources aim to support researchers in their studies of neuronal morphometry.

## Contents

### Shared utilities
- **pylib**. This module is designed to offer fundamental utilities for processing of neuronal images and morphologies, such as parsing and saving 3D images and morphologies, as well as various morphological and topological analyses. Additionally, the project provides tools for mouse brain atlas analyses and visualization. For the most up-to-date versions of this module, please visit the following link: https://github.com/SEU-ALLEN-codebase/pylib.
- **common_lib**. Shared utilities and resources for other scripts within this project.


### Major modules
[motif](./motif/README.md)

[arbors](./arbors/README.md)

[full_morpho](./full_morpho/README.md)

[microenviron](./microenviron/README.md)

[bouton](./bouton/README.md)

[neuron_population](./neuron_population/README.md)

### Cross-scale analyses
[sd_features](./sd_features/README.md)

[sd_matrix](./sd_matrix/README.md)

## Installation
Most of the code is written in Python, and we recommend readers to first set up the environment using Anaconda (version 2023.09). Afterward, you can install or update several dependencies by:

	pip install -r requirements.txt

There are some customized visualization snippets using other standardalone third-party non-Python packages, including R (version 4.2.2), and ParaView (5.11.0) [https://www.paraview.org/download/](https://www.paraview.org/download/) to execute *paraview\_obj.py*,*paraview\_obj.py* in *neuron_population/analyses*.

## Contributing
Refer to our manuscript: [https://www.researchsquare.com/article/rs-3146034/v1](https://www.researchsquare.com/article/rs-3146034/v1).

## License
[MIT License](./LICENSE)

full spectrum sources is Copyright © 2023 SEU-ALLEN-codebase
