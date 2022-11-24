# rbatools

rbatools serves as an interface to resource allocation modelling using Resource Balance Analysis (RBA) models (https://rba.inrae.fr).
It includes methods to solve, modify, and analyse RBA models and to export model information and simulation results into various formats.

## Usage remark

rbatools is based on linear-programming and therefore requires a linear-optimization package.
We strongly suggest using the (propriatary) CPLEX tool, which is available under a free academic license.

Nonetheless we included functionality with the freely available swiglpk solver, however this sometimes is not able to obtain satisfactory solutions.

## Installation

### Install this package locally
After cloning this repository, in your terminal navigate into the base directory of rbatools locally on your computer and execute:

    pip install .

### Install this from PyPI (not available yet)
    pip install rbatools

#### --------------------> Please note that installation of rbatools is required to run the tutorial.

## Tutorials
### Instructional Jupyter notebooks
This repository includes various instructional jupyter notebooks, exemplifying the use of rbatools.

They are located in this repositories' subdirectory:

    /tutorials/jupyter_notebooks

#### --------------------> Please note that rba-models should be obtained (see section "Model availability") and stored in the same directory as this repository, prior to running the notebooks.

##### Requirements to run jupyter notebooks:
The jupyter notebooks use additional functionality provided by libraries, not defined as dependencies of rbatools. The user therefore should install those to be able to run the tutorial-notebooks. Please execute the following commands in your terminal:

    pip install notebook

    pip install ipywidgets

    pip install seaborn

##### Instructions to run jupyter notebooks:
In order to lauch the jupyter notebook interface, please execute the following command in your terminal:

    jupyter notebook

Once a file browser opens in your web browser, please navigate in to the "/tutorials/jupyter_notebooks" directory and doubleclick one of the following notebooks.

#### Example_workflows_rba_tools.ipynb
This is a Jupyter notebook with example applications and workflows of rbatools with B.subtilis model.
#### Model_components_information_access.ipynb

This is a Jupyter notebook with example applications on how to access information on model-components and their relationships in the B.subtilis model.

#### Example_workflows_rba_tools_Ecoli.ipynb
Example applications and workflows with E.coli model.
(Please note significantly longer computation times, than with B.subtilis model)


## Running

### scripts

we provided scripts with basic functionalities of rbatools in the subdirectory "scripts":

#### run_growth_rate_optimization.py

#### generate_sbtab_of_model_for_html.py

### run command-line tool

#### run-growth-rate-optimization

#### generate-sbtab-of-model-for-html

## Documentation

### HTML documentation

## Model availability
RBA models can be obtained from the repository https://github.com/SysBioInra/Bacterial-RBA-models

The Example_workflow_rba_tools notebook requires the previosly mentioned Bacterial-RBA-models repository to be placed in the same location, as this repository.

## Authors

Bodeit, O., Liebermeister, W. and Goelzer A.

## License

Copyright (c) 2022 INRAE - MaIAGE - France.

rbatools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rbatools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rbatools.  If not, see <https://www.gnu.org/licenses/>
