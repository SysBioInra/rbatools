# rbatools

rbatools is a programming interface for Resource Balance Analysis (RBA) models (https://rba.inrae.fr) based on the RBApy package (https://github.com/SysBioInra/RBApy).
It includes methods to solve, modify, and analyse RBA models and to export model information and simulation results into various formats.

## Usage remark

rbatools is based on linear-programming and therefore requires a linear-optimization package.
We strongly suggest using the (propriatary) CPLEX tool, which is available under a free academic license.

Nonetheless we included functionality with the freely available swiglpk solver, however this sometimes is not able to obtain satisfactory solutions.

## Installation

### Install this package locally:

After cloning this repository, in your terminal navigate into the base directory of rbatools locally on your computer and execute:

    pip install .

### Install this from PyPI:

    pip install rbatools

#### --------------------> Please note that installation of rbatools is required to run the tutorial.

## Tutorials
### MyBinder Tutorial
We generated a short tutorial, to be provided online with mybinder.

Please  follow this link:
https://mybinder.org/v2/gh/obodeit/rbatools_mybinder/main?labpath=tutorial_mybinder.ipynb 
    
### Instructional Jupyter notebooks
This repository includes various instructional jupyter notebooks, exemplifying the use of rbatools.

They are located in this repositories' subdirectory:

    /tutorials/jupyter_notebooks

#### --------------------><--------------------
#### We included a sample model in the tutorial folder and added a loader module, to import it into the notebook. However we only included the model of B.Subtilis. In order to run the other models, we suggest to obtain them from our repository. 
#### (see section "Model availability") and stored in the same directory as this repository.

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

We provided scripts with basic functionalities of rbatools in the subdirectory "scripts":

#### run_growth_rate_optimization.py

This script runs a growth-rate optimisation and stores the obtained results as sbtab.
The arguments to provide are an rba-model and the optional arguments for the used lp-solver (default: swiglpk, alternative: cplex) and the path where the results should be stored. To run it, execute the following command in your terminal:

    python run_growth_rate_optimization.py path/to/rba/model --lp-solver swiglpk --output-dir ...

#### generate_sbtab_of_model_for_html.py

This script imports an rba-model and exports its component-structure as tabular sbtab (as presented on https://rba.inrae.fr/models)
The arguments to provide are an rba-model and the optional argument for the path where the results should be stored.
To run it, execute the following command in your terminal:

    python generate_sbtab_of_model_for_html.py path/to/rba/model --output-dir ...

### run command-line tool

If the rbatools packaage is installed properly, the aforementioned scripts can also be executed as command-line tools
The respective commands (with identical arguments as scripts) for the terminal are:

    run-growth-rate-optimization path/to/rba/model --lp-solver swiglpk --output-dir ...
    
and 

    generate-sbtab-of-model-for-html path/to/rba/model --output-dir ...

## Documentation

A full api-reference documentation can be found at: https://sysbioinra.github.io/rbatools/api_ref

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
