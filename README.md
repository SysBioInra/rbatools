# rbatools

rbatools serves as an interface to resource allocation modelling using Resource Balance Analysis (RBA) models (https://rba.inrae.fr).
It includes methods to solve, modify, and analyse RBA models and to export model information and simulation results into various formats.

## Installation

### Install this package locally
Navigate into the base directory of rbatools locally on your computer and execute:

    ```
    pip install .
    ```
### Install this from PyPI


## Repository structure

### rbatools
This is the library with all the necessary classes for the usage of RBA_tools.

### Instructional Jupyter notebooks
Please note that rba-models should be obtained (see section "Model availability") and stored in the same directory as this repository, prior to running the notebooks.
#### Example_workflows_rba_tools.ipynb
This is a Jupyter notebook with example applications and workflows of rbatools with B.subtilis model.
#### Example_workflows_rba_tools_Ecoli.ipynb
Same example applications and workflows with E.coli model.
#### Model_components_information_access.ipynb
This is a Jupyter notebook with example applications on how to access information on model-components and their relationships.

## Model availability
RBA models can be obtained from the repository https://github.com/SysBioInra/Bacterial-RBA-models

The Example_workflow_rba_tools notebook requires the previosly mentioned Bacterial-RBA-models repository to be placed in the same location, as this repository.
