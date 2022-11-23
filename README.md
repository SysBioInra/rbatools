# RBA_tools_WIP

RBAtools serves as an interface to resource allocation modelling using Resource Balance Analysis (RBA) models (https://rba.inrae.fr).
It includes methods to solve, modify, and analyse RBA models and to export model information and simulation results into various formats.

To use RBAtools, you need to install the rba library:
https://github.com/SysBioInra/RBApy

## Repository structure
### rbatools
This is the library with all the necessary classes for the usage of RBA_tools.
### html_documentation
This is a directory with html-files, serving as documentation of the rbatools library.
It can be viewed by navigating into this folder and opening the index.html file in a browser (via double click).
### Instructional Jupyter notebooks
#### Example_workflows_rba_tools.ipynb
This is a Jupyter notebook with example applications and workflows of rbatools with B.subtilis model.
#### Example_workflows_rba_tools_Ecoli.ipynb
Same example applications and workflows with E.coli model.
#### Model_components_information_access.ipynb
This is a Jupyter notebook with example applications on how to access information on model-components and their relationships.
#### Fig_rbatools_Manuscript.ipynb
This is a Jupyter notebook with the necessary code to generate the manuscript-figure.

## Model availability
RBA models can be obtained from the repository https://github.com/SysBioInra/Bacterial-RBA-models

The Example_workflow_rba_tools notebook requires the previosly mentioned Bacterial-RBA-models repository to be placed in the same location, as this repository.
