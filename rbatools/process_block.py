# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import numpy
from rbatools.information_block import InformationBlock


class ProcessBlock(InformationBlock):
    """
    Class holding information on the processes in the model.

    Attributes
    ----------
    Elements : dict
        Each model-process is represented by a key.
        The values, holding information on each process, are dicts with predefined keys:
            'ID' : process ID in model (type str)
            'Name' : name of process (same as key) (type str)
            'Initiation' : reaction-string of initiation (type list)
            'Composition' : Which proteins the process-machinery is composed of and how many (type dict)
            'Components' : Substrates to process (type dict)
            'Capacity_Constraint' :  ID of associated capacity constraint
    """

    def from_files(self, model, Info):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model

        Returns
        -------
        Dictionary with process-info.

        """
        self.Elements = {}
        index = 0
        # for i in range(len(model.processes.processing_maps._elements)):
        for i in range(len(model.processes.processes._elements)):
            # if i <= len(model.processes.processing_maps._elements)-1:
            if i <= len(model.processes.processing_maps._elements)+1:
                pc = _get_process_info(model, i)
                index += 1
                self.Elements[model.processes.processes._elements[i].name.replace(' ', '_')] = {'ID': model.processes.processes._elements[i].id,
                                                                                                'Name': model.processes.processes._elements[i].name,
                                                                                                'Initiation': pc['Initiation'],
                                                                                                'index': index,
                                                                                                'Composition': pc['Machinery'],
                                                                                                'Components': pc['Components'],
                                                                                                'Capacity_Constraint':''}

    def overview(self):
        """
        Derive statistics on processes.

        Returns
        -------
        Dictionary with general numbers on processes.

        """
        nT = len(self.Elements.keys())
        out = {'ProcessesTotal': nT}
        return(out)


def _get_process_info(model, i):
    machinery = {}
    Processes = model.processes.processes._elements
    for j in range(len(Processes[i].machinery.machinery_composition.reactants._elements)):
        spec = Processes[i].machinery.machinery_composition.reactants._elements[j].species
        # int() added#
        stoc = int(
            Processes[i].machinery.machinery_composition.reactants._elements[j].stoichiometry)
        machinery[spec] = stoc

    respectiveMap = ''
    if len(Processes[i].processings.productions._elements) > 0:
        respectiveMap = Processes[i].processings.productions._elements[0].processing_map
    elif len(Processes[i].processings.degradations._elements) > 0:
        respectiveMap = Processes[i].processings.degradations._elements[0].processing_map

    ProcessingMaps = model.processes.processing_maps._elements
    mapIndex = -1
    for j in range(len(ProcessingMaps)):
        if str(ProcessingMaps[j].id) == str(respectiveMap):
            mapIndex = j
    if mapIndex > -1:
        ProcessingMaps[mapIndex].id

    compos = {}
    if mapIndex > -1:
        for component in ProcessingMaps[mapIndex].component_processings:
            ID = component.component
            Cost = component.machinery_cost
            Products = {k.species: k.stoichiometry for k in component.products._elements}
            Reactants = {k.species: k.stoichiometry for k in component.reactants._elements}
            compos.update({ID: {'Products': Products, 'Reactants': Reactants, 'Cost': Cost}})

    initiationString = ''
    if len(ProcessingMaps[mapIndex].constant_processing.reactants._elements) > 0:
        for j in range(len(ProcessingMaps[mapIndex].constant_processing.reactants._elements)):
            initiationString = initiationString + ' ' + \
                str(int(ProcessingMaps[mapIndex].constant_processing.reactants._elements[j].stoichiometry)) + ' ' + \
                ProcessingMaps[mapIndex].constant_processing.reactants._elements[j].species  # int() added#
    initiationString = initiationString + ' -->'
    if len(ProcessingMaps[mapIndex].constant_processing.products._elements) > 0:
        for j in range(len(ProcessingMaps[mapIndex].constant_processing.products._elements)):
            initiationString = initiationString + ' ' + \
                str(int(ProcessingMaps[mapIndex].constant_processing.products._elements[j].stoichiometry)) + ' ' + \
                ProcessingMaps[mapIndex].constant_processing.products._elements[j].species  # int() added#
    if initiationString == ' -->':
        initiationString = 'undefined'
    out = {'Machinery': machinery,
           'Components': compos,
           'Initiation': initiationString}
    return(out)
