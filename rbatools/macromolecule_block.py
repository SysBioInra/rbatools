from __future__ import division, print_function

# global imports
import numpy
import pandas
# package imports
from rbatools.information_block import InformationBlock


class MacromoleculeBlock(InformationBlock):
    """
    Class holding information on the macromolecules, other than proteins in the model.

   Attributes
   ----------
   Elements : dict
       Each model-macromolecule is represented by a key.
       The values, holding information on each macromolecule, are dicts with predefined keys:
           'ID' : macromolecule ID in model (type str).
           'ProtoID' : location independent macromolecule ID in model (type str). Equals ID without compartment-specific ending.
           'Type' : DNA or RNA (type str)
           'Compartment' : Location of the macromolecule (type str)
           'Composition' : Which building-blocs the macromolecule is composed of and what is the stoichiometry (type dict)
           'ProcessRequirements' : Which processes the macromolecule requires for synthesis and maintenance and how much (type dict)
           'SupportsProcess' : Process-machineries, this macromolecule is a subunit of (type list)
           'AssociatedTarget' : Target associated with macromolecule (type str)
    """

    def from_files(self, model):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model
        Dataframe, holding uniprot-file information.
        Dataframe, holding ID-map information.

        Returns
        -------
        Dictionary with macromolecule-info.

        """
        DNAs = _get_dna_list(model)
        RNAs = _get_rna_list(model)
        self.Elements = {}
        index = 0
        for i in range(len(DNAs)):
            RequiredProcess = _get_process_requirements(model, DNAs[i], type='DNA')
            Composition = _get_macromolecule_composition(model, DNAs[i], type='DNA')
            protoid = DNAs[i]
            if DNAs[i].endswith(')'):
                if '_(' in DNAs[i]:
                    protoid = DNAs[i].split('_(')[0]
            else:
                protoid = DNAs[i]
            index += 1
            self.Elements[DNAs[i]] = {'ID': DNAs[i],
                                      'ProtoID': protoid,
                                      'Type': 'DNA',
                                      'index': index,
                                      'Compartment': _get_compartment(model, i, type='DNA'),
                                      'Composition': Composition['Composition'],
                                      'ProcessRequirements': _get_process_cost(model, Composition['Composition'], RequiredProcess),
                                      'SupportsProcess':  RequiredProcess['Part'],
                                      'AssociatedTarget':''}
        for i in range(len(RNAs)):
            RequiredProcess = _get_process_requirements(model, RNAs[i], type='RNA')
            Composition = _get_macromolecule_composition(model, RNAs[i], type='RNA')
            protoid = RNAs[i]
            if RNAs[i].endswith(')'):
                if '_(' in RNAs[i]:
                    protoid = RNAs[i].split('_(')[0]
            else:
                protoid = RNAs[i]
            index += 1
            self.Elements[RNAs[i]] = {'ID': RNAs[i],
                                      'ProtoID': protoid,
                                      'Type': 'RNA',
                                      'index': index,
                                      'Compartment': _get_compartment(model, i, type='RNA'),
                                      'Composition': Composition['Composition'],
                                      'ProcessRequirements': _get_process_cost(model, Composition['Composition'], RequiredProcess),
                                      'SupportsProcess':  RequiredProcess['Part'],
                                      'AssociatedTarget':''}

    def overview(self):
        """
        Derive statistics on macromolecules.

        Returns
        -------
        Dictionary with general numbers on macromolecules.

        """
        nRNA = len([i for i in self.Elements.keys() if self.Elements[i]['Type'] == 'RNA'])
        nDNA = len([i for i in self.Elements.keys() if self.Elements[i]['Type'] == 'DNA'])
        out = {'RNAsTotal': nRNA, 'DNAsTotal': nDNA}
        return(out)


def _get_rna_list(model):
    out = []
    for e in model.rnas.macromolecules._elements:
        out.append(e.id)
    return(out)


def _get_dna_list(model):
    out = []
    for e in model.dna.macromolecules._elements:
        out.append(e.id)
    return(out)


def _get_macromolecule_composition(model, i, type):
    if type == 'RNA':
        MacroMolecules = model.rnas.macromolecules
    elif type == 'DNA':
        MacroMolecules = model.dna.macromolecules
    out = {}
    numberNucleotides = 0
    composition = MacroMolecules.get_by_id(i).composition
    for base_number in range(len(composition)):
        out[composition[base_number].component] = round(
            composition[base_number].stoichiometry, 3)  # round(...,3) added#
        numberNucleotides += composition[base_number].stoichiometry
    out = {'Composition': out,
           'numberMonomers': int(numberNucleotides)}
    return(out)


def _get_compartment(model, i, type):
    if type == 'RNA':
        return(model.rnas.macromolecules._elements[i].__dict__['compartment'])
    elif type == 'DNA':
        return(model.dna.macromolecules._elements[i].__dict__['compartment'])


def _get_process_requirements(model, macromolecule, type):
    if type == 'RNA':
        macromoleculetype = 'rna'
    elif type == 'DNA':
        macromoleculetype = 'dna'
    out1 = []
    out2 = []
    Processes = model.processes.processes._elements
    for p in range(len(Processes)):
        for p2 in range(len(Processes[p].processings.productions._elements)):
            if Processes[p].processings.productions._elements[p2].set == macromoleculetype:
                for inp in range(len(Processes[p].processings.productions._elements[p2].inputs)):
                    if Processes[p].processings.productions._elements[p2].inputs._elements[inp].__dict__['species'] == macromolecule:
                        out1.append(Processes[p].name)
        for p3 in range(len(Processes[p].machinery.machinery_composition.reactants._elements)):
            if Processes[p].machinery.machinery_composition.reactants._elements[p3].species == macromolecule:
                out2.append(Processes[p].name)
    out = {'Req': out1,
           'Part': out2}
    return(out)


def _get_process_cost(model, Monomer, req):
    out = {}
    Processes = model.processes.processes._elements
    ProcessingMaps = model.processes.processing_maps._elements
    for Pref in req['Req']:
        cost = 0
        for p in range(len(Processes)):
            if Pref == Processes[p].name:
                Pmap = Processes[p].processings.productions._elements[0].processing_map
        for pm in range(len(ProcessingMaps)):
            Componentcost = 0
            if ProcessingMaps[pm].id == Pmap:
                for a in Monomer.keys():
                    for ap in range(len(ProcessingMaps[pm].component_processings._elements)):
                        if ProcessingMaps[pm].component_processings._elements[ap].component == a:
                            Componentcost = Monomer[a] * \
                                ProcessingMaps[pm].component_processings._elements[ap].machinery_cost
                            cost += Componentcost
            out[Pref] = round(cost, 3)  # round(...,3) added#
    return(out)
