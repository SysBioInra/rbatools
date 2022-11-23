# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import pandas
# package imports
from rbatools.information_block import InformationBlock


class ProteinBlock(InformationBlock):
    """
    Class holding information on the proteins in the model.

   Attributes
   ----------
   Elements : dict
       Each model-protein is represented by a key.
       The values, holding information on each protein, are dicts with predefined keys:
           'ID' : protein ID in model (type str). May be compartment_isoform.
           'ProtoID' : location independent protein ID in model (type str). Equals ID without compartment-specific ending.
           'ExternalIDs' : identifiers of this protein in other namespaces (Locus-tag, Gene-symbol, ...) (type dict)
           'Function' : associated function, according to Uniprot (type str)
           'Name' : name, according to Uniprot (type list)
           'associatedReactions' : Reactions, this protein is a subunit of (type list)
           'associatedEnzymes' : Enzymes, this protein is a subunit of (type list)
           'Compartment' : Location of the protein (type str)
           'AAcomposition' : Which amino-acids the protein is composed of and what is the stoichiometry(type dict)
           'AAnumber' : Length of protein (type int)
           'Weight' : Weight of protein (type float)
           'ProcessRequirements' : Which processes the protein requires for synthesis and maintenance and how much (type dict)
           'SupportsProcess' : Process-machineries, this protein is a subunit of (type list)
           'AssociatedTarget' : Wheter protein represents a (translation) target.
    """

    def from_files(self, model, IDmap, Info, UniprotData='Not There'):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model
        Dataframe, holding uniprot-file information.
        Dataframe, holding ID-map information.

        Returns
        -------
        Dictionary with protein-info.

        """
        Prots = _get_protein_list(model)
        if type(IDmap) is pandas.core.frame.DataFrame:
            IDmap = pandas.DataFrame(index=Prots)
        self.Elements = {}
        index = 0
        for i in range(len(Prots)):
            RequiredProcess = _get_process_requirements(model, Prots[i])
            Composition = _get_protein_composition(model, i)
            genes = {'ids': {},
                     'desc': ' ',
                     'name': ' ',
                     'weight': ' ',
                     'length': ' '}
            protoid = Prots[i]
            if Prots[i].endswith(')'):
                if '_(' in Prots[i]:
                    protoid = Prots[i].split('_(')[0]
            else:
                if not Prots[i].startswith('average_protein'):
                    protoid = Prots[i].rsplit('_', 1)[0]
                else:
                    protoid = Prots[i]
            if type(UniprotData) is pandas.core.frame.DataFrame:
                genes = _mine_uniprot_file(UniprotData, IDmap, protoid)
            index += 1
            self.Elements[Prots[i]] = {'ID': Prots[i],
                                       'ProtoID': protoid,
                                       'ExternalIDs': genes['ids'],
                                       'Function': genes['desc'],
                                       'Name': genes['name'],
                                       'associatedReactions': [],
                                       'index': index,
                                       'associatedEnzymes': [],
                                       'Compartment': _get_compartment(model, i),
                                       'AAcomposition': Composition['AAcomp'],
                                       'AAnumber': Composition['AAnum'],
                                       'Weight': genes['weight'],
                                       'ProcessRequirements': _get_process_cost(model, Composition['AAcomp'], RequiredProcess),
                                       'SupportsProcess':  RequiredProcess['Part'],
                                       'AssociatedTarget':''}

    def return_protein_iso_form_map(self):
        protoProteins = {}
        for protein in self.Elements.keys():
            if self.Elements[protein]['ProtoID'] in list(protoProteins.keys()):
                protoProteins[self.Elements[protein]['ProtoID']].append(protein)
            else:
                protoProteins[self.Elements[protein]['ProtoID']] = [protein]
        return(protoProteins)

    def overview(self):
        """
        Derive statistics on proteins.

        Returns
        -------
        Dictionary with general numbers on proteins.

        """
        nT = len(self.Elements.keys())
        out = {'ProteinsTotal': nT}
        return(out)


def _get_protein_list(model):
    out = []
    for e in model.proteins.macromolecules._elements:
        out.append(e.id)
    return(out)


def _mine_uniprot_file(UniprotData, IDmap, protein):
    differentIDs = {}
    function = ''
    name = ''
    length = numpy.nan
    mass = numpy.nan
    IDlist = [' ' + m + ' ' + str(n) + ' ' for m,
              n in zip(list(UniprotData['Entry']), list(UniprotData['Gene names']))]
    ProteinRowList = [i for i, s in enumerate(IDlist) if str(' ' + protein + ' ') in s]
    if len(ProteinRowList) > 0:
        ProteinRow = ProteinRowList[0]
        uniprotID = UniprotData['Entry'][ProteinRow]
        differentIDs.update({'UniprotID': uniprotID})
        if UniprotData['Length'][ProteinRow] is not numpy.nan:
            length = UniprotData['Length'][ProteinRow]
        if UniprotData['Mass'][ProteinRow] is not numpy.nan:
            mass = UniprotData['Mass'][ProteinRow]
        if UniprotData['EC number'][ProteinRow] is not numpy.nan:
            differentIDs.update({'ECnumber': 'EC '+str(UniprotData['EC number'][ProteinRow])})
        if UniprotData['Function [CC]'][ProteinRow] is not numpy.nan:
            function = UniprotData['Function [CC]'][ProteinRow]
        if UniprotData['Protein names'][ProteinRow] is not numpy.nan:
            name = UniprotData['Protein names'][ProteinRow]
        if IDmap != 'Not There':
            for j in range(IDmap.shape[1]):
                if '##()##' not in list(IDmap)[j]:
                    differentIDs[list(IDmap)[j]] = IDmap.loc[uniprotID, j]
                if '##()##' in list(IDmap)[j]:
                    differentIDs[list(IDmap)[j].split('##()##')[1]] = list(
                        IDmap)[j].split('##()##')[0] + '###' + IDmap.loc[uniprotID, j]

    out = {'ids': differentIDs,
           'desc': function,
           'name': name,
           'weight': mass,
           'length': length}

    return(out)


def _get_compartment(model, protein):
    return(model.proteins.macromolecules._elements[protein].__dict__['compartment'])


def _get_protein_composition(model, protein):
    out = {}
    numberAA = 0
    MacroMolecules = model.proteins.macromolecules._elements
    for a in range(len(MacroMolecules[protein].composition._elements)):
        out[MacroMolecules[protein].composition._elements[a].component] = int(
            round(MacroMolecules[protein].composition._elements[a].stoichiometry, 3))  # round(...,3) added#
        numberAA += MacroMolecules[protein].composition._elements[a].stoichiometry
    out = {'AAcomp': out,
           'AAnum': int(numberAA)}
    return(out)


def _get_process_requirements(model, protein):
    out1 = []
    out2 = []
    Processes = model.processes.processes._elements
    for p in range(len(Processes)):
        for p2 in range(len(Processes[p].processings.productions._elements)):
            if Processes[p].processings.productions._elements[p2].set == 'protein':
                for inp in range(len(Processes[p].processings.productions._elements[p2].inputs)):
                    if Processes[p].processings.productions._elements[p2].inputs._elements[inp].__dict__['species'] == protein:
                        out1.append(Processes[p].name)
        for p3 in range(len(Processes[p].machinery.machinery_composition.reactants._elements)):
            if Processes[p].machinery.machinery_composition.reactants._elements[p3].species == protein:
                out2.append(Processes[p].name)
    out = {'Req': out1,
           'Part': out2}
    return(out)


def _get_process_cost(model, AminoAcid, req):
    out = {}
    Processes = model.processes.processes._elements
    ProcessingMaps = model.processes.processing_maps._elements
    for Pref in req['Req']:
        cost = 0
        for p in range(len(Processes)):
            if Pref == Processes[p].name:
                Pmap = Processes[p].processings.productions._elements[0].processing_map
        for pm in range(len(ProcessingMaps)):
            AAcost = 0
            if ProcessingMaps[pm].id == Pmap:
                for a in AminoAcid.keys():
                    for ap in range(len(ProcessingMaps[pm].component_processings._elements)):
                        if ProcessingMaps[pm].component_processings._elements[ap].component == a:
                            AAcost = AminoAcid[a] * \
                                ProcessingMaps[pm].component_processings._elements[ap].machinery_cost
                            cost += AAcost
            out[Pref] = round(cost, 3)  # round(...,3) added#
    return(out)
