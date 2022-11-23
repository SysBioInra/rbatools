# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
# package imports
from collections import OrderedDict as ordered_dict
from rba.core.constraint_blocks import ConstraintBlocks
from rbatools.information_block import InformationBlock


class EnzymeBlock(InformationBlock):
    """
    Class holding information on the enzymes in the model.

   Attributes
   ----------
   Elements : dict
       Each model-enzyme is represented by a key.
       The values, holding information on each enzyme, are dicts with predefined keys:
           'ID' : enzyme ID in model (type str)
           'OtherIDs' : identifiers of this enzyme in other namespaces (BiGG, KEGG ...) (type dict)
           'Reaction' : associated metabolic reaction (type str)
           'Isozymes' : Other enzymes catalysing the same metabolic reaction (type list)
           'IdenticalEnzymes' : Other enzymatic activities of this enzyme, cellular location ignored. (type list)
           'EnzymesWithIdenticalSubunitComposition' : Other enzymatic activities of this enzyme IN SAME COMPARTMENT. (type list)
           'Subunits' : Which proteins this enzyme is composed of and how many (type dict)
           'EnzymeCompartment' : Location of enzyme (type str)
           'ForwardCapacity_Constraint' : Id of associated forward capacity constraint.
           'BackwardCapacity_Constraint' : Id of associated backward capacity constraint.
    """

    def from_files(self, model, Info):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model

        Returns
        -------
        Dictionary with enzyme-info.

        """
        blocks = ConstraintBlocks(model)
        self.Elements = {}
        index = 0
        for i in range(len(blocks.enzymes.ids)):
            rx = _find_associated_enzymatic_reaction(i, blocks)
            enzymes = _get_enzyme_list(model)
            proteins = _get_protein_list(model)
            subunits = _find_subunits(i, blocks, model, enzymes)
            compartments = _determine_compartment(subunits, model, proteins)
            index += 1
            self.Elements[blocks.enzymes.ids[i]] = {'ID': blocks.enzymes.ids[i],
                                                    'OtherIDs': {},
                                                    'Reaction': rx,
                                                    'index': index,
                                                    'Isozymes': [],
                                                    'IdenticalEnzymes': numpy.nan,
                                                    'EnzymesWithIdenticalSubunitComposition': numpy.nan,
                                                    'Subunits': subunits,
                                                    'EnzymeCompartment': compartments,
                                                    'ForwardCapacity_Constraint':'',
                                                    'BackwardCapacity_Constraint':''}

        for i in self.Elements.keys():
            Subunits_Composition = ordered_dict(self.Elements[i]['Subunits'])
            identEnzs = []
            identCompoEnzs = []
            proto_subunits = []
            for su in list(Subunits_Composition.keys()):
                if '_' in su:
                    protoid = su.split('_')[0]
                else:
                    protoid = su
                proto_subunits.append(protoid)
            protoSU_dict = dict(zip(proto_subunits, list(Subunits_Composition.values())))

            for j in self.Elements.keys():
                if not i == j:
                    Subunits_Composition_other_enzyme = ordered_dict(self.Elements[j]['Subunits'])
                    if Subunits_Composition == Subunits_Composition_other_enzyme:
                        identEnzs.append(j)
                    proto_subunits_other_enzyme = []
                    for su in list(Subunits_Composition_other_enzyme.keys()):
                        protoid = su
                        if su.endswith(')'):
                            if '_(' in su:
                                protoid = su.split('_(')[0]
                        else:
                            if not su.startswith('average_protein'):
                                protoid = su.rsplit('_', 1)[0]
                            else:
                                protoid = su
                        proto_subunits_other_enzyme.append(protoid)
                    protoSU_dict_other_enzyme = dict(
                        zip(proto_subunits_other_enzyme, list(Subunits_Composition_other_enzyme.values())))
                    if protoSU_dict == protoSU_dict_other_enzyme:
                        identCompoEnzs.append(j)

            self.Elements[i]['IdenticalEnzymes'] = identEnzs
            self.Elements[i]['EnzymesWithIdenticalSubunitComposition'] = identCompoEnzs

    def overview(self):
        """
        Derive statistics on enzymes.

        Returns
        -------
        Dictionary with general numbers on enzymes.

        """
        dups = []
        unis = []
        for e in self.Elements.keys():
            if not e in dups:
                unis.append(e)
                if len(self.Elements[e]['IdenticalEnzymes']) > 0:
                    for eD in self.Elements[e]['IdenticalEnzymes']:
                        dups.append(eD)
        out = {'EnzymesTotal': len(self.Elements.keys()),
               'EnzymesUnique': len(unis)}
        return(out)


def _find_associated_enzymatic_reaction(enzyme_inQuestion, blocks):
    return blocks.enzymes.reaction_catalyzed[enzyme_inQuestion]


def _get_enzyme_list(model):
    out = []
    for e in model.enzymes.enzymes._elements:
        out.append(e.id)
    return(out)


def _get_protein_list(model):
    out = []
    for e in model.proteins.macromolecules._elements:
        out.append(e.id)
    return(out)


def _find_subunits(enzyme_inQuestion, blocks, model, enzymes):
    out = {}
    E = blocks.enzymes.ids[enzyme_inQuestion]
    ind = enzymes.index(E)
    if len(model.enzymes.enzymes._elements[ind].machinery_composition.reactants._elements) > 0:
        for p in model.enzymes.enzymes._elements[ind].machinery_composition.reactants._elements:
            out[p.__dict__['species']] = int(p.__dict__['stoichiometry'])
    return(out)


def _determine_compartment(subunits, model, enzymes):
    out = []
    if len(subunits.keys()) > 0:
        for s in subunits.keys():
            ind = enzymes.index(s)
            out.append(model.proteins.macromolecules._elements[ind].__dict__['compartment'])
    return(list(numpy.unique(out)))
