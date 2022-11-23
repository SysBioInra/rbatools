# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import numpy
from rbatools.information_block import InformationBlock


class CompartmentBlock(InformationBlock):
    """
    Class holding information on the compartments in the model.

    Attributes
    ----------
    Elements : dict
        Each model-compartment is represented by a key.
        The values, holding information on each compartment, are dicts with predefined keys:
            'ID' : compartment ID in model (type str)
            'associatedProteins' : proteins localised to compartment (type list)
            'associatedEnzymes' : enzyme located in compartment (type list)
            'associatedReactions' : metabolic reactions located in compartment (type list)
            'associatedMacromolecules' : other macromolecules than proteins (DNA,RNA), located in compartment (type list)
            'Capacity_Constraint' : id of corresponding density constraint.
    """

    def from_files(self, model, Info):

        Compartments = _get_compartment_list(model)
        self.Elements = {}
        index = 0
        for i in Compartments:
            index += 1
            self.Elements[i] = {'ID': i,
                                'associatedProteins': [],
                                'index': index,
                                'associatedReactions': [],
                                'associatedEnzymes': [],
                                'associatedMacromolecules':[],
                                'Capacity_Constraint':''}

    def overview(self):
        """
        Derive statistics on compartments.

        Returns
        -------
        Dictionary with general numbers on compartments.

        """
        nT = len(self.Elements.keys())
        out = {'CompartmentsTotal': nT}
        return(out)


def _get_compartment_list(model):
    out = []
    for c in range(len(model.metabolism.compartments._elements)):
        out.append(model.metabolism.compartments._elements[c].id)
    return(out)
