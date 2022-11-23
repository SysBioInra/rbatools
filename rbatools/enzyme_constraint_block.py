# python 2/3 compatibility
from __future__ import division, print_function

# package imports
from rbatools.information_block import InformationBlock
import rbatools._auxiliary_functions as _auxiliary_functions


class EnzymeConstraintBlock(InformationBlock):
    """
    Class holding information on the constraints regarding the enzymes in the model.

    Attributes
    ----------
    Elements : dict
        Each model enzyme-constraint is represented by a key.
        The values, holding information on each enzyme-constraint, are dicts with predefined keys:
           'ID' : enzyme-constraint ID in model (type str)
           'AssociatedEnzyme' : ID of enzyme this constraint relates to (type str)
           'AssociatedReaction' : ID of reaction this constraint relates to (type str)
           'Direction': forward or backward efficiency (type str)
           'Type' : Equality or inequality (type str)
           'CapacityParameterID' :  ID of capacity parameter (type str)
           'CapacityParameter' : Definition of capacity parameter (type dict)
                See doc string of rbatools.rba_Session.RBA_Session.get_parameter_definition method
    """

    def from_files(self, model, Cs, matrix):
        index = 0
        self.Elements = {}
        for i in Cs['EnzymeConsts'].keys():
            index += 1
            if matrix.row_signs[Cs['EnzymeConsts'][i]] == 'L':
                cSign = '<='
            if matrix.row_signs[Cs['EnzymeConsts'][i]] == 'E':
                cSign = '='
            rx = i.split('_enzyme_')[0]
            enzy = rx+'_enzyme'
            dire = i.split('_enzyme_')[1].split('_')[0]
            effPar = _get_efficiency_parameter(model, enzy, dire)
            ParsedParam=_auxiliary_functions.return_parameter_definition(model=model,parameter=effPar)
            self.Elements[i] = {'ID': i,
                                'AssociatedEnzyme': enzy,
                                'Direction': dire,
#                                'index': index,
                                'AssociatedReaction': rx,
                                'CapacityParameterID':effPar,
                                'CapacityParameter': ParsedParam[effPar],
                                'Generic parameter definition':ParsedParam[effPar]["Generic_latex"],
                                'Specific parameter definition':ParsedParam[effPar]["Specific_latex"],
                                'Type': cSign}


def _get_efficiency_parameter(model, enzy, direction):
    if direction == 'forward':
        o = model.enzymes.enzymes.__dict__['_elements_by_id'][enzy].__dict__['forward_efficiency']
        return(o)
    if direction == 'backward':
        o = model.enzymes.enzymes.__dict__['_elements_by_id'][enzy].__dict__['backward_efficiency']
        return(o)
