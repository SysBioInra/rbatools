# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import rbatools._auxiliary_functions as _auxiliary_functions

from rbatools.information_block import InformationBlock

class DensityConstraintBlock(InformationBlock):
    """
    Class holding information on the constraints regarding the compartment-densities in the model.

   Attributes
   ----------
   Elements : dict
       Each model metabolite-constraint is represented by a key.
       The values, holding information on each process-constraint, are dicts with predefined keys:
           'ID' : density constraint ID in model (type str)
           'AssociatedCompartment' : ID of compartment this constraint defines the capacity for (type str)
           'Type' : Equality or inequality (type str)
           'CapacityParameterID' :  ID of density parameter (type str)
           'CapacityParameter' : Definition of capacity parameter (type dict)
                See doc string of rbatools.rba_Session.RBA_Session.get_parameter_definition method
    """

    def from_files(self, model, Cs, matrix):
        index = 0
        self.Elements = {}
        for i in Cs['DensityConsts'].keys():
            index += 1
            if matrix.row_signs[Cs['DensityConsts'][i]] == 'L':
                cSign = '<='
            if matrix.row_signs[Cs['DensityConsts'][i]] == 'E':
                cSign = '='
            compartment = i.split('_density')[0]
            sizePar = _get_size_parameter(model, compartment)
            ParsedParam=_auxiliary_functions.return_parameter_definition(model=model,parameter=sizePar)
            Generic_Latex=""
            Generic_Latex=""
            self.Elements[i] = {'ID': i,
                                'AssociatedCompartment': compartment,
#                                'index': index,
                                'CapacityParameterID':sizePar,
                                'CapacityParameter': ParsedParam[sizePar],
                                'Generic parameter definition':ParsedParam[sizePar]["Generic_latex"],
                                'Specific parameter definition':ParsedParam[sizePar]["Specific_latex"],
                                'Type': cSign}


def _get_size_parameter(model, compartment):
    s = ''
    for c in model.density.target_densities._elements:
        if c.compartment==compartment:
            s=c.upper_bound
    return(s)
