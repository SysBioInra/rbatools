# python 2/3 compatibility
from __future__ import division, print_function

# package imports
import rbatools._auxiliary_functions as _auxiliary_functions
from rbatools.information_block import InformationBlock


class ProcessConstraintBlock(InformationBlock):
    """
    Class holding information on the constraints regarding the processes in the model.

    Attributes
    ----------
    Elements : dict
        Each model process-constraint is represented by a key.
        The values, holding information on each process-constraint, are dicts with predefined keys:
            'ID' : process-constraint ID in model (type str)
            'AssociatedProcess' : Name of process this constraint relates to (type str)
            'AssociatedProcessID' : ID of process this constraint relates to (type str)
            'Type' : Equality or inequality (type str)
            'CapacityParameterID' :  ID of capacity parameter (type str)
            'CapacityParameter' : Definition of capacity parameter (type dict)
                See doc string of rbatools.rba_Session.RBA_Session.get_parameter_definition method
    """

    def from_files(self, model, Cs, matrix):
        self.Elements = {}
        index = 0
        for i in Cs['ProcessConsts'].keys():
            index += 1
            if matrix.row_signs[Cs['ProcessConsts'][i]] == 'L':
                cSign = '<='
            if matrix.row_signs[Cs['ProcessConsts'][i]] == 'E':
                cSign = '='
            effPar = _get_efficiency_parameter(model, i)
            if effPar is not None:
                ParsedParam=_auxiliary_functions.return_parameter_definition(model=model,parameter=effPar)
                self.Elements[i] = {'ID': i,
#                                    'index': index,
                                    'AssociatedProcess': None,
                                    'AssociatedProcessID': i.rsplit('_capacity')[0],
                                    'CapacityParameterID':effPar,
                                    'CapacityParameter': ParsedParam[effPar],
                                    'Generic parameter definition':ParsedParam[effPar]["Generic_latex"],
                                    'Specific parameter definition':ParsedParam[effPar]["Specific_latex"],
                                    'Type': cSign}
            else:
                self.Elements[i] = {'ID': i,
#                                    'index': index,
                                    'AssociatedProcess': None,
                                    'AssociatedProcessID': i.rsplit('_capacity')[0],
                                    'CapacityParameterID':"",
                                    'CapacityParameter': {},
                                    'Generic parameter definition':"",
                                    'Specific parameter definition':"",
                                    'Type': ""}

def _get_efficiency_parameter(model, process):
    x = model.processes.__dict__['processes'].__dict__['_elements_by_id']
    return(x[process.split('_capacity')[0]].__dict__['machinery'].__dict__['capacity'].__dict__['value'])
