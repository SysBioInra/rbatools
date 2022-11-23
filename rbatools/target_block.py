# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
# package imports
import rbatools._auxiliary_functions as _auxiliary_functions
from rbatools.information_block import InformationBlock


class TargetBlock(InformationBlock):
    """
    Class holding information on the targets in the model.

   Attributes
   ----------
   Elements : dict
       Each model-target is represented by a key.
       The values, holding information on each target, are dicts with predefined keys:
           'ID' : Target ID in model (type str)
           'Group' : Class of targets (type str):
                Eg:
                'macrocomponent_production','maintenance_atp_target',
                'metabolite_production','replication_targets',
                'rna_degradation' or 'translation_targets'
           'Type' : Type of targets (type str):
                Eg:
                'metabolite_production', 'degradation_fluxes',
                'concentrations'  or 'reaction_fluxes'
            'TargetEntity' : For which entity the target is defined (type str)
            'TargetParameterID' : Model-parameter ID which defines target value (type str)
            'TargetParameter' : Model-parameter definition of target value (type dict)
            'TargetConstraint' : 'Value' if target is defined as explicit value or 'upper/lower Bound'.
        """

    def from_files(self, model):
        self.Elements = {}
        self.GroupList = []
        nGroups = len(list(model.targets.target_groups._elements))
        for i in range(nGroups):
            typeList = list(model.targets.target_groups._elements[i].__dict__.keys())
            typeList.remove('id')
            for j in typeList:
                for k in model.targets.target_groups._elements[i].__dict__[j].__dict__['_elements']:
                    self.GroupList.append(
                        str(model.targets.target_groups._elements[i].__dict__['id']))
                    spec = 'Other'
                    if 'species' in list(k.__dict__.keys()):
                        spec = str(k.__dict__['species'])
                    if 'reaction' in list(k.__dict__.keys()):
                        spec = str(k.__dict__['reaction'])
                    if str(j) == 'degradation_fluxes':
                        targettype_specifier = 'degradation'
                    elif str(j) == 'production_fluxes':
                        targettype_specifier = 'production'
                    elif str(j) == 'reaction_fluxes':
                        targettype_specifier = 'flux'
                    elif str(j) == 'concentrations':
                        targettype_specifier = 'concentration'
                    else:
                        targettype_specifier = ''
                    target_param=str(k.__dict__['value'])
                    if target_param != "None":
                        constraint_type="Value"
                    else:
                        target_param=str(k.__dict__['lower_bound'])
                        if target_param != "None":
                            constraint_type="LowerBound"
                        else:
                            target_param=str(k.__dict__['upper_bound'])
                            constraint_type="UpperBound"
                    ParsedParam=_auxiliary_functions.return_parameter_definition(model=model,parameter=target_param)
                    Edict = {'ID': 'Target_'+targettype_specifier+'_'+spec,
                             'Group': str(model.targets.target_groups._elements[i].__dict__['id']),
                             'Type': str(j),
                             'TargetEntity': spec,
                             'TargetParameterID': target_param,
                             'TargetParameter': ParsedParam[target_param],
                             'Generic parameter definition':ParsedParam[target_param]["Generic_latex"],
                             'Specific parameter definition':ParsedParam[target_param]["Specific_latex"],
                             'TargetConstraint': constraint_type
                             }
                    self.Elements.update({Edict['ID']: Edict})

    def overview(self):
        """
        Derive statistics on targets.

        Returns
        -------
        Dictionary with general numbers on targets.

        """
        out = {}
        for i in numpy.unique(self.GroupList):
            nI = len(list(numpy.where(numpy.array(self.GroupList) == i)[0]))
            out.update({'Targets_'+i: nI})
        return(out)
