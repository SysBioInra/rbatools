# python 2/3 compatibility
from __future__ import division, print_function

# package imports
from rbatools.information_block import InformationBlock


class MetaboliteConstraintBlock(InformationBlock):
    """
    Class holding information on the constraints regarding the metabolites in the model.

   Attributes
   ----------
   Elements : dict
       Each model metabolite-constraint is represented by a key.
       The values, holding information on each process-constraint, are dicts with predefined keys:
           'ID' : metabolite ID in model (type str)
           'AssociatedMetabolite' : ID of metabolite this constraint defines the mass-balance for (same as key) (type str)
           'Type' : Equality or inequality (type dict)
    """

    def from_files(self, Cs, matrix):
        index = 0
        self.Elements = {}
        for i in Cs['MetaboliteConsts'].keys():
            index += 1
            if matrix.row_signs[Cs['MetaboliteConsts'][i]] == 'L':
                cSign = '<='
            if matrix.row_signs[Cs['MetaboliteConsts'][i]] == 'E':
                cSign = '='
            self.Elements[i] = {'ID': i+'_mass_balance',
                                'AssociatedMetabolite': i,
                                'Type': cSign
                                }
