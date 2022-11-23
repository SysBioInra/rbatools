# python 2/3 compatibility
from __future__ import division, print_function

import copy
import json
import pandas

# package imports
from rbatools.information_block import InformationBlock


class DescriptionBlock(InformationBlock):
    """
    Class holding general information on model.

    Author, metabolic reconstruction etc...

    Attributes
    ----------
    Elements : Dictionary information on model.


    """

    def __init__(self):
        self.Elements = {}

    def from_files(self, File):
        for i in File.index.tolist():
            self.Elements[i] = File.loc[i][1]

    def to_sbtab_compatible_data_frame(self, NameList=None, Col_list=None):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            if Col_list is None:
                fields = list(Block[list(Block.keys())[0]].keys())
            else:
                fields = Col_list
            if NameList is not None:
                if len(fields) == len(NameList):
                    colNames = NameList
                else:
                    colNames = fields
            else:
                colNames = fields

            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = i
                if isinstance(Block[i], str):
                    value_entry = Block[i].replace("'", "")
                Val = value_entry
                TableOut.loc[i, 'Measure'] = Var
                TableOut.loc[i, 'Value'] = Val
            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()
