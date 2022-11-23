# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import json
import copy
import pandas
import numpy
from sbtab import SBtab
from rbatools.information_block import InformationBlock


class DataBlock(InformationBlock):
    """
    Class holding data from model-simulations.

    """

    def add_entries(self, Dict):
        for i in Dict.keys():
            self.Elements[i] = Dict[i]

    def from_dict(self, Dict):
        self.Elements = Dict

    def jsonize(self):
        Block = self.Elements
        block2 = copy.deepcopy(Block)
        for i in list(Block.keys()):
            if type(Block[i]) is dict:
                for j in list(Block[i].keys()):
                    block2[i][j] = json.dumps(Block[i][j], default=_json_int64_compensation)
            else:
                block2[i] = json.dumps(Block[i], default=_json_int64_compensation)
        return(block2)

    def to_data_frame(self, Col_list=None):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            if Col_list is None:
                fields = list(Block[list(Block.keys())[0]].keys())
            else:
                fields = Col_list
            TableOut = pandas.DataFrame(index=list(Block.keys()), columns=fields)
            for i in list(Block.keys()):
                for j in fields:
                    TableOut.loc[i, j] = Block[i][j]
            return TableOut
        else:
            return pandas.DataFrame()

    def to_sbtab(self, table_id, table_type, document_name=None, table_name=None, document=None, unit=None, Col_list=None, NameList=None):
        DF = self.to_data_frame(Col_list=Col_list)
        DF.reset_index(drop=False, inplace=True)
        DF.rename(columns={'index': 'VariableID'}, inplace=True)
        return(SBtab.SBtabTable.from_data_frame(df=DF, table_id=table_id, table_type=table_type, document_name=document_name, table_name=table_name, document=document, unit=unit, sbtab_version='1.0'))


def _json_int64_compensation(n):
    if isinstance(n, numpy.int64):
        return int(n)
    raise TypeError
