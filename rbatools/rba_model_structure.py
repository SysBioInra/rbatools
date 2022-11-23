# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import sys
import os.path
import numpy
import pandas
import copy
import json
import jxmlease
import warnings
import libsbml

from sbtab import SBtab
from rba import RbaModel , ConstraintMatrix
from rba.core.constraint_blocks import ConstraintBlocks

from rbatools.info_matrix import InfoMatrix
from rbatools.description_block import DescriptionBlock
from rbatools.metabolite_block import MetaboliteBlock
from rbatools.module_block import ModuleBlock
from rbatools.process_block import ProcessBlock
from rbatools.reaction_block import ReactionBlock
from rbatools.enzyme_block import EnzymeBlock
from rbatools.protein_block import ProteinBlock
from rbatools.macromolecule_block import MacromoleculeBlock
from rbatools.compartment_block import CompartmentBlock
from rbatools.metabolite_constraint_block import MetaboliteConstraintBlock
from rbatools.density_constraint_block import DensityConstraintBlock
from rbatools.process_constraint_block import ProcessConstraintBlock
from rbatools.enzyme_constraint_block import EnzymeConstraintBlock
from rbatools.statistics_block import StatisticsBlock
from rbatools.target_block import TargetBlock
from rbatools._auxiliary_functions import return_parameter_definition



class ModelStructureRBA(object):
    """
    Class holding information on model-structure.

    Attributes
    ----------
    GeneralInfo : rbatools.description_block.DescriptionBlock
         Model description
    MetaboliteInfo : rbatools.metabolite_block.MetaboliteBlock
         Metabolite information
    ModuleInfo : rbatools.module_block.ModuleBlock
         Module information
    ProcessInfo : rbatools.process_block.ProcessBlock
         Process information
    ReactionInfo : rbatools.reaction_block.ReactionBlock
         Reaction information
    EnzymeInfo : rbatools.enzyme_block.EnzymeBlock
         Enzyme information
    ProteinInfo : rbatools.protein_block.ProteinBlock
         Protein information
    MacromoleculeInfo : rbatools.macromolecule_block.MacromoleculeBlock
         Macromolecule (other than protein) information
    CompartmentInfo : rbatools.compartment_block.CompartmentBlock
         Compartment information
    MetaboliteConstraintsInfo : rbatools.metabolite_constraint_block.MetaboliteConstraintBlock
         Metabolite-constraint information
    DensityConstraintsInfo : rbatools.density_constraint_block.DensityConstraintBlock
         Density-constraint information
    ProcessConstraintsInfo : rbatools.process_constraint_block.ProcessConstraintBlock
         Process-constraint information
    EnzymeConstraintsInfo : rbatools.enzyme_constraint_block.EnzymeConstraintBlock
         Enzyme-constraint information
    ModelStatistics : rbatools.statistics_block.StatisticsBlock
         Statistics on Model
    TargetInfo : rbatools.target_block.TargetBlock
         Target information
    ProteinMatrix : numpy.array
        Matrix mapping proteins to consumers (enzymes and process-machineries)
    MediumDependencies : dict
        Dictionary with boundary metabolites as keys
        and a list of constraint_IDs, which are affected by its concentration, as values.
    MuDependencies : list
        List of constraint_IDs, which are affected by the growth-rate
    """

    def from_files(self, xml_dir:str):
        """
        Generates model-structure object from model-files and provided auxilliary information.

        Parameters
        ----------
        xml_dir : str
            Directory, where RBA-model is located
        """

        UniprotFile = _import_uniprot_file(xml_dir)
        GeneMap = _import_gene_annotations(xml_dir)
        Info = _import_model_info(xml_dir)
        SBMLfile = str('Not There')
        if Info['Value']['SBML-file'] != 'Not Provided':
            SBMLfile = _import_sbml_file(xml_dir, str(Info['Value']['SBML-file']))
        MetaboliteAnnotations = _import_metabolite_annotations(xml_dir)
        ReactionAnnotations = _import_reaction_annotations(xml_dir)

        model = RbaModel.from_xml(xml_dir)
        Zero_matrix = ConstraintMatrix(model)
        Zero_matrix.build_matrices(0)
        constraints = _sort_constraints(Zero_matrix, model)

        MetaboliteInfo = MetaboliteBlock()
        ModuleInfo = ModuleBlock()
        ProcessInfo = ProcessBlock()
        ReactionInfo = ReactionBlock()
        EnzymeInfo = EnzymeBlock()
        ProteinInfo = ProteinBlock()
        MacromoleculeInfo = MacromoleculeBlock()
        CompartmentInfo = CompartmentBlock()
        TargetInfo = TargetBlock()

        TargetInfo.from_files(model)
        MetaboliteInfo.from_files(model, Info, MetaboliteAnnotations, SBMLfile)
        ModuleInfo.from_files(model, SBMLfile)
        ProcessInfo.from_files(model, Info)
        ReactionInfo.from_files(model, Info, ReactionAnnotations, SBMLfile, MetaboliteInfo)
        EnzymeInfo.from_files(model, Info)
        ProteinInfo.from_files(model, GeneMap, Info, UniprotFile)
        MacromoleculeInfo.from_files(model)
        CompartmentInfo.from_files(model, Info)

        self.MetaboliteConstraintsInfo = MetaboliteConstraintBlock()
        self.DensityConstraintsInfo = DensityConstraintBlock()
        self.ProcessConstraintsInfo = ProcessConstraintBlock()
        self.EnzymeConstraintsInfo = EnzymeConstraintBlock()

        self.MetaboliteConstraintsInfo.from_files(constraints, Zero_matrix)
        self.DensityConstraintsInfo.from_files(model, constraints, Zero_matrix)
        self.ProcessConstraintsInfo.from_files(model, constraints, Zero_matrix)
        self.EnzymeConstraintsInfo.from_files(model, constraints, Zero_matrix)

        self.GeneralInfo = DescriptionBlock()
        self.GeneralInfo.from_files(Info)

        for constraint in self.MetaboliteConstraintsInfo.Elements.keys():
            MetaboliteInfo.Elements[self.MetaboliteConstraintsInfo.Elements[constraint]['AssociatedMetabolite']]['MassBalance_Constraint']=self.MetaboliteConstraintsInfo.Elements[constraint]['ID']
        for constraint in self.DensityConstraintsInfo.Elements.keys():
            CompartmentInfo.Elements[self.DensityConstraintsInfo.Elements[constraint]['AssociatedCompartment']]['Capacity_Constraint']=constraint
        for constraint in self.ProcessConstraintsInfo.Elements.keys():
            for i in list(ProcessInfo.Elements.keys()):
                if ProcessInfo.Elements[i]['ID']==self.ProcessConstraintsInfo.Elements[constraint]['AssociatedProcessID']:
                    ProcessInfo.Elements[i]['Capacity_Constraint']=constraint
                    self.ProcessConstraintsInfo.Elements[constraint]['AssociatedProcess']=i
        for constraint in self.EnzymeConstraintsInfo.Elements.keys():
            dir=self.EnzymeConstraintsInfo.Elements[constraint]['Direction']
            enz=self.EnzymeConstraintsInfo.Elements[constraint]['AssociatedEnzyme']
            if dir == 'forward':
                EnzymeInfo.Elements[enz]['ForwardCapacity_Constraint']=constraint
            elif dir == 'backward':
                EnzymeInfo.Elements[enz]['BackwardCapacity_Constraint']=constraint

        for target in TargetInfo.Elements.keys():
            target_species=TargetInfo.Elements[target]["TargetEntity"]
            if target_species in ProteinInfo.Elements.keys():
                ProteinInfo.Elements[target_species]["AssociatedTarget"]=target
            elif target_species in MetaboliteInfo.Elements.keys():
                MetaboliteInfo.Elements[target_species]["AssociatedTarget"]=target
            elif target_species in ReactionInfo.Elements.keys():
                ReactionInfo.Elements[target_species]["AssociatedTarget"]=target
            elif target_species in MacromoleculeInfo.Elements.keys():
                MacromoleculeInfo.Elements[target_species]["AssociatedTarget"]=target

        for protein in ProteinInfo.Elements.keys():
            AssociatedEnzyme = _find_associated_enzyme(EnzymeInfo.Elements, protein)
            ProteinInfo.Elements[protein]['associatedEnzymes'] = AssociatedEnzyme['Enz']
            ProteinInfo.Elements[protein]['associatedReactions'] = AssociatedEnzyme['Rx']

        CB = ConstraintBlocks(model)
        for enzyme in EnzymeInfo.Elements.keys():
            EnzymeInfo.Elements[enzyme]['Isozymes'] = _find_iso_enzymes(
                enzyme, CB, ReactionInfo.Elements, EnzymeInfo.Elements[enzyme]['Reaction'])

        for rx in ReactionInfo.Elements.keys():
            if ReactionInfo.Elements[rx]['Enzyme'] is not '':
                ReactionInfo.Elements[rx]['Compartment_Machinery'] = EnzymeInfo.Elements[ReactionInfo.Elements[rx]
                                                                                         ['Enzyme']]['EnzymeCompartment']

        for comp in CompartmentInfo.Elements.keys():
            ContainedMacromolecules = _find_contained_macromolecules(comp, MacromoleculeInfo.Elements)
            ContainedProteins = _find_contained_proteins(comp, ProteinInfo.Elements)
            ContainedEnzymes = _find_contained_enzymes(
                ContainedProteins, ProteinInfo.Elements, EnzymeInfo.Elements)
            CompartmentInfo.Elements[comp]['associatedMacromolecules'] = ContainedMacromolecules
            CompartmentInfo.Elements[comp]['associatedProteins'] = ContainedProteins
            CompartmentInfo.Elements[comp]['associatedReactions'] = ContainedEnzymes['containedReactions']
            CompartmentInfo.Elements[comp]['associatedEnzymes'] = ContainedEnzymes['containedEnzymes']

        self.ProcessInfo = ProcessInfo
        self.ReactionInfo = ReactionInfo
        self.EnzymeInfo = EnzymeInfo
        self.ProteinInfo = ProteinInfo
        self.MacromoleculeInfo = MacromoleculeInfo
        self.CompartmentInfo = CompartmentInfo
        self.MetaboliteInfo = MetaboliteInfo
        self.ModuleInfo = ModuleInfo
        self.TargetInfo = TargetInfo


        BioConstraintStats = _stats_constraints_biological(self.MetaboliteConstraintsInfo.Elements,
                                                        self.EnzymeConstraintsInfo.Elements,
                                                        self.DensityConstraintsInfo.Elements,
                                                        self.ProcessConstraintsInfo.Elements)

        MathConstraintStats = _stats_constraints_mathematical(self.MetaboliteConstraintsInfo.Elements,
                                                           self.EnzymeConstraintsInfo.Elements,
                                                           self.DensityConstraintsInfo.Elements,
                                                           self.ProcessConstraintsInfo.Elements,
                                                           self.EnzymeInfo.Elements,
                                                           self.ReactionInfo.Elements,
                                                           self.ProcessInfo.Elements)

        FullOverview = _generate_overview(self.ReactionInfo.overview(),
                                        self.MetaboliteInfo.overview(),
                                        self.ModuleInfo.overview(),
                                        self.EnzymeInfo.overview(),
                                        self.ProteinInfo.overview(),
                                        self.MacromoleculeInfo.overview(),
                                        self.ProcessInfo.overview(),
                                        self.CompartmentInfo.overview(),
                                        self.TargetInfo.overview(),
                                        BioConstraintStats,
                                        MathConstraintStats)

        self.ModelStatistics = StatisticsBlock()
        self.ModelStatistics.derive(FullOverview)

        self.ProteinMatrix = _generate_protein_matrix(self)
        self.ProteinGeneMatrix = _generate_protein_gene_matrix(self)
        self.MediumDependencies, self.MuDependencies = _find_parameter_dependencies(self,model=model)

    def from_json(self, inputString:str):
        """
        Generates model-structure object from model-structure in JSON-format.

        Parameters
        ----------
        inputString : json-str
            JSON-string to be parsed in to ModelStructure-object.
        """
        Block = json.loads(inputString)
        self.ModelStatistics = StatisticsBlock()
        self.GeneralInfo = DescriptionBlock()
        self.ProcessInfo = ProcessBlock()
        self.CompartmentInfo = CompartmentBlock()
        self.MetaboliteInfo = MetaboliteBlock()
        self.TargetInfo = TargetBlock()
        self.ModuleInfo = ModuleBlock()
        self.EnzymeInfo = EnzymeBlock()
        self.ProteinInfo = ProteinBlock()
        self.MacromoleculeInfo = MacromoleculeBlock()
        self.ReactionInfo = ReactionBlock()
        self.DensityConstraintsInfo = DensityConstraintBlock()
        self.ProcessConstraintsInfo = ProcessConstraintBlock()
        self.MetaboliteConstraintsInfo = MetaboliteConstraintBlock()
        self.EnzymeConstraintsInfo = EnzymeConstraintBlock()

        self.ModelStatistics.from_dict(Block['ModelStatistics'])
        self.GeneralInfo.from_dict(Block['ModelInformation'])
        self.ProcessInfo.from_dict(Block['Process'])
        self.CompartmentInfo.from_dict(Block['Compartment'])
        self.MetaboliteInfo.from_dict(Block['Metabolite'])
        self.ModuleInfo.from_dict(Block['Module'])
        self.EnzymeInfo.from_dict(Block['Enzyme'])
        self.ProteinInfo.from_dict(Block['Protein'])
        self.MacromoleculeInfo.from_dict(Block['Macromolecule'])
        self.ReactionInfo.from_dict(Block['Reaction'])
        self.DensityConstraintsInfo.from_dict(Block['DensityConstraint'])
        self.ProcessConstraintsInfo.from_dict(Block['ProcessConstraint'])
        self.MetaboliteConstraintsInfo.from_dict(Block['MetaboliteConstraint'])
        self.EnzymeConstraintsInfo.from_dict(Block['EnzymeConstraint'])
        self.TargetInfo.from_dict(Block['Target'])
        self.ProteinMatrix = Block['ProteinMatrix']
        self.ProteinMatrix['Matrix'] = numpy.array(self.ProteinMatrix['Matrix'])
        self.ProteinGeneMatrix = Block['ProteinGeneMatrix']
        self.ProteinGeneMatrix['Matrix'] = numpy.array(self.ProteinGeneMatrix['Matrix'])
        self.MediumDependencies = Block['MediumDependencies']
        self.MuDependencies = Block['MuDependencies']

    def export_json(self, path:str):
        """
        Saves model-structure object in JSON-format in specified directory, under the name ModelStructure.json.

        Parameters
        ----------
        path : str
            Directory, where to save JSON-file
        """
        Block = {'ModelInformation': self.GeneralInfo.Elements,
                 'ModelStatistics': self.ModelStatistics.Elements,
                 'Process': self.ProcessInfo.Elements,
                 'Compartment': self.CompartmentInfo.Elements,
                 'Metabolite': self.MetaboliteInfo.Elements,
                 'Target': self.TargetInfo.Elements,
                 'Module': self.ModuleInfo.Elements,
                 'Enzyme': self.EnzymeInfo.Elements,
                 'Protein': self.ProteinInfo.Elements,
                 'Macromolecule': self.MacromoleculeInfo.Elements,
                 'Reaction': self.ReactionInfo.Elements,
                 'DensityConstraint': self.DensityConstraintsInfo.Elements,
                 'ProcessConstraint': self.ProcessConstraintsInfo.Elements,
                 'MetaboliteConstraint': self.MetaboliteConstraintsInfo.Elements,
                 'EnzymeConstraint': self.EnzymeConstraintsInfo.Elements,
                 'ProteinMatrix': self.ProteinMatrix,
                 'ProteinGeneMatrix': self.ProteinGeneMatrix,
                 'MediumDependencies': self.MediumDependencies,
                 'MuDependencies': self.MuDependencies}
        Block['ProteinMatrix']['Matrix'] = Block['ProteinMatrix']['Matrix'].tolist()
        Block['ProteinGeneMatrix']['Matrix'] = Block['ProteinGeneMatrix']['Matrix'].tolist()
        JSONstring = json.dumps(Block, default=_json_int64_compensation)
        filename = path + '/ModelStructure.json'
        f = open(filename, 'w')
        f.write(JSONstring)
        f.close()
        return(JSONstring)

    def export_sbtab(self, add_links:bool=False, filename:str=""):
        """
        Saves model-structure object in SBtab-format under the specified name (filename).

        Parameters
        ----------
        filename : str
            Name, under which to save SBtab-file
        add_links : str
            Wheter to implement entry-format, which allows links between table-elements.
        """
        GeneralInfoTable = self.GeneralInfo.to_sbtab(table_id='ModelMetadata', table_type='Quantity', table_name='Model Information', Col_list=[
                                                    'Measure', 'Value'], NameList=['Element', 'Value'])
        GeneralInfoTable.filename = 'ModelMetadata.tsv'
        GeneralInfoTable.change_attribute('Text', 'Model metadata')
        try:
            GeneralInfoTable.unset_attribute('Date')
        except:
            pass
        GeneralInfoTable.unset_attribute('SBtabVersion')

        StatsTable = self.ModelStatistics.to_sbtab(
            table_id='ModelElements', table_type='Quantity', table_name='Model Size', Col_list=['Measure', 'Value'], NameList=['Element', 'Number'])
        StatsTable.filename = 'ModelElements.tsv'
        StatsTable.change_attribute('Text', 'Numbers of cell elements and constraints')
        try:
            StatsTable.unset_attribute('Date')
        except:
            pass
        StatsTable.unset_attribute('SBtabVersion')

        MetaboliteBlock_forChanges = copy.deepcopy(self.MetaboliteInfo)
        if add_links:
            for k in list(MetaboliteBlock_forChanges.Elements.keys()):
                if MetaboliteBlock_forChanges.Elements[k]['AssociatedTarget'] != '':
                    MetaboliteBlock_forChanges.Elements[k]['AssociatedTarget'] = '(!' + 'CellTarget' + \
                        '/'+MetaboliteBlock_forChanges.Elements[k]['AssociatedTarget']+'!)'
                if MetaboliteBlock_forChanges.Elements[k]['MassBalance_Constraint'] != '':
                    MetaboliteBlock_forChanges.Elements[k]['MassBalance_Constraint'] = '(!' + 'MetaboliteConstraint' + \
                        '/'+MetaboliteBlock_forChanges.Elements[k]['MassBalance_Constraint']+'!)'
                if MetaboliteBlock_forChanges.Elements[k]['Compartment'] in self.CompartmentInfo.Elements.keys():
                    MetaboliteBlock_forChanges.Elements[k]['Compartment'] = '(!' + 'Compartment' + \
                        '/'+MetaboliteBlock_forChanges.Elements[k]['Compartment']+'!)'
                for index, id in enumerate(MetaboliteBlock_forChanges.Elements[k]['ReactionsInvolvedWith']):
                    MetaboliteBlock_forChanges.Elements[k]['ReactionsInvolvedWith'][
                        index] = '(!' + 'Reaction'+'/'+id+'!)'
                for id in list(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    if not pandas.isna(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id]):
                        MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id] = str('(!identifiers:'+id+'/'+str(
                            MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id])+'|'+str(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id])+'!)')

        MetaboliteTable = MetaboliteBlock_forChanges.to_sbtab(table_id='Compound', table_type='Quantity', table_name='Metabolites', Col_list=[
            'ID', 'Name', 'Compartment', 'Type', 'ReactionsInvolvedWith', 'OtherIDs','AssociatedTarget','MassBalance_Constraint'], NameList=['ID', 'Name', 'Compartment', 'Type', 'Reactions', 'Annotation','IsCellTarget','MassBalance'])
        MetaboliteTable.filename = 'Compound.tsv'
        MetaboliteTable.change_attribute(
            'Text', 'Metabolite species are localised in cell compartments and are associated with metabolite mass balance constraints.')
        try:
            MetaboliteTable.unset_attribute('Date')
        except:
            pass
        MetaboliteTable.unset_attribute('SBtabVersion')

        ReactionBlock_forChanges = copy.deepcopy(self.ReactionInfo)
        if add_links:
            for k in list(ReactionBlock_forChanges.Elements.keys()):
                if ReactionBlock_forChanges.Elements[k]['AssociatedTarget'] != '':
                    ReactionBlock_forChanges.Elements[k]['AssociatedTarget'] = '(!' + 'CellTarget' + \
                        '/'+ReactionBlock_forChanges.Elements[k]['AssociatedTarget']+'!)'
                oldEnz = ReactionBlock_forChanges.Elements[k]['Enzyme']
                if len(oldEnz) > 0:
                    ReactionBlock_forChanges.Elements[k]['Enzyme'] = '(!' + \
                        'Enzyme'+'/'+oldEnz+'!)'

                replacements = {}
                for i in ReactionBlock_forChanges.Elements[k]['Formula'].split(' '):
                    if i in self.MetaboliteInfo.Elements.keys():
                        replacements[i] = '(!'+'Compound' + '/'+i+'!)'
                for i in replacements.keys():
                    ReactionBlock_forChanges.Elements[k]['Formula'] = ReactionBlock_forChanges.Elements[k]['Formula'].replace(
                        str(' '+i+' '), str(' '+replacements[i]+' '))

                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Compartment_Species']):
                    if id in self.CompartmentInfo.Elements.keys():
                        ReactionBlock_forChanges.Elements[k]['Compartment_Species'][
                            index] = '(!' + 'Compartment'+'/'+id+'!)'

                for id in list(ReactionBlock_forChanges.Elements[k]['Reactants'].keys()):
                    ReactionBlock_forChanges.Elements[k]['Reactants']['(!'+'Compound' +
                                                                      '/'+id+'!)'] = ReactionBlock_forChanges.Elements[k]['Reactants'].pop(id)
                for id in list(ReactionBlock_forChanges.Elements[k]['Products'].keys()):
                    ReactionBlock_forChanges.Elements[k]['Products']['(!'+'Compound' +
                                                                     '/'+id+'!)'] = ReactionBlock_forChanges.Elements[k]['Products'].pop(id)
                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Twins']):
                    ReactionBlock_forChanges.Elements[k]['Twins'][index] = '(!' + \
                        'Reaction'+'/'+id+'!)'
                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Compartment_Machinery']):
                    ReactionBlock_forChanges.Elements[k]['Compartment_Machinery'][
                        index] = '(!' + 'Compartment'+'/'+id+'!)'
                if 'ProtoID' in list(ReactionBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    ReactionBlock_forChanges.Elements[k]['OtherIDs'].pop('ProtoID')
                for id in list(ReactionBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    if not pandas.isna(ReactionBlock_forChanges.Elements[k]['OtherIDs'][id]):
                        ReactionBlock_forChanges.Elements[k]['OtherIDs'][id] = str('(!identifiers:'+id+'/'+str(
                            ReactionBlock_forChanges.Elements[k]['OtherIDs'][id])+'|'+str(ReactionBlock_forChanges.Elements[k]['OtherIDs'][id])+'!)')

        ReactionTable = ReactionBlock_forChanges.to_sbtab(table_id='Reaction', table_type='Quantity', table_name='Reactions', Col_list=['ID', 'Name', 'Type', 'Reversible', 'Formula', 'Enzyme', 'Compartment_Machinery', 'Twins', 'Compartment_Species','AssociatedTarget', 'OtherIDs'], NameList=[
            'ID', 'Name', 'Type', 'IsReversible', 'ReactionFormula', 'Enzyme', 'EnzymeCompartment', 'IsoenzymeReactions', 'ReactionCompartment','IsCellTarget', 'Annotation'])
        ReactionTable.filename = 'Reaction.tsv'
        ReactionTable.change_attribute(
            'Text', 'Chemical reactions are localised in cell compartments. All reactions are enzyme catalysed.')
        try:
            ReactionTable.unset_attribute('Date')
        except:
            pass
        ReactionTable.unset_attribute('SBtabVersion')

        EnzymeBlock_forChanges = copy.deepcopy(self.EnzymeInfo)
        if add_links:
            for k in list(EnzymeBlock_forChanges.Elements.keys()):
                if EnzymeBlock_forChanges.Elements[k]['ForwardCapacity_Constraint'] != '':
                    EnzymeBlock_forChanges.Elements[k]['ForwardCapacity_Constraint'] = '(!' + 'EnzymeCapacityConstraint' + \
                        '/'+EnzymeBlock_forChanges.Elements[k]['ForwardCapacity_Constraint']+'!)'
                if EnzymeBlock_forChanges.Elements[k]['BackwardCapacity_Constraint'] != '':
                    EnzymeBlock_forChanges.Elements[k]['BackwardCapacity_Constraint'] = '(!' + 'EnzymeCapacityConstraint' + \
                        '/'+EnzymeBlock_forChanges.Elements[k]['BackwardCapacity_Constraint']+'!)'
                oldRx = EnzymeBlock_forChanges.Elements[k]['Reaction']
                EnzymeBlock_forChanges.Elements[k]['Reaction'] = '(!'+'Reaction'+'/'+oldRx+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['Isozymes']):
                    EnzymeBlock_forChanges.Elements[k]['Isozymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['IdenticalEnzymes']):
                    EnzymeBlock_forChanges.Elements[k]['IdenticalEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['EnzymeCompartment']):
                    EnzymeBlock_forChanges.Elements[k]['EnzymeCompartment'][index] = '(!' + \
                        'Compartment'+'/' + id+'!)'
                for id in list(EnzymeBlock_forChanges.Elements[k]['Subunits'].keys()):
                    EnzymeBlock_forChanges.Elements[k]['Subunits']['(!'+'Protein'+'/' +
                                                                   id+'!)'] = EnzymeBlock_forChanges.Elements[k]['Subunits'].pop(id)

        EnzymeTable = EnzymeBlock_forChanges.to_sbtab(table_id='Enzyme', table_type='Quantity', table_name='Enzymes', Col_list=[
            'ID', 'Reaction', 'IdenticalEnzymes', 'Subunits', 'EnzymeCompartment', 'OtherIDs', 'Isozymes','ForwardCapacity_Constraint','BackwardCapacity_Constraint'], NameList=['ID', 'CatalysedReaction', 'OtherEnzymaticActivities', 'Subunits', 'Compartment', 'Annotation', 'Isoenzymes','ForwardCapacity','BackwardCapacity'])
        EnzymeTable.filename = 'Enzyme.tsv'
        EnzymeTable.change_attribute(
            'Text', 'Enzymes are localised in cell compartments and catalyse specific reactions. To describe multi-functional enzyme, RBA uses multiple enzyme entries. Enzymes are associated with enzyme capacity constraints.')
        try:
            EnzymeTable.unset_attribute('Date')
        except:
            pass
        EnzymeTable.unset_attribute('SBtabVersion')

        ProteinBlock_forChanges = copy.deepcopy(self.ProteinInfo)
        if add_links:
            for k in list(ProteinBlock_forChanges.Elements.keys()):
                if ProteinBlock_forChanges.Elements[k]['AssociatedTarget'] != '':
                    ProteinBlock_forChanges.Elements[k]['AssociatedTarget'] = '(!' + 'CellTarget' + \
                        '/'+ProteinBlock_forChanges.Elements[k]['AssociatedTarget']+'!)'
                oldComp = ProteinBlock_forChanges.Elements[k]['Compartment']
                ProteinBlock_forChanges.Elements[k]['Compartment'] = '(!' + \
                    'Compartment'+'/'+oldComp+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['associatedReactions']):
                    ProteinBlock_forChanges.Elements[k]['associatedReactions'][
                        index] = '(!' + 'Reaction'+'/'+id+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['associatedEnzymes']):
                    ProteinBlock_forChanges.Elements[k]['associatedEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['SupportsProcess']):
                    ProteinBlock_forChanges.Elements[k]['SupportsProcess'][index] = '(!' + \
                        'Process'+'/'+id+'!)'
                for id in list(ProteinBlock_forChanges.Elements[k]['ProcessRequirements'].keys()):
                    ProteinBlock_forChanges.Elements[k]['ProcessRequirements'][
                        '(!'+'Process'+'/' + id+'!)'] = ProteinBlock_forChanges.Elements[k]['ProcessRequirements'].pop(id)
                for id in list(ProteinBlock_forChanges.Elements[k]['ExternalIDs'].keys()):
                    if not pandas.isna(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id]):
                        if id == 'UniprotID':
                            ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] = str('(!identifiers:uniprot/'+str(
                                ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'|'+str(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'!)')
                        if id == 'ECnumber':
                            if ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] != 'nan':
                                ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] = str('(!identifiers:ec-code/'+str(
                                    ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id]).split('EC ')[1]+'|'+str(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'!)')

        ProteinTable = ProteinBlock_forChanges.to_sbtab(table_id='Protein', table_type='Quantity', table_name='Proteins', Col_list=['ID', 'ProtoID', 'Name', 'Compartment', 'SupportsProcess', 'associatedEnzymes', 'associatedReactions', 'AAnumber', 'ProcessRequirements','AssociatedTarget', 'Weight', 'ExternalIDs', 'Function'], NameList=[
            'ID', 'ProtoID', 'Name', 'Compartment', 'ContributesToProcess', 'EnzymaticActivity', 'CatalysedReaction', 'ChainLength', 'RequiresProcess','IsCellTarget', 'MolecularWeight', 'Annotation', 'Function'])
        ProteinTable.filename = 'Protein.tsv'
        ProteinTable.change_attribute(
            'Text', 'Proteins are localised in cell compartments and can be (or be part of) metabolic enzymes.')
        try:
            ProteinTable.unset_attribute('Date')
        except:
            pass
        ProteinTable.unset_attribute('SBtabVersion')

        MacromoleculeBlock_forChanges = copy.deepcopy(self.MacromoleculeInfo)
        if add_links:
            for k in list(MacromoleculeBlock_forChanges.Elements.keys()):
                if MacromoleculeBlock_forChanges.Elements[k]['AssociatedTarget'] != '':
                    MacromoleculeBlock_forChanges.Elements[k]['AssociatedTarget'] = '(!' + 'CellTarget' + '/'+MacromoleculeBlock_forChanges.Elements[k]['AssociatedTarget']+'!)'
                oldComp = MacromoleculeBlock_forChanges.Elements[k]['Compartment']
                MacromoleculeBlock_forChanges.Elements[k]['Compartment'] = '(!' + \
                    'Compartment'+'/'+oldComp+'!)'
                for index, id in enumerate(MacromoleculeBlock_forChanges.Elements[k]['SupportsProcess']):
                    MacromoleculeBlock_forChanges.Elements[k][
                        'SupportsProcess'][index] = '(!' + 'Process'+'/'+id+'!)'
                for id in list(MacromoleculeBlock_forChanges.Elements[k]['ProcessRequirements'].keys()):
                    MacromoleculeBlock_forChanges.Elements[k]['ProcessRequirements'][
                        '(!'+'Process'+'/' + id+'!)'] = MacromoleculeBlock_forChanges.Elements[k]['ProcessRequirements'].pop(id)

        MacromoleculeTable = MacromoleculeBlock_forChanges.to_sbtab(table_id='Macromolecule', table_type='Quantity', table_name='Macromolecules', Col_list=[
                                                                   'ID', 'ProtoID', 'Type', 'Compartment', 'SupportsProcess', 'ProcessRequirements','AssociatedTarget'], NameList=['ID', 'ProtoID', 'Type', 'Compartment', 'SupportedtedProcess', 'RequiresProcess','IsCellTarget'])
        MacromoleculeTable.filename = 'Macromolecules.tsv'
        MacromoleculeTable.change_attribute(
            'Text', 'Macromolecule are localised in cell compartments and can be part of cellular machinery or serve another function inside the cell.')
        try:
            MacromoleculeTable.unset_attribute('Date')
        except:
            pass
        MacromoleculeTable.unset_attribute('SBtabVersion')

        CompartmentBlock_forChanges = copy.deepcopy(self.CompartmentInfo)
        if add_links:
            for k in list(CompartmentBlock_forChanges.Elements.keys()):
                if CompartmentBlock_forChanges.Elements[k]['Capacity_Constraint'] != '':
                    CompartmentBlock_forChanges.Elements[k]['Capacity_Constraint'] = '(!' + 'DensityConstraint' + \
                        '/'+CompartmentBlock_forChanges.Elements[k]['Capacity_Constraint']+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedReactions']):
                    CompartmentBlock_forChanges.Elements[k]['associatedReactions'][
                        index] = '(!' + 'Reaction'+'/'+id+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedEnzymes']):
                    CompartmentBlock_forChanges.Elements[k]['associatedEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedProteins']):
                    CompartmentBlock_forChanges.Elements[k]['associatedProteins'][
                        index] = '(!' + 'Protein'+'/'+id+'!)'

        CompartmentTable = CompartmentBlock_forChanges.to_sbtab(table_id='Compartment', table_type='Quantity', table_name='Cell Compartments', Col_list=[
                                                               'ID', 'associatedProteins', 'associatedEnzymes', 'associatedReactions','Capacity_Constraint'], NameList=['ID', 'Proteins', 'Enzymes', 'Reactions','Density'])

        CompartmentTable.filename = 'Compartment.tsv'
        CompartmentTable.change_attribute(
            'Text', 'Cell compartments are used to describe the localisation of proteins, enzymes, and reactions and are associated with density constraints.')
        try:
            CompartmentTable.unset_attribute('Date')
        except:
            pass
        CompartmentTable.unset_attribute('SBtabVersion')

        ModuleBlock_forChanges = copy.deepcopy(self.ModuleInfo)
        if add_links:
            for k in list(ModuleBlock_forChanges.Elements.keys()):
                for index, id in enumerate(ModuleBlock_forChanges.Elements[k]['Members']):
                    if id in self.ReactionInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Reaction'+'/'+id+'!)'
                    if id in self.ProteinInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Protein'+'/'+id+'!)'
                    if id in self.MetaboliteInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Compound'+'/'+id+'!)'
                    if id in self.CompartmentInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Compartment'+'/'+id+'!)'

        ModuleTable = ModuleBlock_forChanges.to_sbtab(table_id='CellModule', table_type='Quantity', table_name='Cell Modules', Col_list=[
                                                     'ID', 'Name', 'Members'], NameList=['ID', 'Name', 'Contains'])
        ModuleTable.filename = 'Module.tsv'
        ModuleTable.change_attribute('Text', 'Information on Modules in RBAmodel')
        try:
            ModuleTable.unset_attribute('Date')
        except:
            pass
        ModuleTable.unset_attribute('SBtabVersion')

        ProcessBlock_forChanges = copy.deepcopy(self.ProcessInfo)
        if add_links:
            for k in list(ProcessBlock_forChanges.Elements.keys()):
                if ProcessBlock_forChanges.Elements[k]['Capacity_Constraint'] != '':
                    ProcessBlock_forChanges.Elements[k]['Capacity_Constraint'] = '(!' + 'MachineryCapacityConstraint' + \
                        '/'+ProcessBlock_forChanges.Elements[k]['Capacity_Constraint']+'!)'
                for component in ProcessBlock_forChanges.Elements[k]['Components'].keys():
                    ProductKeys = list(
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Products'].keys())
                    for species in ProductKeys:
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Products'][
                            '(!'+'Compound'+'/'+species+'!)'] = ProcessBlock_forChanges.Elements[k]['Components'][component]['Products'].pop(species)
                    ReactantKeys = list(
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Reactants'].keys())
                    for species in ReactantKeys:
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Reactants'][
                            '(!'+'Compound'+'/'+species+'!)'] = ProcessBlock_forChanges.Elements[k]['Components'][component]['Reactants'].pop(species)
                for id in list(ProcessBlock_forChanges.Elements[k]['Composition'].keys()):
                    if id in self.ProteinInfo.Elements.keys():
                        ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Protein'+'/' +
                                                                           id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)
                    elif id in self.MacromoleculeInfo.Elements.keys():
                        ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Macromolecule' +
                                                                           '/' + id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)
                    elif id in self.MetaboliteInfo.Elements.keys():
                        ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Metabolite'+'/' +
                                                                           id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)

        ProcessTable = ProcessBlock_forChanges.to_sbtab(table_id='Process', table_type='Quantity', table_name='Macromolecular Processes', Col_list=[
            'ID', 'Name', 'Composition', 'Components', 'Initiation','Capacity_Constraint'], NameList=['ID', 'Name', 'MachineSubunits', 'MachineComponents', 'InitiationCofactors','Capacity'])
        ProcessTable.filename = 'Process.tsv'
        ProcessTable.change_attribute(
            'Text', 'Macromolecular machines catalyse the biochemical reactions that produce, modify, and degrade macromolecules. They are catalysed by macromolecular machines and associated with process capacity constraints.')
        try:
            ProcessTable.unset_attribute('Date')
        except:
            pass
        ProcessTable.unset_attribute('SBtabVersion')

        MetaboliteConstraintsBlock_forChanges = copy.deepcopy(self.MetaboliteConstraintsInfo)
        if add_links:
            for k in list(MetaboliteConstraintsBlock_forChanges.Elements.keys()):
                oldMet = MetaboliteConstraintsBlock_forChanges.Elements[k]['AssociatedMetabolite']
                MetaboliteConstraintsBlock_forChanges.Elements[k][
                    'AssociatedMetabolite'] = '(!' + 'Compound'+'/'+oldMet+'!)'

        MetaboliteConstraintTable = MetaboliteConstraintsBlock_forChanges.to_sbtab(table_id='MetaboliteConstraint', table_type='Quantity', table_name='Metabolite Mass Balance Constraints', Col_list=[
                                                                                  'ID', 'AssociatedMetabolite', 'Type'], NameList=['ID', 'Metabolite', 'Type'])
        MetaboliteConstraintTable.filename = 'MetaboliteConstraint.tsv'
        MetaboliteConstraintTable.change_attribute(
            'Text', 'Metabolite mass balance constraints ensure mass balance and stationary fluxes in metabolism.')
        try:
            MetaboliteConstraintTable.unset_attribute('Date')
        except:
            pass
        MetaboliteConstraintTable.unset_attribute('SBtabVersion')

        CompartmentConstraintsBlock_forChanges = copy.deepcopy(self.DensityConstraintsInfo)
        if add_links:
            for k in list(CompartmentConstraintsBlock_forChanges.Elements.keys()):
                oldComp = CompartmentConstraintsBlock_forChanges.Elements[k]['AssociatedCompartment']
                CompartmentConstraintsBlock_forChanges.Elements[k]['AssociatedCompartment'] = '(!' + 'Compartment'+'/'+oldComp+'!)'

        DensityConstraintTable = CompartmentConstraintsBlock_forChanges.to_sbtab(table_id='DensityConstraint', table_type='Quantity', table_name='Density Constraints', Col_list=['ID', 'AssociatedCompartment', 'Type','CapacityParameterID', 'Generic parameter definition','Specific parameter definition'], NameList=['ID', 'Compartment', 'Type', 'CapacityParameter','Formula', 'Formula (parameterized)'])
        DensityConstraintTable.filename = 'DensityConstraint.tsv'
        DensityConstraintTable.change_attribute(
            'Text', 'Density constraints put an upper bound on the sum of macromolecule concentrations in a given compartment. The capacity parameter defines this bound in units corresponding to amino acids (contained in proteins), or one third of nucleotides (contained in RNA).')
        try:
            DensityConstraintTable.unset_attribute('Date')
        except:
            pass
        DensityConstraintTable.unset_attribute('SBtabVersion')

        ProcessConstraintsBlock_forChanges = copy.deepcopy(self.ProcessConstraintsInfo)
        if add_links:
            for k in list(ProcessConstraintsBlock_forChanges.Elements.keys()):
                oldComp = ProcessConstraintsBlock_forChanges.Elements[k]['AssociatedProcess']
                ProcessConstraintsBlock_forChanges.Elements[k]['AssociatedProcess'] = '(!' + \
                    'Process'+'/'+oldComp+'!)'

        ProcessConstraintTable = ProcessConstraintsBlock_forChanges.to_sbtab(table_id='MachineryCapacityConstraint', table_type='Quantity', table_name='Machinery Capacity Constraints', Col_list=[
            'ID', 'AssociatedProcess', 'Type','CapacityParameterID', 'Generic parameter definition','Specific parameter definition'], NameList=['ID', 'Process', 'Type', 'CapacityParameter', 'Formula', 'Formula (parameterized)'])
        ProcessConstraintTable.filename = 'MachineryConstraint.tsv'
        ProcessConstraintTable.change_attribute(
            'Text', 'A machinery capacity constraint states that the rate of a macromolecular process is proportional to the concentration of the catalysing machine. The proportionality constant  (capacity parameter) can be a constant or a function of model parameters such as the growth rate.')
        try:
            ProcessConstraintTable.unset_attribute('Date')
        except:
            pass
        ProcessConstraintTable.unset_attribute('SBtabVersion')

        EnzymeConstraintsBlock_forChanges = copy.deepcopy(self.EnzymeConstraintsInfo)
        if add_links:
            for k in list(EnzymeConstraintsBlock_forChanges.Elements.keys()):
                oldComp = EnzymeConstraintsBlock_forChanges.Elements[k]['AssociatedEnzyme']
                EnzymeConstraintsBlock_forChanges.Elements[k]['AssociatedEnzyme'] = '(!' + \
                    'Enzyme'+'/'+oldComp+'!)'

        EnzymeConstraintTable = EnzymeConstraintsBlock_forChanges.to_sbtab(table_id='EnzymeCapacityConstraint', table_type='Quantity', table_name='Enzyme Capacity Constraints', Col_list=[
            'ID', 'AssociatedEnzyme', 'AssociatedReaction', 'Direction', 'Type','CapacityParameterID' ,'Generic parameter definition','Specific parameter definition'], NameList=['ID', 'Enzyme', 'Reaction', 'Direction', 'Type', 'CapacityParameter', 'Formula', 'Formula (parameterized)'])
        EnzymeConstraintTable.filename = 'EnzymeConstraint.tsv'
        EnzymeConstraintTable.change_attribute(
            'Text', 'An enzyme capacity constraint states that a reaction rate is proportional to the concentration of the catalysing enzyme. The proportionality constant (capacity parameter) can be a constant or a function of model parameters such as the growth rate.')
        try:
            EnzymeConstraintTable.unset_attribute('Date')
        except:
            pass
        EnzymeConstraintTable.unset_attribute('SBtabVersion')

        TargetBlock_forChanges = copy.deepcopy(self.TargetInfo)
        if add_links:
            for k in list(TargetBlock_forChanges.Elements.keys()):
                oldTargSpec = TargetBlock_forChanges.Elements[k]['TargetEntity']
                if oldTargSpec in list(self.MetaboliteInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Compound'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.ReactionInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Reaction'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.ProteinInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Protein'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.MacromoleculeInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Macromolecule'+'/'+oldTargSpec+'!)'

        TargetTable = TargetBlock_forChanges.to_sbtab(table_id='CellTarget', table_type='Quantity', table_name='Cell Targets', Col_list=[
            'ID', 'Group', 'Type', 'TargetEntity', 'TargetParameterID', 'Generic parameter definition','Specific parameter definition'], NameList=['ID', 'TargetGroup', 'TargetType', 'TargetEntity', 'TargetParameter', 'Formula', 'Formula (parameterized)'])
        TargetTable.filename = 'Target.tsv'
        TargetTable.change_attribute(
            'Text', 'Cell targets are biochemical variables (fluxes or concentrations) that need to remain below or above a certain threshold value to ensure a viable cell. They define constraints in the RBA model.')
        try:
            TargetTable.unset_attribute('Date')
        except:
            pass
        TargetTable.unset_attribute('SBtabVersion')

        if filename != "":
            filename_SBtab = filename
        else:
            if add_links:
                filename_SBtab = 'RBA_model_withLinks'
            else:
                filename_SBtab = 'RBA_model'

        Out = SBtab.SBtabDocument(name='rbatools_withLinks', sbtab_init=None,
                                  filename=str(filename_SBtab+'.tsv'))
        Out.add_sbtab(GeneralInfoTable)
        Out.add_sbtab(StatsTable)
        Out.add_sbtab(MetaboliteTable)
        Out.add_sbtab(ReactionTable)
        Out.add_sbtab(EnzymeTable)
        Out.add_sbtab(ProteinTable)
        Out.add_sbtab(MacromoleculeTable)
        Out.add_sbtab(CompartmentTable)
        Out.add_sbtab(ModuleTable)
        Out.add_sbtab(ProcessTable)
        Out.add_sbtab(TargetTable)
        Out.add_sbtab(MetaboliteConstraintTable)
        Out.add_sbtab(DensityConstraintTable)
        Out.add_sbtab(ProcessConstraintTable)
        Out.add_sbtab(EnzymeConstraintTable)

        Out.name = filename_SBtab
        Out.change_attribute('DocumentName', self.GeneralInfo.Elements['Name'])
        Out.change_attribute('DocumentType', 'rba-model-synopsis')
        Out.write()


def _find_parameter_dependencies(ModelStructure,model):
    MedDepts = {}
    MuDepts = []
    for tA in ModelStructure.TargetInfo.Elements.keys():
        if len(list(ModelStructure.TargetInfo.Elements[tA]['TargetParameter'].keys()))>0:
            if ModelStructure.TargetInfo.Elements[tA]['TargetParameter']['Type'] not in ['constant']:
                if ModelStructure.TargetInfo.Elements[tA]['TargetParameter']['Type'] not in ['Aggregate']:
                    for var in ModelStructure.TargetInfo.Elements[tA]['TargetParameter']['Variables']:
                        if var == 'growth_rate':
                            MuDepts.append(tA)
                        else :
                            if var in MedDepts.keys():
                                MedDepts[var].append(tA)
                            else:
                                MedDepts[var]=[tA]
                else:
                    for term in ModelStructure.TargetInfo.Elements[tA]['TargetParameter']['Multiplicative Terms']:
                        term_definition=return_parameter_definition(model=model,parameter=term)
                        for var in term_definition[term]['Variables']:
                            if var == 'growth_rate':
                                if tA not in MuDepts:
                                    MuDepts.append(tA)
                            else :
                                if var in MedDepts.keys():
                                    if tA not in MedDepts[var]:
                                        MedDepts[var].append(tA)
                                else:
                                    MedDepts[var]=[tA]


    for eK in ModelStructure.EnzymeConstraintsInfo.Elements.keys():
        if len(list(ModelStructure.EnzymeConstraintsInfo.Elements[eK]['CapacityParameter'].keys()))>0:
            if ModelStructure.EnzymeConstraintsInfo.Elements[eK]['CapacityParameter']['Type'] not in ['constant']:
                if ModelStructure.EnzymeConstraintsInfo.Elements[eK]['CapacityParameter']['Type'] not in ['Aggregate']:
                    for var in ModelStructure.EnzymeConstraintsInfo.Elements[eK]['CapacityParameter']['Variables']:
                        if var == 'growth_rate':
                            MuDepts.append(eK)
                        else :
                            if var in MedDepts.keys():
                                MedDepts[var].append(eK)
                            else:
                                MedDepts[var]=[eK]
                else:
                    for term in ModelStructure.EnzymeConstraintsInfo.Elements[eK]['CapacityParameter']['Multiplicative Terms']:
                        term_definition=return_parameter_definition(model=model,parameter=term)
                        for var in term_definition[term]['Variables']:
                            if var == 'growth_rate':
                                if eK not in MuDepts:
                                    MuDepts.append(eK)
                            else :
                                if var in MedDepts.keys():
                                    if eK not in MedDepts[var]:
                                        MedDepts[var].append(eK)
                                else:
                                    MedDepts[var]=[eK]
    for dK in ModelStructure.DensityConstraintsInfo.Elements.keys():
        if len(list(ModelStructure.DensityConstraintsInfo.Elements[dK]['CapacityParameter'].keys()))>0:
            if ModelStructure.DensityConstraintsInfo.Elements[dK]['CapacityParameter']['Type'] not in ['constant']:
                if ModelStructure.DensityConstraintsInfo.Elements[dK]['CapacityParameter']['Type'] not in ['Aggregate']:
                    for var in ModelStructure.DensityConstraintsInfo.Elements[dK]['CapacityParameter']['Variables']:
                        if var == 'growth_rate':
                            MuDepts.append(dK)
                        else :
                            if var in MedDepts.keys():
                                MedDepts[var].append(dK)
                            else:
                                MedDepts[var]=[dK]
                else:
                    for term in ModelStructure.DensityConstraintsInfo.Elements[dK]['CapacityParameter']['Multiplicative Terms']:
                        term_definition=return_parameter_definition(model=model,parameter=term)
                        for var in term_definition[term]['Variables']:
                            if var == 'growth_rate':
                                if dK not in MuDepts:
                                    MuDepts.append(dK)
                            else :
                                if var in MedDepts.keys():
                                    if dK not in MedDepts[var]:
                                        MedDepts[var].append(dK)
                                else:
                                    MedDepts[var]=[dK]
    for pK in ModelStructure.ProcessConstraintsInfo.Elements.keys():
        if len(list(ModelStructure.ProcessConstraintsInfo.Elements[pK]['CapacityParameter'].keys()))>0:
            if ModelStructure.ProcessConstraintsInfo.Elements[pK]['CapacityParameter']['Type'] not in ['constant']:
                if ModelStructure.ProcessConstraintsInfo.Elements[pK]['CapacityParameter']['Type'] not in ['Aggregate']:
                    for var in ModelStructure.ProcessConstraintsInfo.Elements[pK]['CapacityParameter']['Variables']:
                        if var == 'growth_rate':
                            MuDepts.append(pK)
                        else :
                            if var in MedDepts.keys():
                                MedDepts[var].append(pK)
                            else:
                                MedDepts[var]=[pK]
                else:
                    for term in ModelStructure.ProcessConstraintsInfo.Elements[pK]['CapacityParameter']['Multiplicative Terms']:
                        term_definition=return_parameter_definition(model=model,parameter=term)
                        if term_definition[term]['Type'] not in ['constant']:
                            for var in term_definition[term]['Variables']:
                                if var == 'growth_rate':
                                    if pK not in MuDepts:
                                        MuDepts.append(pK)
                                else :
                                    if var in MedDepts.keys():
                                        if pK not in MedDepts[var]:
                                            MedDepts[var].append(pK)
                                    else:
                                        MedDepts[var]=[pK]

    return(MedDepts, MuDepts)


def _generate_protein_gene_matrix(ModelStructure):
    uniqueProteins = []
    uniqueProteinMap = {}
    Proteins = list(ModelStructure.ProteinInfo.Elements.keys())
    for i in ModelStructure.ProteinInfo.Elements.keys():
        if ModelStructure.ProteinInfo.Elements[i]['ProtoID'] not in list(uniqueProteinMap.keys()):
            uniqueProteinMap.update({ModelStructure.ProteinInfo.Elements[i]['ProtoID']: []})
            uniqueProteins.append(ModelStructure.ProteinInfo.Elements[i]['ProtoID'])
        uniqueProteinMap[ModelStructure.ProteinInfo.Elements[i]['ProtoID']].append(i)
    ProteinProteinMatrix = numpy.zeros((len(list(uniqueProteinMap.keys())), len(list(ModelStructure.ProteinInfo.Elements.keys()))))
    for u in list(uniqueProteinMap.keys()):
        row_ind = uniqueProteins.index(u)
        for i in uniqueProteinMap[u]:
            col_ind = Proteins.index(i)
            ProteinProteinMatrix[row_ind, col_ind] = 1
    return({'Matrix': numpy.array(ProteinProteinMatrix), 'Proteins': Proteins, 'ProtoProteins': uniqueProteins})


def _generate_protein_matrix(ModelStructure):
    Proteins = list(ModelStructure.ProteinInfo.Elements.keys())
    Processes = [ModelStructure.ProcessInfo.Elements[i]['ID'] +
                 '_machinery' for i in list(ModelStructure.ProcessInfo.Elements.keys())]
    Enzymes = list(ModelStructure.EnzymeInfo.Elements.keys())
    Consumers = list(set(list(Enzymes+Processes)))
    ProteinMatrix = numpy.zeros((len(Proteins), len(Consumers)))
    for p in Proteins:
        if len(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']) > 0:
            for pc in list(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']):
                coeff = 0
                row_ind = Proteins.index(p)
                col_ind = Consumers.index(
                    ModelStructure.ProcessInfo.Elements[pc]['ID']+'_machinery')
                coeff = ModelStructure.ProcessInfo.Elements[pc]['Composition'][p]
                ProteinMatrix[row_ind, col_ind] += coeff
        if len(ModelStructure.ProteinInfo.Elements[p]['associatedEnzymes']) > 0:
            for ez in list(ModelStructure.ProteinInfo.Elements[p]['associatedEnzymes']):
                coeff = 0
                row_ind = Proteins.index(p)
                col_ind = Consumers.index(ez)
                coeff = ModelStructure.EnzymeInfo.Elements[ez]['Subunits'][p]
                ProteinMatrix[row_ind, col_ind] += coeff
    return({'Matrix': numpy.array(ProteinMatrix), 'Consumers': Consumers, 'Proteins': Proteins})


def _import_uniprot_file(xml_dir):
    if os.path.isfile(str(xml_dir+'/uniprot.csv')):
        return(pandas.read_csv(str(xml_dir+'/uniprot.csv'), sep='\t'))
    elif os.path.isfile(str(xml_dir+'/data/uniprot.csv')):
        return(pandas.read_csv(str(xml_dir+'/data/uniprot.csv'), sep='\t'))
    else:
        print('WARNING: Uniprot-file "uniprot.csv" not found.\n' + ' Continuing without additional information...\n')
        return(str('Not There'))


def _import_sbml_file(xml_dir, filename):
    if os.path.isfile(str(xml_dir+'/'+filename)):
        SBfile = libsbml.readSBML(str(xml_dir+'/'+filename))
        if SBfile.getNumErrors() > 0:
            SBfile.printErrors()
            print('WARNING: Invalid SBML')
            return(str('Not There'))
        else:
            sbml = SBfile
            return(sbml)
    elif os.path.isfile(str(xml_dir+'/data/'+filename)):
        SBfile = libsbml.readSBML(str(xml_dir+'/data/'+filename))
        if SBfile.getNumErrors() > 0:
            SBfile.printErrors()
            print('WARNING: Invalid SBML')
            return(str('Not There'))
        else:
            sbml = SBfile
            return(sbml)
    else:
        print('WARNING: SBML-file {} not found.\n' + ' Continuing without additional information...\n'.format(filename))
        return(str('Not There'))


def _import_gene_annotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/GeneAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/GeneAnnotations.csv'), sep=',', index_col=0)
        if not list(out):
            print('WARNING: File "GeneAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(str('Not There'))
        return(out)
    elif os.path.isfile(str(xml_dir+'/data/GeneAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/data/GeneAnnotations.csv'), sep=',', index_col=0)
        if not list(out):
            print('WARNING: File "GeneAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(str('Not There'))
        return(out)
    else:
        print('WARNING: No Gene-annotation file "GeneAnnotations.csv" provided.\n' + ' Continuing without additional information...\n')
        return(str('Not There'))


def _import_reaction_annotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/ReactionAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/ReactionAnnotations.csv'), sep=',', index_col=0)
        if not list(out):
            print('WARNING: File "ReactionAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(str('Not There'))
        return(out)
    elif os.path.isfile(str(xml_dir+'/data/ReactionAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/data/ReactionAnnotations.csv'), sep=',', index_col=0)
        if not list(out):
            print('WARNING: File "ReactionAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(str('Not There'))
        return(out)
    else:
        print('WARNING: No Reaction-annotation file "ReactionAnnotations.csv" provided.\n' + ' Continuing without additional information...\n')
        return(str('Not There'))


def _import_metabolite_annotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/MetaboliteAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/MetaboliteAnnotations.csv'), sep=',', index_col=0)
        if not list(out):
            print('WARNING: File "MetaboliteAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(str('Not There'))
        return(out)
    elif os.path.isfile(str(xml_dir+'/data/MetaboliteAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/data/MetaboliteAnnotations.csv'), sep=',', index_col=0)
        if not list(out):
            print('WARNING: File "MetaboliteAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(str('Not There'))
        return(out)
    else:
        print('WARNING: No Reaction-annotation file "MetaboliteAnnotations.csv" provided.\n' + ' Continuing without additional information...\n')
        return(str('Not There'))


def _import_model_info(xml_dir):
    if os.path.isfile(str(xml_dir+'/ModelInformation.csv')):
        out = pandas.read_csv(str(xml_dir+'/ModelInformation.csv'),sep=',', header=0)
        out.index = list(out['Key'])
        if list(out):
            print('WARNING: File "ModelInformation.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(pandas.DataFrame([['Name', 'ModelName'], ['Author', 'John Doe'], ['Organism', 'Life'], ['Reconstruction', 'GSMM'], ['SBML-file', 'Not Provided']], index=['Name', 'Author', 'Organism', 'Reconstruction', 'SBML-file'], columns=['Key', 'Value']))
        return(out)
    elif os.path.isfile(str(xml_dir+'/data/ModelInformation.csv')):
        out = pandas.read_csv(str(xml_dir+'/data/ModelInformation.csv'),sep=',', header=0)
        out.index = list(out['Key'])
        if list(out):
            print('WARNING: File "ModelInformation.csv" seems to be empty or has the wrong delimiter (comma required).')
            return(pandas.DataFrame([['Name', 'ModelName'], ['Author', 'John Doe'], ['Organism', 'Life'], ['Reconstruction', 'GSMM'], ['SBML-file', 'Not Provided']], index=['Name', 'Author', 'Organism', 'Reconstruction', 'SBML-file'], columns=['Key', 'Value']))
        return(out)
    else:
        print('WARNING: No model-info file "ModelInformation.csv" provided.\n' + ' Using dummy-information\n')
        return(pandas.DataFrame([['Name', 'ModelName'], ['Author', 'John Doe'], ['Organism', 'Life'], ['Reconstruction', 'GSMM'], ['SBML-file', 'Not Provided']], index=['Name', 'Author', 'Organism', 'Reconstruction', 'SBML-file'], columns=['Key', 'Value']))


def _generate_overview(StatsReactions, StatsMetabolites, StatsModules, StatsEnzymes, StatsProteins, StatsMacromolecules, StatsProcesses, StatsCompartments, StatsTargets, StatsConstraintsBiological, StatsConstraintsMathematical):
    out = {'Reactions Total': StatsReactions['ReactionsTotal'],
           'Reactions Unique': StatsReactions['ReactionsUnique'],
           'Reactions Spontaneous': StatsReactions['ReactionsSpontaneous'],
           'Reactions Enzymatic': StatsReactions['ReactionsEnzymatic'],
           'Reactions Internal': StatsReactions['ReactionsInternal'],
           'Reactions Exchange': StatsReactions['ReactionsExchange'],
           'Reactions CompartmentTransport': StatsReactions['ReactionsCompartmentTransport'],
           'Reactions Reversible': StatsReactions['ReactionsReversible'],
           'Reactions Irreversible': StatsReactions['ReactionsIrreversible'],
           'Metabolites Total': StatsMetabolites['MetabolitesTotal'],
           'Metabolites Internal': StatsMetabolites['MetabolitesInternal'],
           'Metabolites External': StatsMetabolites['MetabolitesExternal'],
           'Metabolites GrowthRelevant': StatsMetabolites['MetabolitesGrowthRelevant'],
           'Boundary Metabolites': StatsMetabolites['BoundaryMetabolites'],
           'Enzymes Total': StatsEnzymes['EnzymesTotal'],
           'Enzymes Unique': StatsEnzymes['EnzymesUnique'],
           'Proteins Total': StatsProteins['ProteinsTotal'],
           'RNAs Total': StatsMacromolecules['RNAsTotal'],
           'DNAs Total': StatsMacromolecules['DNAsTotal'],
           'Processes Total': StatsProcesses['ProcessesTotal'],
           'Modules Total': StatsModules['ModulesTotal'],
           'Compartments Total': StatsCompartments['CompartmentsTotal'],
           'Biological constraints metabolite': StatsConstraintsBiological['BioConstraintsMetabolite'],
           'Biological constraints capacity': StatsConstraintsBiological['BioConstraintsCapacity'],
           'Biological constraints process': StatsConstraintsBiological['BioConstraintsProcess'],
           'Biological constraints density': StatsConstraintsBiological['BioConstraintsDensity'],
           'Mathematical constraints variables': StatsConstraintsMathematical['MathConstraintsVariables'],
           'Mathematical constraints constraints': StatsConstraintsMathematical['MathConstraintsConstraints'],
           'Mathematical constraints equalities': StatsConstraintsMathematical['MathConstraintsEqualities'],
           'Mathematical constraints inequalities': StatsConstraintsMathematical['MathConstraintsInequalities']}
    out.update(StatsTargets)
    return(out)


def _stats_constraints_biological(MetCs, CapCs, DenCs, EffCs):
    out = {'BioConstraintsMetabolite': len(MetCs.keys()),
           'BioConstraintsCapacity': len(CapCs.keys()),
           'BioConstraintsProcess': len(EffCs.keys()),
           'BioConstraintsDensity': len(DenCs.keys())}
    return(out)


def _stats_constraints_mathematical(MetCs, CapCs, DenCs, EffCs, Enzymes, Reactions, Processes):
    nVars = len(Enzymes.keys())+len(Reactions.keys())+len(Processes.keys())
    nConsts = len(MetCs.keys())+len(CapCs.keys())+len(DenCs.keys())+len(EffCs.keys())
    nEqC = 0
    nInC = 0
    for i in MetCs.keys():
        if MetCs[i]['Type'] == '<=':
            nInC += 1
        if MetCs[i]['Type'] == '=':
            nEqC += 1
    for i in CapCs.keys():
        if CapCs[i]['Type'] == '<=':
            nInC += 1
        if CapCs[i]['Type'] == '=':
            nEqC += 1
    for i in DenCs.keys():
        if DenCs[i]['Type'] == '<=':
            nInC += 1
        if DenCs[i]['Type'] == '=':
            nEqC += 1
    for i in EffCs.keys():
        if EffCs[i]['Type'] == '<=':
            nInC += 1
        if EffCs[i]['Type'] == '=':
            nEqC += 1
    out = {'MathConstraintsVariables': nVars,
           'MathConstraintsConstraints': nConsts,
           'MathConstraintsEqualities': nEqC,
           'MathConstraintsInequalities': nInC}
    return(out)


def _find_contained_proteins(comp, Prots):
    out = []
    for k in Prots.keys():
        if Prots[k]['Compartment'] == comp:
            out.append(k)
    return(out)


def _find_contained_macromolecules(comp, Macromolecules):
    out = []
    for k in Macromolecules.keys():
        if Macromolecules[k]['Compartment'] == comp:
            out.append(k)
    return(out)


def _find_contained_enzymes(ContainedProteins, Prots, Enzymes):
    rs = []
    enzs = []
    for k in ContainedProteins:
        rs = rs+Prots[k]['associatedReactions']
        enzs = enzs+Prots[k]['associatedEnzymes']
    out = {'containedEnzymes': list(numpy.unique(enzs)),
           'containedReactions': list(numpy.unique(rs))}
    return(out)


def _find_associated_enzyme(Enzymes, protein):
    out1 = []
    out2 = []
    for e in Enzymes.keys():
        if protein in Enzymes[e]['Subunits'].keys():
            out1.append(e)
            out2.append(Enzymes[e]['Reaction'])
    out = {'Enz': out1,
           'Rx': out2}
    return(out)


def _find_iso_enzymes(ez, blocks, Reactions, rx):
    out = []
    twins = Reactions[rx]['Twins']
    if len(twins) > 0:
        for r in twins:
            if not type(r) == list:
                if not type(Reactions[r]['Enzyme']) == list:
                    out.append(Reactions[r]['Enzyme'])
    return(out)


def _sort_constraints(matrix, model):
    metaboliteList = [model.metabolism.species._elements[m].id for m in range(
        len(model.metabolism.species._elements))]
    mets = {}
    capacity = {}
    efficiency = {}
    density = {}
    for j in range(len(matrix.row_names)):
        i = matrix.row_names[j]
        if i in metaboliteList:
            mets[i] = j
        if 'enzyme' in i:
            capacity[i] = j
        if i.startswith('P_'):
            efficiency[i] = j
        if '_density' in i:
            density[i] = j
    out = {'MetaboliteConsts': mets,
           'ProcessConsts': efficiency,
           'EnzymeConsts': capacity,
           'DensityConsts': density}
    return(out)


def _json_int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
