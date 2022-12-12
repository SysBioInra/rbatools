# python 2/3 compatibility
from __future__ import division, print_function

import os.path
import warnings
import numpy
import pandas
import copy
import difflib
import scipy
import rba
import rbatools._auxiliary_functions as _auxiliary_functions
from rba.core import metabolism as rba_core_metabolism
from typing import Union
from rbatools.rba_simulation_data import SimulationDataRBA
from rbatools.rba_simulation_parameters import SimulationParametersRBA
from rbatools.rba_model_structure import ModelStructureRBA
from rbatools.rba_problem import ProblemRBA
from rbatools.fba_problem import ProblemFBA
from rbatools.rba_lp import LinearProblem
from rbatools.rba_problem_matrix import ProblemMatrix
from rbatools._warnings_and_errors import *

warnings.simplefilter('ignore', FutureWarning)

class SessionRBA(object):
    """
    User interface with high level functions to import model, change model,
    different solving methods and results export.

    Attributes
    ----------
    xml_dir : str
        Directory of imported xml model files
    model : rba.RbaModel
        RBA model (parsed from xml files), from which matrices are built
    Problem : rbatools.rba_problem.ProblemRBA
        RBA Problem
    Medium : dict
        Dictionary with external metabolites and corresponding concentrations
    ModelStructure : rbatools.rba_model_structure.ModelStructureRBA
        Model structure representation
    Results : dict
        Simulation results, added if record_results method has been called
    Parameters : dict
        Simulation parameters, added if record_parameters method has been called
    SimulationData : rbatools.rba_simulation_data.SimulationDataRBA
        SimulationData object, added if write_results method has been called
    SimulationParameters : rbatools.rba_simulation_parameters.SimulationParametersRBA
        SimulationParameters object, added if write_results method has been called
    ExchangeMap :  dict
        Map of all metabolites, the corresponding transport-reactions and stoichiometires
    ExchangeReactions : bool --> default:False
        True if exchange reactions have been added via add_exchange_reactions method
    """

    def __init__(self, xml_dir: str,lp_solver="cplex"):
        """
        Creates SessionRBA object from files

        Parameters
        ----------
        xml_dir : str
            Path to the directory where rba-model files are located.
        lp_solver : str
            Specifies which LP-solver should be used ('cplex' or 'swiglpk')
            Default: 'cplex'
        """
        self.xml_dir = xml_dir
        self.lp_solver=lp_solver

        if not hasattr(self, 'ModelStructure'):
            if os.path.isfile(str(self.xml_dir+'/ModelStructure.json')):
                self.ModelStructure = ModelStructureRBA()
                with open(str(self.xml_dir+'/ModelStructure.json'), 'r') as myfile:
                    data = myfile.read()
                self.ModelStructure.from_json(inputString=data)
            else:
                self.build_model_structure()
        self.model = rba.RbaModel.from_xml(input_dir=xml_dir)
        self.Problem = ProblemRBA(matrix=rba.ConstraintMatrix(model=self.model),ModelStructure=self.ModelStructure,lp_solver=self.lp_solver)

        medium = pandas.read_csv(xml_dir+'/medium.tsv', sep='\t')
        self.Medium = dict(zip(list(medium.iloc[:, 0]), [float(i) for i in list(medium.iloc[:, 1])]))

        self.Mu = self.Problem.Mu
        self.ExchangeMap = _auxiliary_functions.build_exchange_map(RBA_Session=self)
        self.ExchangeReactions=False

    def build_model_structure(self):
        """
        Builds model structure object from model xml-files, adds it as
        ModelStructure attribute and stores as json file.
        """
        self.ModelStructure = ModelStructureRBA()
        self.ModelStructure.from_files(xml_dir=self.xml_dir)
        self.ModelStructure.export_json(path=self.xml_dir)
        with open(str(self.xml_dir+'/ModelStructure.json'), 'r') as myfile:
            data = myfile.read()
        self.ModelStructure.from_json(inputString=data)

    def add_exchange_reactions(self):
        """
        Adds explicit exchange-reactions of boundary-metabolites to RBA-problem,
        named R_EX_ followed by metabolite name (without M_ prefix).
        """
        Mets_external = [m.id for m in self.model.metabolism.species if m.boundary_condition]
        Mets_internal = [m.id for m in self.model.metabolism.species if not m.boundary_condition]
        Reactions = [r.id for r in self.model.metabolism.reactions]
        full_S = rba_core_metabolism.build_S(Mets_external+Mets_internal, self.model.metabolism.reactions)
        S_M_ext = full_S[:len(Mets_external), ].toarray()
        col_indices_toremove = []
        for i in range(S_M_ext.shape[1]):
            s_col_uniques = list(set(list(S_M_ext[:, i])))
            if len(s_col_uniques) == 1:
                if s_col_uniques[0] == 0:
                    col_indices_toremove.append(i)
        RemainingReactions = [i for i in Reactions if Reactions.index(
            i) not in col_indices_toremove]
        S_ext = numpy.delete(S_M_ext, col_indices_toremove, axis=1)
        A = numpy.concatenate((S_ext, numpy.eye(len(Mets_external))), axis=1, out=None)
        ColNames = RemainingReactions+[str('R_EX_'+i.split('M_')[-1]) for i in Mets_external]
        # print(str('R_EX_'+i.split('M_')[-1]))
        LBs = list([self.Problem.LP.LB[self.Problem.LP.col_names.index(i)]
                    for i in RemainingReactions]+[-10000]*len(Mets_external))
        UBs = list([self.Problem.LP.UB[self.Problem.LP.col_names.index(i)]
                    for i in RemainingReactions]+[10000]*len(Mets_external))
        b = [0]*len(Mets_external)
        f = list([self.Problem.LP.f[self.Problem.LP.col_names.index(i)]
                  for i in RemainingReactions]+[0]*len(Mets_external))

        ExchangeMatrix = ProblemMatrix()
        ExchangeMatrix.A = scipy.sparse.coo_matrix(A)
        ExchangeMatrix.b = numpy.array([0]*len(Mets_external))
        ExchangeMatrix.f = numpy.array(f)
        ExchangeMatrix.LB = numpy.array(LBs)
        ExchangeMatrix.UB = numpy.array(UBs)
        ExchangeMatrix.row_signs = ['E']*len(Mets_external)
        ExchangeMatrix.row_names = Mets_external
        ExchangeMatrix.col_names = ColNames
        ExchangeMatrix.map_indices()
        self.Problem.LP.add_matrix(matrix=ExchangeMatrix)

        self.ExchangeReactionMap = dict(zip([_auxiliary_functions.find_exchange_metabolite_in_medium(metabolite=i,Medium=self.Medium) for i in Mets_external], [str('R_EX_'+i.split('M_')[-1]) for i in Mets_external]))
        self.ExchangeReactions=True

    def rebuild_from_model(self):
        """
        Rebuilds computational model-representation from own attribute "model" (rba.RbaModel-object).
        """
        self.Problem = ProblemRBA(matrix=rba.ConstraintMatrix(model=self.model),ModelStructure=self.ModelStructure,lp_solver=self.lp_solver)
        self.set_medium(changes=self.Medium)
        if self.ExchangeReactions:
            self.add_exchange_reactions()

    def reload_model(self):
        """
        Reloads model from xml-files and then rebuild computational model-representation.
        """
        self.model = rba.RbaModel.from_xml(input_dir=self.xml_dir)
        self.rebuild_from_model()

    def record_results(self, run_name: str):
        """
        Records Simulation output for further use.
        and stores them in own 'Results'-attribute as pandas.DataFrames in a dictionary
        with the respective run-name being a column in all DataFrames.

        Parameters
        ----------
        run_name : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        if not hasattr(self, 'Results'):
            self.Results = {'Reactions': pandas.DataFrame(index=list(self.ModelStructure.ReactionInfo.Elements.keys())),
                            'Enzymes': pandas.DataFrame(index=list(self.ModelStructure.EnzymeInfo.Elements.keys())),
                            'Processes': pandas.DataFrame(index=[self.ModelStructure.ProcessInfo.Elements[i]['ID']+'_machinery' for i in self.ModelStructure.ProcessInfo.Elements.keys()]),
                            'Proteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinMatrix['Proteins'])),
                            'ProtoProteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinGeneMatrix['ProtoProteins'])),
                            'Constraints': pandas.DataFrame(index=self.Problem.LP.row_names),
                            'SolutionType': pandas.DataFrame(index=['SolutionType']),
                            'ObjectiveFunction': pandas.DataFrame(index=self.Problem.LP.col_names),
                            'Mu': pandas.DataFrame(index=['Mu']),
                            'ObjectiveValue': pandas.DataFrame(index=['ObjectiveValue']),
                            'ExchangeFluxes': pandas.DataFrame(index=list(self.ExchangeMap.keys()))}

        Exchanges = self.return_exchange_fluxes()
        for i in Exchanges.keys():
            self.Results['ExchangeFluxes'].loc[i, run_name] = _auxiliary_functions.check_solution_feasibility(
                Value=Exchanges[i], RBA_Session=self)

        self.Results['Reactions'][run_name] = [_auxiliary_functions.check_solution_feasibility(Value=self.Problem.SolutionValues[i], RBA_Session=self) for i in list(self.Results['Reactions'].index)]
        self.Results['Enzymes'][run_name] = [_auxiliary_functions.check_solution_feasibility(
            Value=self.Problem.SolutionValues[i], RBA_Session=self) for i in list(self.Results['Enzymes'].index)]
        self.Results['Processes'][run_name] = [_auxiliary_functions.check_solution_feasibility(
            Value=self.Problem.SolutionValues[i], RBA_Session=self) for i in list(self.Results['Processes'].index)]
        self.Results['Constraints'][run_name] = [_auxiliary_functions.check_solution_feasibility(
            Value=self.Problem.DualValues[i], RBA_Session=self) for i in self.Problem.LP.row_names]
        self.Results['Proteins'][run_name] = _auxiliary_functions.record_proteome(RBA_Session=self, run=run_name)
        self.Results['ProtoProteins'][run_name] = _auxiliary_functions.record_proto_proteome(RBA_Session=self, run=run_name, Proteinlevels=self.Results['Proteins'])
        self.Results['SolutionType'][run_name] = _auxiliary_functions.check_solution_feasibility(Value=self.Problem.SolutionType, RBA_Session=self)
        self.Results['Mu'][run_name] = _auxiliary_functions.check_solution_feasibility(Value=self.Problem.Mu, RBA_Session=self)
        self.Results['ObjectiveValue'][run_name] = _auxiliary_functions.check_solution_feasibility(
            Value=self.Problem.ObjectiveValue, RBA_Session=self)
        self.Results['ObjectiveFunction'][run_name] = [0.0]*self.Results['ObjectiveFunction'].shape[0]
        objective_function=self.Problem.get_objective()
        for i in objective_function.keys():
            self.Results['ObjectiveFunction'].loc[i,run_name]=objective_function[i]

    def record_parameters(self, run_name: str):
        """
        Records Simulation parameters (LP-coefficients etc.) for further use.
        and stores them in own 'Parameters'-attribute as pandas.DataFrames in a dictionary
        with the respective run-name being a column in all DataFrames.

        Parameters
        ----------
        run_name : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        enzymes=self.get_enzymes()
        processes=[i for i in self.get_processes() if self.get_process_information(i)["Capacity_Constraint"] in self.get_process_constraints()]
        compartments=[i for i in self.get_compartments() if self.get_compartment_information(i)["Capacity_Constraint"] in self.get_density_constraints()]

        if not hasattr(self, 'Parameters'):
            self.Parameters = {'EnzymeEfficiencies_FW': pandas.DataFrame(index=enzymes),
                               'EnzymeEfficiencies_BW': pandas.DataFrame(index=enzymes),
                               'NetProcessEfficiencies': pandas.DataFrame(index=processes),
                               'CompartmentCapacities': pandas.DataFrame(index=compartments)}
        self.Parameters['EnzymeEfficiencies_FW'][run_name] = [self.get_current_parameter_value(self.get_enzyme_constraint_information(self.get_enzyme_information(enzyme)["ForwardCapacity_Constraint"])['CapacityParameterID']) for enzyme in self.Parameters['EnzymeEfficiencies_FW'].index]
        self.Parameters['EnzymeEfficiencies_BW'][run_name] = [self.get_current_parameter_value(self.get_enzyme_constraint_information(self.get_enzyme_information(enzyme)["BackwardCapacity_Constraint"])['CapacityParameterID']) for enzyme in self.Parameters['EnzymeEfficiencies_BW'].index]
        self.Parameters['NetProcessEfficiencies'][run_name] =[self.get_current_parameter_value(self.get_process_constraint_information(self.get_process_information(process)["Capacity_Constraint"])['CapacityParameterID']) for process in self.Parameters['NetProcessEfficiencies'].index]
        self.Parameters['CompartmentCapacities'][run_name] = [self.get_current_parameter_value(self.get_density_constraint_information(self.get_compartment_information(compartment)['Capacity_Constraint'])["CapacityParameterID"]) for compartment in self.Parameters['CompartmentCapacities'].index]

    def clear_results_and_parameters(self):
        """
        Removes all previosly recorded results and deletes own 'Results'-attribute.
        """
        if hasattr(self, 'Results'):
            delattr(self, 'Results')
        if hasattr(self, 'Parameters'):
            delattr(self, 'Parameters')

    def write_results(self, session_name: str, digits: int = 10):
        """
        Creates SimulationData and SimulationParameters objects from recordings ('Results'.'Parameters').

        Stores them as rbatools.rba_simulation_data.SimulationDataRBA
        and rbatools.rba_simulation_parameters.SimulationParametersRBA objects as attributes.
        Access via attributes .SimulationData and SimulationParameters respectively.

        Parameters
        ----------
        digits : int
            Number of decimal places in the numeric results
            Default: 10
        session_name : str
            Name of Simulation session.
        """
        if hasattr(self, 'Results'):
            self.Results['uniqueReactions'] = _auxiliary_functions.map_iso_reactions(RBA_Session=self)
            self.Results['Mu'] = self.Results['Mu'].round(digits)
            self.Results['ObjectiveValue'] = self.Results['ObjectiveValue'].round(digits)
            self.Results['Proteins'] = self.Results['Proteins'].round(digits)
            self.Results['uniqueReactions'] = self.Results['uniqueReactions'].round(digits)
            self.Results['Reactions'] = self.Results['Reactions'].round(digits)
            self.Results['Enzymes'] = self.Results['Enzymes'].round(digits)
            self.Results['Processes'] = self.Results['Processes'].round(digits)
            self.Results['Constraints'] = self.Results['Constraints'].round(digits)
            self.Results['ExchangeFluxes'] = self.Results['ExchangeFluxes'].round(digits)

            self.SimulationData = SimulationDataRBA(rbaModelStructure=self.ModelStructure)
            self.SimulationData.from_simulation_results(rbaSession=self, session_name=session_name)

        if hasattr(self, 'Parameters'):
            self.Parameters['EnzymeEfficiencies_FW'] = self.Parameters['EnzymeEfficiencies_FW'].round(
                digits)
            self.Parameters['EnzymeEfficiencies_BW'] = self.Parameters['EnzymeEfficiencies_BW'].round(
                digits)
            self.Parameters['NetProcessEfficiencies'] = self.Parameters['NetProcessEfficiencies'].round(
                digits)
            self.Parameters['CompartmentCapacities'] = self.Parameters['CompartmentCapacities'].round(
                digits)
            self.SimulationParameters = SimulationParametersRBA()
            self.SimulationParameters.from_simulation_results(rbaSession=self)

    def return_exchange_fluxes(self) -> dict:
        """
        Generates a dictonary with the exchang-rates of boundary-metabolites.

        Returns
        -------
        Dictonary with exchange-keys and respective -rates.
        """
        out = {}
        for j in self.ExchangeMap.keys():
            netflux = 0
            for k in self.ExchangeMap[j].keys():
                netflux += self.ExchangeMap[j][k]*self.Problem.SolutionValues[k]
            out[j] = netflux
        return(out)

    def set_growth_rate(self, Mu: float):
        """
        Sets growth-rate to desired value.

        Parameters
        ----------
        Mu : float
            Growth rate
        """
        self.Problem.set_growth_rate(Mu=float(Mu), ModelStructure=self.ModelStructure)
        self.Mu = float(Mu)

    def solve(self, run_name: str = 'DontSave', feasible_stati: list = ["optimal","feasible"], try_unscaling_if_sol_status_is_feasible_only_before_unscaling: bool = True) -> bool:
        """
        Solves problem to find a solution.

        Parameters
        ----------
        run_name : str
            Name of observation.
            Serves as ID for all data, originating from this run.
            Special values :
                'DontSave' : Results are not recorded
                'Auto' : Results are automatically recorded
                         and appended to existing ones.
                    Named with number.
                Any other string: Results are recorded under this name.
            Default: 'DontSave'
        feasible_stati : list of str
            List with acceptable solution statuses.
            Default: ["optimal","feasible"]
        try_unscaling_if_sol_status_is_feasible_only_before_unscaling : bool
            If true; the problem will be attempted to be solved without scaling,
            if the scaled problem is feasible but the solution is not feasible
            after unscaling (CPLEX solution-status 5).
            Default: True
        Returns
        -------
        Bool, indicating if problem could be solved (feasible) or not (infeasible).
        """

        self.Problem.solve_lp(feasible_stati=feasible_stati, try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        if self.Problem.Solved:
            if run_name is not 'DontSave':
                if run_name is 'Auto':
                    if hasattr(self, 'Results'):
                        name = str(self.Results['Reactions'].shape[1]+1)
                    if not hasattr(self, 'Results'):
                        name = '1'
                if run_name is not 'Auto':
                    name = run_name
                self.record_results(run_name=name)
            return(True)
        else:
            return(False)

    def rbapy_mumax_search(self):
        solver=rba.Solver(matrix=self.Problem.ClassicRBAmatrix)
        solver.solve()
        return(solver.mu_opt)

    def find_max_growth_rate(self, precision: float = 0.001, max: float = 4.0, start_value: float = numpy.nan, recording: bool = False, omit_objective: bool = False, feasible_stati: list = ["optimal","feasible"], try_unscaling_if_sol_status_is_feasible_only_before_unscaling: bool = True) -> float:
        """
        Applies dichotomy-search to find the maximal feasible growth-rate.

        Parameters
        ----------
        precision : float
            Numberic precision with which maximum is approximated.
            Default: 0.00001
        max : float
            Defines the highest growth rate to be screened for.
            Default: 4.0
        start_value : float
            Defines the first growth-rate to test during the dichotomy search.
            Default: numpy.nan --> then the middle between 0 and max is used.
        recording : bool
            Records intermediate feasible solutions
            while approaching the maximum growth-rate.
            Default: False
        feasible_stati : list of str
            List with acceptable solution statuses.
            Default: ["optimal","feasible"]
        try_unscaling_if_sol_status_is_feasible_only_before_unscaling : bool
            If true; the problem will be attempted to be solved without scaling,
            if the scaled problem is feasible but the solution is not feasible
            after unscaling (CPLEX solution-status 5).
            Default: True

        Returns
        -------
        maximum feasible growth rate as float.
        """

        minMu = 0
        maxMu = max
        if numpy.isfinite(start_value):
            testMu = start_value
        else:
            testMu =0
            #testMu = max/2
        iteration = 0

        if omit_objective:
            old_Obj = self.Problem.get_objective()
            self.Problem.clear_objective()
        optMu=testMu
        while (maxMu - minMu) > precision:
            self.set_growth_rate(Mu=testMu)
            self.Problem.solve_lp(feasible_stati=feasible_stati,try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
            if self.Problem.Solved:
                iteration += 1
                if recording:
                    self.record_results('DichotomyMu_iteration_'+str(iteration))
                minMu = testMu
                optMu = testMu
            else:
                maxMu = testMu
            testMu = (maxMu+minMu)/2
        effective_Mu_difference=maxMu - minMu
        if optMu == max:
            print('WARNING: Maximum growth rate might exceed specified range. Try rerunning this method with larger "max" argument.')
        if omit_objective:
            self.Problem.set_objective(old_Obj)
        self.set_growth_rate(Mu=optMu)
        self.Problem.solve_lp(feasible_stati=feasible_stati,try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
        return(optMu)

    def find_min_substrate_concentration(self, metabolite: str, precision: float =0.00001, max: float =100.0, recording:bool =False, feasible_stati: list = ["optimal","feasible"], try_unscaling_if_sol_status_is_feasible_only_before_unscaling: bool = True):
        """
        Applies dichotomy-search to find the minimal feasible concentration of
        growth-substrate in medium, at a previously set growth-rate.

        Parameters
        ----------
        metabolite : str
            ID of metabolite in medium.
        precision : float
            Numberic precision with which minimum is approximated.
            Default: 0.00001
        max : float
            Defines the highest concentration rate to be screened for.
            Default: 100
        recording : bool
            Records intermediate feasible solutions
            while approaching the minimum concentration.
            Default: False
        feasible_stati : list of str
            List with acceptable solution statuses.
            Default: ["optimal","feasible"]
        try_unscaling_if_sol_status_is_feasible_only_before_unscaling : bool
            If true; the problem will be attempted to be solved without scaling,
            if the scaled problem is feasible but the solution is not feasible
            after unscaling (CPLEX solution-status 5).
            Default: True

        Returns
        -------
        minimum feasible growth-substrate concentration as float.
        """

        minConc = 0.0
        maxConc = max
        testConc = minConc
        iteration = 0
        oldConc = self.Medium[metabolite]
        while (maxConc - minConc) > precision:
            self.set_medium(changes={metabolite: testConc})
            self.Problem.solve_lp(feasible_stati=feasible_stati,try_unscaling_if_sol_status_is_feasible_only_before_unscaling=try_unscaling_if_sol_status_is_feasible_only_before_unscaling)
            if self.Problem.Solved:
                iteration += 1
                if recording:
                    run_name = 'Dichotomy_'+metabolite+'_' + \
                        str(testConc)+'_iteration_'+str(iteration)
                    self.record_results(run_name)
                maxConc = testConc
            else:
                minConc = testConc
            testConc = numpy.mean([maxConc, minConc])
        self.set_medium(changes={metabolite: oldConc})
        return(maxConc)

    def get_medium(self) -> dict:
        """
        Returns dictionary of current medium-composition, 
        with external metabolites as keys and their concentrations as values.

        Returns
        ----------
        Dict: medium composition
        """
        return(self.Medium)

    def set_medium(self, changes: dict):
        """
        Sets the concentration of specified growth-substrate(s) in medium.

        Parameters
        ----------
        changes : dict
            Keys : ID of metabolite(s) in medium.
            Values : New concention(s)
        """

        for species in (changes.keys()):
            if species in self.Medium.keys():
                self.Medium[species] = float(changes[species])
            else:
                raise InputError('{} not a medium component'.format(species))


        self.Problem.ClassicRBAmatrix.set_medium(self.Medium)
        self.set_growth_rate(Mu=self.Mu)

    def undo_gene_knock_out(self, gene: Union[list,str]):
        """
        Undoes a gene knock out.
            Removes all constraints, on knocking out the provided genes

        Parameters
        ----------
        gene : str or list of strings
            ID(s) of model-proteins to be un-knocked out.
            Can either be gene-identifier, represented as ID or ProtoID of
            proteins in rbatools.protein_bloc.ProteinBlock.Elements class
            (depends on whether protein-isoforms are considered).
        """
        if type(gene) is str:
            genes = [gene]
        if type(gene) is list:
            genes = gene
        isoform_genes = [g for g in genes if g in list(self.ModelStructure.ProteinInfo.Elements.keys())]+[i for g in genes for i in self.ModelStructure.ProteinInfo.Elements.keys() if self.ModelStructure.ProteinInfo.Elements[i]['ProtoID'] == g]
        for g in isoform_genes:
            ConsumersEnzymes = self.ModelStructure.ProteinInfo.Elements[g]['associatedEnzymes']
            for i in ConsumersEnzymes:
                LikeliestVarName = difflib.get_close_matches(i, self.Problem.LP.col_names, 1)[0]
                self.Problem.MuDependencies['FromMatrix']['LB'].append(LikeliestVarName)
                self.Problem.MuDependencies['FromMatrix']['UB'].append(LikeliestVarName)
            ConsumersProcess = self.ModelStructure.ProteinInfo.Elements[g]['SupportsProcess']
            for i in ConsumersProcess:
                LikeliestVarName = difflib.get_close_matches(str(self.ModelStructure.ProcessInfo.Elements[i]['ID']+'_machinery'), self.Problem.LP.col_names, 1)[0]
                self.Problem.MuDependencies['FromMatrix']['LB'].append(LikeliestVarName)
                self.Problem.MuDependencies['FromMatrix']['UB'].append(LikeliestVarName)

    def apply_gene_knock_out(self, gene: Union[list,str]):
        """
        Simulates a gene knock out.
        Constrains all variables in the LP-problem (enzymes, other machineries),
        which require this gene(s), to zero.

        Parameters
        ----------
        gene : str or list of strings
            ID(s) of model-proteins to be knocked out.
            Can either be gene-identifier, represented as ID or ProtoID of
            proteins in rbatools.protein_bloc.ProteinBlock.Elements class
            (depends on whether protein-isoforms are considered).
        """

        if type(gene) is str:
            genes = [gene]
        if type(gene) is list:
            genes = gene
        isoform_genes = [g for g in genes if g in list(self.ModelStructure.ProteinInfo.Elements.keys())]+[i for g in genes for i in self.ModelStructure.ProteinInfo.Elements.keys() if self.ModelStructure.ProteinInfo.Elements[i]['ProtoID'] == g]
        for g in isoform_genes:
            ConsumersEnzymes = self.ModelStructure.ProteinInfo.Elements[g]['associatedEnzymes']
            for i in ConsumersEnzymes:
                LikeliestVarName = difflib.get_close_matches(i, self.Problem.LP.col_names, 1)[0]
                self.Problem.set_lb(inputDict={LikeliestVarName: 0.0})
                self.Problem.set_ub(inputDict={LikeliestVarName: 0.0})
                if LikeliestVarName in self.Problem.MuDependencies['FromMatrix']['LB']:
                    self.Problem.MuDependencies['FromMatrix']['LB'].remove(LikeliestVarName)
                if LikeliestVarName in self.Problem.MuDependencies['FromMatrix']['UB']:
                    self.Problem.MuDependencies['FromMatrix']['UB'].remove(LikeliestVarName)
            ConsumersProcess = self.ModelStructure.ProteinInfo.Elements[g]['SupportsProcess']
            for i in ConsumersProcess:
                LikeliestVarName = difflib.get_close_matches(str(self.ModelStructure.ProcessInfo.Elements[i]['ID']+'_machinery'), self.Problem.LP.col_names, 1)[0]
                self.Problem.set_lb(inputDict={LikeliestVarName: 0.0})
                self.Problem.set_ub(inputDict={LikeliestVarName: 0.0})
                if LikeliestVarName in self.Problem.MuDependencies['FromMatrix']['LB']:
                    self.Problem.MuDependencies['FromMatrix']['LB'].remove(LikeliestVarName)
                if LikeliestVarName in self.Problem.MuDependencies['FromMatrix']['UB']:
                    self.Problem.MuDependencies['FromMatrix']['UB'].remove(LikeliestVarName)

    def get_feasible_range(self, variables: Union[list,str,dict]=[]) -> dict:
        """
        Determines the feasible range of model variables or linear combinations of model variables.

        Parameters
        ----------
        variables : str or list of str or dict
            Specifies variable(s) for which the feasible range is to be determined.
            If str: Feasible range of the respective variable is determined
            If list of str: Feasible range for each variable in list is determined
                            (If empty list, all model variables are chosen)
            If dict: Dictionary with user-specified ID of objective as key and dictionary 
            with variable IDs as values and stoichiometric coefficients in objective as values.
            {'ID_objective':{'model_variable_1':coeff,'model_variable_2':coeff, ...}}
            Default: []
                All model-variables are taken one by one

        Returns
        -------
        Dictionary with variable-names or objective-names as keys and other dictionaries as values.
        The 'inner' dictionaries hold keys 'Min' and 'Max'
        with values representing lower and upper bound of feasible range respectively.
        E.g. : {'variableA':{'Min':42 , 'Max':9000},'variableB':{'Min':-9000 , 'Max':-42}...}
        or {'ID_objective':{'Min':42 , 'Max':9000}}

        """
        out = {}
        if isinstance(variables, dict):
            for objective in variables.keys():
                min_value = numpy.nan
                max_value = numpy.nan
                self.Problem.clear_objective()
                self.Problem.set_objective(inputDict=variables[objective])
                self.Problem.solve_lp(feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False)
                if self.Problem.Solved:
                    min_value = self.Problem.ObjectiveValue
                self.Problem.set_objective(inputDict={i:-variables[objective][i] for i in variables[objective].keys()})
                self.Problem.solve_lp(feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False)
                if self.Problem.Solved:
                    max_value = -self.Problem.ObjectiveValue
                out.update({objective: {'Min': min_value, 'Max': max_value}})
        else:
            if isinstance(variables, list):
                if len(variables) > 0:
                    VariablesInQuestion = variables
                else:
                    VariablesInQuestion = self.Problem.LP.col_names
            elif isinstance(variables, str):
                VariablesInQuestion = [variables]
            for i in VariablesInQuestion:
                min_value = numpy.nan
                max_value = numpy.nan
                self.Problem.clear_objective()
                self.Problem.set_objective(inputDict={i: 1.0})
                self.Problem.solve_lp(feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False)
                if self.Problem.Solved:
                    min_value = self.Problem.ObjectiveValue
                self.Problem.set_objective(inputDict={i: -1.0})
                self.Problem.solve_lp(feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False)
                if self.Problem.Solved:
                    max_value = -self.Problem.ObjectiveValue
                out.update({i: {'Min': min_value, 'Max': max_value}})

        return(out)

    def get_pareto_front(self, variable_X: str, variable_Y: str, N: int = 10, sign_VY: str = 'max') -> pandas.core.frame.DataFrame:
        """
        Determine Pareto front of two model variables or linear combinations of model variables.
        Initially the feasible range of the X variable is determined and the Y variable is optimised 
        at different values of the X variable.

        Parameters
        ----------
        variable_X : str or dictionary
            If str : ID of variable, representing the X-coordinate of the Pareto-front
            If dict: Dictionary with user-specified ID of objective as key and dictionary 
            with variable IDs as values and stoichiometric coefficients in objective as values.
            {'ID_objective_X':{'model_variable_1':coeff,'model_variable_2':coeff, ...}}
            IMPORTANT: ID_objective_X must not equal any model variable ID.
        variable_Y : str
            If str : ID of variable, representing the Y-coordinate of the Pareto-front
            If dict: Dictionary with user-specified ID of objective as key and dictionary 
            with variable IDs as values and stoichiometric coefficients in objective as values.
            {'ID_objective_Y':{'model_variable_1':coeff,'model_variable_2':coeff, ...}}
            IMPORTANT: ID_objective_Y must not equal any model variable ID.
        N : int
            Number of intervals within the feasible range of X variable.
            Default: 10.
        sign_VY : str
            'max': Y variable is maximised
            'min': Y variable is minimised
            Default: 'max'.
        Returns
        -------
        Pandas DataFrame with columns named after the two input variables
        and 'N' rows. Each row represents an interval on the Pareto front.
        Entries on each row are the X and Y coordinate on the Pareto front,
        representing the values of the two variables.
        """
        
        if isinstance(variable_X, dict):
            x_varname=str(list(variable_X.keys())[0])
            variables_in_linear_combi=list(variable_X[x_varname].keys())
            added_matrix = ProblemMatrix()
            added_matrix.A = scipy.sparse.coo_matrix(numpy.matrix([variable_X[x_varname][i] for i in variables_in_linear_combi]+[-1.0]))
            added_matrix.b = numpy.array([0.0])
            added_matrix.f = numpy.array([0.0]*(len(variables_in_linear_combi)+1))
            added_matrix.LB = numpy.array([self.Problem.get_lb(i)[i] for i in variables_in_linear_combi]+[-10000.0])
            added_matrix.UB = numpy.array([self.Problem.get_ub(i)[i] for i in variables_in_linear_combi]+[10000.0])
            added_matrix.row_signs = ['E']
            added_matrix.row_names = ["constraint_{}".format(x_varname)]
            added_matrix.col_names = variables_in_linear_combi+[x_varname]
            added_matrix.map_indices()
            self.Problem.LP.add_matrix(matrix=added_matrix)
        elif isinstance(variable_X, str):
            x_varname=variable_X
        if isinstance(variable_Y, dict):
            y_varname=str(list(variable_Y.keys())[0])
            variables_in_linear_combi=list(variable_Y[y_varname].keys())
            added_matrix = ProblemMatrix()
            added_matrix.A = scipy.sparse.coo_matrix(numpy.matrix([variable_Y[y_varname][i] for i in variables_in_linear_combi]+[-1.0]))
            added_matrix.b = numpy.array([0.0])
            added_matrix.f = numpy.array([0.0]*(len(variables_in_linear_combi)+1))
            added_matrix.LB = numpy.array([self.Problem.get_lb(i)[i] for i in variables_in_linear_combi]+[-10000.0])
            added_matrix.UB = numpy.array([self.Problem.get_ub(i)[i] for i in variables_in_linear_combi]+[10000.0])
            added_matrix.row_signs = ['E']
            added_matrix.row_names = ["constraint_{}".format(y_varname)]
            added_matrix.col_names = variables_in_linear_combi+[y_varname]
            added_matrix.map_indices()
            self.Problem.LP.add_matrix(matrix=added_matrix)
        elif isinstance(variable_Y, str):
            y_varname=variable_Y

        FR_x = self.get_feasible_range(x_varname)
        x_min = FR_x[x_varname]['Min']
        x_max = FR_x[x_varname]['Max']
        x_intervals = [float(x_min+(x_max-x_min)*i/N) for i in range(N+1)]+[x_max]

        Out = pandas.DataFrame(columns=[x_varname, y_varname])
        oldLB = self.Problem.get_lb(x_varname)
        oldUB = self.Problem.get_ub(x_varname)
        iteration = -1

        for x_val in x_intervals:
            iteration += 1
            self.Problem.set_lb(inputDict={x_varname: x_val})
            self.Problem.set_ub(inputDict={x_varname: x_val})
            self.Problem.clear_objective()
            if sign_VY == 'max':
                self.Problem.set_objective(inputDict={y_varname: -1})
            elif sign_VY == 'min':
                self.Problem.set_objective(inputDict={y_varname: 1})
            self.Problem.solve_lp(feasible_stati=["optimal","feasible","feasible_only_before_unscaling"],try_unscaling_if_sol_status_is_feasible_only_before_unscaling=False)
            if self.Problem.Solved:
                if sign_VY == 'max':
                    optimum_y = -self.Problem.ObjectiveValue
                elif sign_VY == 'min':
                    optimum_y = self.Problem.ObjectiveValue
            else:
                optimum_y = numpy.nan
            self.Problem.set_lb(inputDict=oldLB)
            self.Problem.set_ub(inputDict=oldUB)
            Out.loc[iteration, x_varname] = x_val
            Out.loc[iteration, y_varname] = optimum_y
        self.rebuild_from_model()
        return(Out)

    def get_constraint_saturation(self, constraints: Union[list,str]=[]) -> pandas.core.frame.DataFrame:
        """
        Determines the saturation of model constraints at current solution.

        Parameters
        ----------
        constraints : str or list of str
            Specifies constraints(s) for which the saturation is to be determined.
            Optional input:
                If not provided all model-constraints are taken

        Returns
        -------
        Pandas DataFrame with constraint-names as indices and the columns
        'LHS', 'RHS', and 'Saturation'.
            'LHS': The sum over the respoctive constraint-row multiplied
                   elementwise with the solution vector.
            'RHS': The value of the problem's righthand side, correesponding
                   to the respective constraint.
            'Saturation': The saturation of the respective constraint ('LHS'/'RHS').
        """
        if isinstance(constraints, list):
            if len(constraints) > 0:
                ConstraintsInQuestion=constraints
            else:
                ConstraintsInQuestion = self.Problem.LP.row_names
        elif isinstance(constraints, str):
            ConstraintsInQuestion=[constraints]

        rhs = self.Problem.get_right_hand_side(ConstraintsInQuestion)
        lhs = self.Problem.calculate_left_hand_side(ConstraintsInQuestion)
        Out = pandas.DataFrame(columns=['LHS', 'RHS', 'Saturation'], index=ConstraintsInQuestion)
        for i in ConstraintsInQuestion:
            lhval = lhs[i]
            rhval = rhs[i]
            sat = numpy.nan
            if rhval != 0:
                sat = lhval/rhval
            Out.loc[i, 'LHS'] = lhval
            Out.loc[i, 'RHS'] = rhval
            Out.loc[i, 'Saturation'] = sat
        return(Out)

    def add_parameter_multiplier(self, model_parameter: str,rebuild_model: bool = True):
        """
        Adds a multiplicative coeffiecient for specified model parameter,
        named after provided model_parameter, followed by '_multiplier'.
        This coefficient might be changed from 1 (neutral) to any multiplicative
        term of parameter value.

        Parameters
        ----------
        model_parameter : str
            Model parameter, for which a constant multiplier should be added
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """
        if str(model_parameter+'_multiplier') not in list(self.model.parameters.functions._elements_by_id.keys()):
            if model_parameter in self.model.parameters.functions._elements_by_id.keys():
                self.model.parameters.functions._elements_by_id[model_parameter].id = str(model_parameter+'_original_definition')
                self.model.parameters.functions._elements_by_id[str(model_parameter+'_original_definition')] = self.model.parameters.functions._elements_by_id.pop(model_parameter)
                self.model.parameters.functions.append(rba.xml.parameters.Function(str(model_parameter+'_multiplier'), 'constant', parameters={'CONSTANT': 1.0}, variable=None))
                self.model.parameters.functions._elements = list(self.model.parameters.functions._elements_by_id.values())
                if model_parameter not in self.model.parameters.aggregates._elements_by_id.keys():
                    self.model.parameters.aggregates.append(rba.xml.parameters.Aggregate(model_parameter, 'multiplication'))
                    self.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(rba.xml.parameters.FunctionReference(str(model_parameter+'_original_definition')))
                    self.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(rba.xml.parameters.FunctionReference(str(model_parameter+'_multiplier')))
                else:
                    for funct in self.model.parameters.aggregates._elements_by_id[model_parameter].function_references:
                        if funct.function == model_parameter:
                            funct.function=str(model_parameter+'_original_definition')
                    self.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(rba.xml.parameters.FunctionReference(str(model_parameter+'_multiplier')))
            elif model_parameter in self.model.parameters.aggregates._elements_by_id.keys():
                self.model.parameters.functions.append(rba.xml.parameters.Function(str(model_parameter+'_multiplier'), 'constant', parameters={'CONSTANT': 1.0}, variable=None))
                self.model.parameters.aggregates._elements_by_id[model_parameter].function_references.append(rba.xml.parameters.FunctionReference(str(model_parameter+'_multiplier')))
            if rebuild_model:
                self.rebuild_from_model()

    def add_parameter_multipliers_for_enzyme_efficiencies(self,enzymes_list: list = [],default_parameter_ids: list = ["default_efficiency","default_transporter_efficiency"]) -> list:
        """
        Adds a multiplicative coeffiecient for  forward- and backward efficiencies
        of provided enzymes. This coefficient might be changed from 1 (neutral)
        to any multiplicative term of parameter value.

        Parameters
        ----------
        enzymes_list : list
            Model parameter, for which a constant multiplier should be added
            Default: [] --> then all enzymes are used
        default_parameter_ids : list
            Specifies the ids of default efficiencies
            Default: ["default_efficiency","default_transporter_efficiency"]
        Returns
        -------
        List of enzyme IDs, for which the multipliers have been added.
        """
        out=[]
        if not enzymes_list:
            enzyme_ids_to_add=[i.id for i in self.model.enzymes.enzymes._elements]
        else:
            enzyme_ids_to_add=enzymes_list
        for enzyme_id in enzyme_ids_to_add:
            if (self.model.enzymes.enzymes._elements_by_id[enzyme_id].forward_efficiency in default_parameter_ids) or (self.model.enzymes.enzymes._elements_by_id[enzyme_id].backward_efficiency in default_parameter_ids):
                self.add_specific_kapps_for_default_kapp_enzyme(enzyme=enzyme_id, default_parameter_ids= default_parameter_ids,rebuild_model= False)
            enzyme=self.model.enzymes.enzymes._elements_by_id[enzyme_id]
            self.add_parameter_multiplier(enzyme.forward_efficiency,rebuild_model=False)
            self.add_parameter_multiplier(enzyme.backward_efficiency,rebuild_model=False)
            out.append(enzyme_id)
        self.rebuild_from_model()
        return(out)

    def set_parameter_value(self, model_parameter: str, function_parameter: str, parameter_type: str='', change_message: str='', new_value : float = numpy.nan,rebuild_model: bool = True):
        """
        Sets function parameter of model parameter to specified value.

        Parameters
        ----------
        model_parameter : str
            Model parameter for which a function parameter should be changed
        function_parameter : str
            Function parameter to be changed
        new_value : float
            New parameter value to be applied
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """
        if numpy.isfinite(new_value):
            if model_parameter in self.model.parameters.functions._elements_by_id.keys():
                old_value = self.model.parameters.functions._elements_by_id[model_parameter].parameters._elements_by_id[function_parameter].value
                self.model.parameters.functions._elements_by_id[model_parameter].parameters._elements_by_id[function_parameter].value = new_value
                if rebuild_model:
                    self.rebuild_from_model()

    def set_parameter_multiplier(self, model_parameter: str, parameter_type: str='', new_value: float = numpy.nan,rebuild_model: bool = True):
        """
        Sets multiplicative factor of model parameter to specified value.
        Requires previous addition of multiplier to parameter, via add_parameter_multiplier method.

        Parameters
        ----------
        model_parameter : str
            Model parameter for which a function parameter should be changed
        new_value : float
            New parameter value to be applied
            Default: numpy.nan
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """
        self.set_parameter_value(model_parameter=str(model_parameter+'_multiplier'),
                                 function_parameter='CONSTANT',
                                 parameter_type=parameter_type,
                                 change_message='multiplicative change',
                                 new_value=new_value,rebuild_model=False)
        if rebuild_model:
            self.rebuild_from_model()

    def add_specific_kapps_for_default_kapp_enzymes(self,enzymes_list: list = [],default_parameter_ids: list = ["default_efficiency","default_transporter_efficiency"],rebuild_model: bool = True) -> list:
        """
        Adds  a specific efficiency parameter for all enzymes provided as input list of IDs,
        which have  the default efficiency assigned as efficiency parameter.
        Parameters
        ----------
        enzymes_list : list
            Model parameter, for which a constant multiplier should be added
            Default: [] --> then all enzymes are used
        default_parameter_ids : list
            Specifies the ids of default efficiencies
            Default: ["default_efficiency","default_transporter_efficiency"]
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        Returns
        -------
        List of enzyme IDs, for specific efficiencies have been added.
        """
        out=[]
        if not enzymes_list:
            enzyme_ids_to_add=[i.id for i in self.model.enzymes.enzymes._elements]
        else:
            enzyme_ids_to_add=enzymes_list
        for enzyme_id in enzyme_ids_to_add:
            if (self.model.enzymes.enzymes._elements_by_id[enzyme_id].forward_efficiency in default_parameter_ids) or (self.model.enzymes.enzymes._elements_by_id[enzyme_id].backward_efficiency in default_parameter_ids):
                self.add_specific_kapps_for_default_kapp_enzyme(enzyme=enzyme_id, default_parameter_ids=default_parameter_ids,rebuild_model= False)
                out.append(enzyme_id)
        if rebuild_model:
            self.rebuild_from_model()
        return(out)

    def add_specific_kapps_for_default_kapp_enzyme(self,enzyme: str ="", default_parameter_ids: list = ["default_efficiency","default_transporter_efficiency"],rebuild_model: bool = True):
        """
        Adds  a specific efficiency  parameter for enzymes, which are parameterised with the default efficiency.
        Specific parameter is named 'ENZYME_ID_forward_efficiency'/'ENZYME_ID_backward_efficiency' and initially
        defined exactly like the  default parameter.
        Parameters
        ----------
        enzyme : str
            ID of model enzyme
        default_parameter_ids : list
            Specifies the ids of default efficiencies
            Default: ["default_efficiency","default_transporter_efficiency"]
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """
        enzyme_structure=self.model.enzymes.enzymes._elements_by_id[enzyme]
        forward_efficiency_param=enzyme_structure.forward_efficiency
        backward_efficiency_param=enzyme_structure.backward_efficiency
        if forward_efficiency_param in default_parameter_ids:
            if forward_efficiency_param in self.model.parameters.functions._elements_by_id.keys():
                param_function=copy.deepcopy(self.model.parameters.functions._elements_by_id[forward_efficiency_param])
                param_function.id=enzyme+"_forward_efficiency"
                self.model.parameters.functions.append(param_function)
            elif forward_efficiency_param in self.model.parameters.aggregates._elements_by_id.keys():
                param_aggregate=copy.deepcopy(self.model.parameters.aggregates._elements_by_id[forward_efficiency_param])
                param_aggregate.id=enzyme+"_forward_efficiency"
                self.model.parameters.aggregates.append(param_function)
            enzyme_structure.forward_efficiency=enzyme+"_forward_efficiency"
        if backward_efficiency_param in default_parameter_ids:
            if backward_efficiency_param in self.model.parameters.functions._elements_by_id.keys():
                param_function=copy.deepcopy(self.model.parameters.functions._elements_by_id[backward_efficiency_param])
                param_function.id=enzyme+"_backward_efficiency"
                self.model.parameters.functions.append(param_function)
            elif backward_efficiency_param in self.model.parameters.aggregates._elements_by_id.keys():
                param_aggregate=copy.deepcopy(self.model.parameters.aggregates._elements_by_id[backward_efficiency_param])
                param_aggregate.id=enzyme+"_backward_efficiency"
                self.model.parameters.aggregates.append(param_function)
            enzyme_structure.backward_efficiency=enzyme+"_backward_efficiency"
        if rebuild_model:
            self.rebuild_from_model()

    def set_enzyme_forward_efficiency_multipliers(self,efficiency_multiplier_values: pandas.core.frame.DataFrame = pandas.DataFrame(),rebuild_model: bool = True):
        """
        Sets multiplicative factors of forward efficiency parameters of enzymes
        to specified values. Requires previous addition of multipliers to enzyme
        forward efficiencies, via add_parameter_multiplier method.

        Parameters
        ----------
        efficiency_multiplier_values : pandas.DataFrame
            Dataframe with mandatory columns 'Enzyme' and 'Forward_efficiency_multiplier'
            Default: {}
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """
        changes={}
        for enzyme_id in efficiency_multiplier_values["Enzyme"]:
            forward_efficiency_parameter=self.model.enzymes.enzymes._elements_by_id[enzyme_id].forward_efficiency
            self.set_parameter_multiplier(model_parameter=forward_efficiency_parameter,
                                          new_value= efficiency_multiplier_values.loc[efficiency_multiplier_values["Enzyme"]==enzyme_id,"Forward_efficiency_multiplier"].values[0],rebuild_model= False)
        if rebuild_model:
            self.rebuild_from_model()

    def set_enzyme_backward_efficiency_multipliers(self,efficiency_multiplier_values: pandas.core.frame.DataFrame = pandas.DataFrame(),rebuild_model: bool = True):
        """
        Sets multiplicative factors of backward efficiency parameters of enzymes
        to specified values. Requires previous addition of multipliers to enzyme
        backward efficiencies, via add_parameter_multiplier method.

        Parameters
        ----------
        efficiency_multiplier_values : pandas.DataFrame
            Dataframe with mandatory columns 'Enzyme' and 'Backward_efficiency_multiplier'
            Default: {}
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """
        changes={}
        for enzyme_id in efficiency_multiplier_values["Enzyme"]:
            backward_efficiency_parameter=self.model.enzymes.enzymes._elements_by_id[enzyme_id].backward_efficiency
            self.set_parameter_multiplier(model_parameter=backward_efficiency_parameter,
                                          new_value= efficiency_multiplier_values.loc[efficiency_multiplier_values["Enzyme"]==enzyme_id,"Backward_efficiency_multiplier"].values[0],rebuild_model= False)
        if rebuild_model:
            self.rebuild_from_model()

    def get_constant_enzyme_efficiencies(self,enzymes_list: list =[],exclude_default: bool = False, default_parameter_ids: list = ["default_efficiency","default_transporter_efficiency"]) -> pandas.core.frame.DataFrame:
        """
        Returns pandas DataFrame with efficiency values of enzymes, which are of parameter type constant.

        Parameters
        ----------
        enzymes_list : list
            Model parameter, for which a constant multiplier should be added
            Default: [] --> then all enzymes are used
        exclude_default : bool
            If True enzymes without a specific enzyme capacity are ignored.
            Default: False
        default_parameter_ids : list
            Specifies the ids of default efficiencies
            Default: ["default_efficiency","default_transporter_efficiency"]
        Returns
        -------
        pandas.DataFrame with columns: 'Enzyme','Forward_efficiency','Backward_efficiency'
        """

        out=pandas.DataFrame(columns=["Enzyme","Forward_efficiency","Backward_efficiency"])
        if not enzymes_list:
            for enz_id in self.model.enzymes.enzymes._elements_by_id.keys():
                enzyme=self.model.enzymes.enzymes._elements_by_id[enz_id]
                enzyme_ID=enzyme.id
                forward_efficiency_parameter=enzyme.forward_efficiency
                backward_efficiency_parameter=enzyme.backward_efficiency
                write_forward=True
                if forward_efficiency_parameter in self.model.parameters.functions._elements_by_id.keys():
                    if exclude_default:
                        if forward_efficiency_parameter in default_parameter_ids:
                            write_forward=False
                    if write_forward:
                        if self.model.parameters.functions._elements_by_id[forward_efficiency_parameter].type=="constant":
                            val=self.model.parameters.functions._elements_by_id[forward_efficiency_parameter].parameters._elements[0].value
                            out.loc[enzyme_ID,"Enzyme"]=enzyme_ID
                            out.loc[enzyme_ID,"Forward_efficiency"]=val
                write_backward=True
                if backward_efficiency_parameter in self.model.parameters.functions._elements_by_id.keys():
                    if exclude_default:
                        if backward_efficiency_parameter in default_parameter_ids:
                            write_backward=False
                    if write_backward:
                        if self.model.parameters.functions._elements_by_id[backward_efficiency_parameter].type=="constant":
                            val=self.model.parameters.functions._elements_by_id[backward_efficiency_parameter].parameters._elements[0].value
                            out.loc[enzyme_ID,"Enzyme"]=enzyme_ID
                            out.loc[enzyme_ID,"Backward_efficiency"]=val
        else:
            for enz_id in enzymes_list:
                enzyme=self.model.enzymes.enzymes._elements_by_id[enz_id]
                enzyme_ID=enzyme.id
                forward_efficiency_parameter=enzyme.forward_efficiency
                backward_efficiency_parameter=enzyme.backward_efficiency
                write_forward=True
                if forward_efficiency_parameter in self.model.parameters.functions._elements_by_id.keys():
                    if exclude_default:
                        if forward_efficiency_parameter in default_parameter_ids:
                            write_forward=False
                    if write_forward:
                        if self.model.parameters.functions._elements_by_id[forward_efficiency_parameter].type=="constant":
                            val=self.model.parameters.functions._elements_by_id[forward_efficiency_parameter].parameters._elements[0].value
                            out.loc[enzyme_ID,"Enzyme"]=enzyme_ID
                            out.loc[enzyme_ID,"Forward_efficiency"]=val
                write_backward=True
                if backward_efficiency_parameter in self.model.parameters.functions._elements_by_id.keys():
                    if exclude_default:
                        if backward_efficiency_parameter in default_parameter_ids:
                            write_backward=False
                    if write_backward:
                        if self.model.parameters.functions._elements_by_id[backward_efficiency_parameter].type=="constant":
                            val=self.model.parameters.functions._elements_by_id[backward_efficiency_parameter].parameters._elements[0].value
                            out.loc[enzyme_ID,"Enzyme"]=enzyme_ID
                            out.loc[enzyme_ID,"Backward_efficiency"]=val
        return(out)

    def set_constant_enzyme_efficiencies(self,efficiencies: pandas.core.frame.DataFrame):
        """
        Sets constant enzyme efficiencies (forward and backward) to specified values

        Parameters
        ----------
        efficiencies : pandas.core.frame.DataFrame
            Input DataFrame with columns: 'Enzyme','Forward_efficiency','Backward_efficiency'
        """
        for enz_id in efficiencies["Enzyme"]:
            enzyme=self.model.enzymes.enzymes._elements_by_id[enz_id]
            forward_efficiency_parameter=enzyme.forward_efficiency
            backward_efficiency_parameter=enzyme.backward_efficiency
            newFW=efficiencies.loc[efficiencies["Enzyme"]==enz_id,"Forward_efficiency"].values[0]
            newBW=efficiencies.loc[efficiencies["Enzyme"]==enz_id,"Backward_efficiency"].values[0]
            if forward_efficiency_parameter in self.model.parameters.functions._elements_by_id.keys():
                if self.model.parameters.functions._elements_by_id[forward_efficiency_parameter].type=="constant":
                    if not pandas.isna(newFW):
                        self.set_parameter_value(model_parameter=forward_efficiency_parameter, function_parameter="CONSTANT", new_value=newFW,rebuild_model=False)
            if backward_efficiency_parameter in self.model.parameters.functions._elements_by_id.keys():
                if self.model.parameters.functions._elements_by_id[backward_efficiency_parameter].type=="constant":
                    if not pandas.isna(newBW):
                        self.set_parameter_value(model_parameter=backward_efficiency_parameter, function_parameter="CONSTANT", new_value=newBW,rebuild_model=False)
        self.rebuild_from_model()

    def screen_multipliers(self,parameter: Union[list,str],factors: list = [],precision: float = 0.001, max: float = 4.0, start_value: float = 4.0,Variables_to_record: list =[]) -> dict:
        """
        Applies several multiplicative  factors to model parameters.
        Records corresponding maximum groth-rate and specified solution values.
        Requires previous addition of multipliers to parameters, via add_parameter_multiplier method.

        Parameters
        ----------
        parameter : str or list
            parameter ID or list of parameters to screen
        factors : list
            List of multiplicative factors to screen
            Default: []
        Variables_to_record : list
            List of problem variables' solution values to be exported .
            Default: []
        Returns
        -------
        Dictionary: with factor values as keys and dictionary with growth-rate ["Mu"] and Variables_to_record as values.
        """
        out={}
        for factor in factors:
            add_dict={i:[] for i in list(["Mu"]+Variables_to_record)}
            if type(parameter) is str:
                self.set_parameter_multiplier(model_parameter=parameter, new_value=factor)
            elif type(parameter) is list:
                for par in parameter:
                    self.set_parameter_multiplier(model_parameter=par, new_value=factor)
            add_dict["Mu"]=self.find_max_growth_rate(precision=precision, max=max, start_value=start_value)
            for j in Variables_to_record:
                if j in self.Problem.SolutionValues.keys():
                    add_dict[j].append(self.Problem.SolutionValues[j])
            out[factor]=add_dict
        return(out)

    def sample_kapp_multipliers(self,n: int = 1,mean: float = 0, stdev:float = 0.561,enzymes: list = [],Variables_to_record: list =[],wt_growth_rate: float=0.0) -> dict:
        """
        Applies n random multiplicative factors to each individual enzymes
        forward and backward capacity and maximises growth-rate.
        The multiplicative factors  are represented as e^^x,
        where x is drawn from a normal distribution. Records corresponding
        maximum groth-rate and specified solution values.
        Requires previous addition of multipliers to parameters, via add_parameter_multiplier method.

        Parameters
        ----------
        n : int
            Number of samples
            Default: 1
        mean : float
            Mean of normal distribution to draw from
            Default: 0
        stdev : float
            Stdev of normal distribution to draw from
            Default: 0.561
        enzymes : list
            List of enzymes to sample efficiencies for.
            Default: [] --> all model enzymes are used.
        Variables_to_record : list
            List of problem variables' solution values to be exported .
            Default: []
        wt_growth_rate : float
            Wild-type growth-rate to use as starting value for mumax-search.
            Default: 0.0 --> Then not set as starting value.

        Returns
        -------
        Dictionary: With growth-rate ["Mu"] and Variables_to_record as keys and list of sampled solution values as values.
        """
        out={i:[] for i in list(["Mu"]+Variables_to_record)}
        if enzymes:
            multipliers=pandas.DataFrame(index=enzymes)
            multipliers["Enzyme"]=enzymes
        else:
            multipliers=pandas.DataFrame(index=self.get_enzymes())
            multipliers["Enzyme"]=self.get_enzymes()
        for i in range(n):
            multipliers["Forward_efficiency_multiplier"]=[numpy.e**j for j in list(numpy.random.normal(loc=mean,scale=stdev,size=multipliers.shape[0]))]
            multipliers["Backward_efficiency_multiplier"]=[numpy.e**j for j in list(numpy.random.normal(loc=mean,scale=stdev,size=multipliers.shape[0]))]
            self.set_enzyme_forward_efficiency_multipliers(efficiency_multiplier_values=multipliers,rebuild_model=False)
            self.set_enzyme_backward_efficiency_multipliers(efficiency_multiplier_values=multipliers,rebuild_model=False)
            self.rebuild_from_model()
            if wt_growth_rate==0.0:
                mumax=self.find_max_growth_rate()
            else:
                mumax=self.find_max_growth_rate(start_value=wt_growth_rate)
            out["Mu"].append(mumax)
            for j in Variables_to_record:
                if j in self.Problem.SolutionValues.keys():
                    out[j].append(self.Problem.SolutionValues[j])
                else:
                    out[j].append(0.0)
        return(out)

    def local_sensitivity(self,parameters: Union[list,str] = [],relative_parameter_difference: float = 0.01,muapprox_precision: float = 0.0000001,musearch_max: float = 1.0) -> pandas.core.frame.DataFrame:
        """
        Determines local sensitivity of (max) growth-rate towards small parameter changes.
        Parameter value is increased and reduced by small relative values
        and the respective maximum growth-rate is determined. From the
        resulting differences in parameter value and growth-rate,
        the local sensitivity is inferred. In theory the local sensitivity
        represents the partial derivative of growth-rate vs parameter, at the wild-type.
        Beware that numerical precision is an issue  with this method and can influence results.

        Parameters
        ----------
        parameters : list or string
            Parameters to apply local sensitivity ananlysis to.
            If list with model parameter IDs is provided, those are used.
            If str 'Enzymes' all enzyme capacities are used.
            If str 'Processes' all process capacities are used.
            If str 'Targets' all cellular targets are used.
            If str 'Densities' all compartment capacities are used.
            If empty list, all model cellular targets, enzyme- process- and compartment capacities are used.
            Default: []
        relative_parameter_difference : float
            Parameter is multplied once with factor (1+relative_parameter_difference) and once with factor (1-relative_parameter_difference).
            The resulting maximum growth rates are determined.
            Default: e-2
        muapprox_precision : float
            Numeric precision to which the maximum growth rate is determined. (See parameter 'precision' in method find_max_growth_rate)
            Default: e-8
        muapprox_precision : float
            Upper bound of search space for maximum growth rate. (See parameter 'max' in method find_max_growth_rate)
            Default: 1.0
        Returns
        -------
        Pandas DataFrame with results.
        Columns:
            'Mu_WT' : Maximum wild-type growth-rate.
            'WT_param_value' : original parameter value at maximum wild-type growth-rate.
            'Relative_param_change' : 'relative_parameter_difference' input parameter
            'Absolute_param_change' : 'relative_parameter_difference' input parameter times original parameter value at maximum wild-type growth-rate.
            'Upper_Mu' : Maximum growth rate at parameter multiplied with factor (1+relative_parameter_difference)
            'Lower_Mu' : Maximum growth rate at parameter multiplied with factor (1-relative_parameter_difference)
            'Mu_change' : Difference between Upper- and Lower_Mu
            'Absolute_Sensitivity' : Absolute difference in maximum growth-rate per absolute change in parameter-value ((d_mu)/(2*out.loc[param,"Absolute_param_change"])).
            'Scaled_Sensitivity' : Relative change in maximum growth-rate (to wild-type) per relative change in parameter-value ((d_mu/WT_mu)/(2*relative_parameter_difference)).
        """
        WT_mu=self.find_max_growth_rate(precision=muapprox_precision,max=musearch_max)
        out=pandas.DataFrame()
        compartment_densities=[i.upper_bound for i in self.model.density.target_densities._elements if i.upper_bound is not None]
        process_efficiencies=[i.machinery.capacity.value for i in self.model.processes.processes._elements if i.machinery.capacity.value is not None]
        targets=[self.ModelStructure.TargetInfo.Elements[i]["TargetParameterID"] for i in self.ModelStructure.TargetInfo.Elements.keys()]
        densities=[i.upper_bound for i in self.model.density.target_densities._elements]
        enzyme_efficiencies=[]
        for i in self.model.enzymes.enzymes._elements:
            enzyme_efficiencies+=[i.forward_efficiency,i.backward_efficiency]

        accepted_parameters=list(set(compartment_densities+process_efficiencies+enzyme_efficiencies+targets))

        if isinstance(parameters, list):
            if len(parameters)>0:
                parameters_to_test=parameters
            else:
                parameters_to_test=accepted_parameters
        elif isinstance(parameters, str):
            if parameters == "Enzymes":
                parameters_to_test=list(set(enzyme_efficiencies))
            elif parameters == "Targets":
                parameters_to_test=list(set(targets))
            elif parameters == "Processes":
                parameters_to_test=list(set(process_efficiencies))
            elif parameters == "Densities":
                parameters_to_test=list(set(densities))

        for param in parameters_to_test:
            self.add_parameter_multiplier(model_parameter=param,rebuild_model=True)
            out.loc[param,"Mu_WT"]=WT_mu
            out.loc[param,"WT_param_value"]=self.get_current_parameter_value(parameter=param)
            out.loc[param,"Relative_param_change"]=relative_parameter_difference
            out.loc[param,"Absolute_param_change"]=relative_parameter_difference*out.loc[param,"WT_param_value"]
            self.add_parameter_multiplier(model_parameter=param,rebuild_model=False)
            screen=self.screen_multipliers(parameter=param,factors=[1-relative_parameter_difference,1+relative_parameter_difference],precision=muapprox_precision, max=musearch_max, start_value=WT_mu)
            lower_mu=screen[1-relative_parameter_difference]["Mu"]
            upper_mu=screen[1+relative_parameter_difference]["Mu"]
            d_mu=upper_mu-lower_mu
            out.loc[param,"Upper_Mu"]=upper_mu
            out.loc[param,"Lower_Mu"]=lower_mu
            out.loc[param,"Mu_change"]=d_mu
            out.loc[param,"Absolute_Sensitivity"]=(d_mu)/(2*out.loc[param,"Absolute_param_change"])
            out.loc[param,"Scaled_Sensitivity"]=(d_mu/WT_mu)/(2*relative_parameter_difference)
            self.set_parameter_multiplier(model_parameter=param, new_value=1.0)
        return(out)

    def get_current_parameter_value(self, parameter: str) -> float:
        """
        Return current value of model parameter (depending on currently
        set growth-rate and medium-composition).

        Parameters
        ----------
        parameter : str
            ID of model parameter

        Returns
        -------
        Float: Parameter value
        """
        try:
            parameter_definition=self.get_parameter_definition(parameter=parameter)
            if parameter_definition[parameter]['Type']=='Aggregate':
                list_param_definitions_to_evaluate=[self.get_parameter_definition(parameter=i) for i in parameter_definition[parameter]['Multiplicative Terms']]
            else:
                list_param_definitions_to_evaluate=[parameter_definition]
            term_results=[_auxiliary_functions.evaluate_current_function_value(expression_dictionary=i,growth_rate=self.Mu , medium=self.Medium) for i in list_param_definitions_to_evaluate]
            return(numpy.prod(term_results))
        except:
            return(numpy.nan)

    def get_parameter_definition(self, parameter: str) -> dict:
        """
        Returns definition of model parameter.

        Parameters
        ----------
        parameter : str
            ID of model parameter

        Returns
        -------
        Dictionary with parameter definition information: {parameter:info}
        The info value is a dictionary itself, which depends on the whether
        the parameter is among aggrgates or functions.
        info={'Type' : str, --> type of parameter  "Aggrgate" if aggregate and function type (michaelisMenten, linear, constant...) if function
              'Equation' : str, --> mathematical equation of parameter definition
              'Variables' : list, --> list of independent variables of paramter definition (growth-rate or medium components)
              'Function_parameters' : list, --> list of function parameters in function definition.
              'Multiplicative Terms' : list, --> list of multiplicative terms (function IDs) composing the parameter (several function ids if aggregate, only own id if function)
              'Generic_latex' : str, --> equation in LaTex format (with variable ids for function parameters)
              'Specific_latex' : str, --> equation in LaTex format (with actual values for function parameters)
            }
        """
        return(_auxiliary_functions.return_parameter_definition(model=self.model,parameter=parameter))

    def get_parameter_evolution(self,model_parameter:str,x_values:dict={"growth_rate":[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]}) -> pandas.core.frame.DataFrame:
        """
        Returns the evolution of model parameter value over range of provided independent variable values.

        Parameters
        ----------
        model_parameter : str
            ID of model parameter, whos values should be determined over x_values.
        x_values : dict
            Dictionary with possible independent variables ('growth_rate' or medium components) as keys and list of respective values as values
            Default: {"growth_rate":[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]}
        Returns
        -------
        pandas.core.frame.DataFrame with columns for independent variable (growth-rate or medium component concentration) 
        and respective parameter values.
        """
        parameter_definition=self.get_parameter_definition(model_parameter)
        medium=self.Medium.copy()
        if len(parameter_definition[model_parameter]['Variables'])==1:
            if parameter_definition[model_parameter]['Variables'][0] in x_values.keys():
                if parameter_definition[model_parameter]['Variables'][0]=="growth_rate":
                    y_values=[]
                    for mu in x_values["growth_rate"]:
                        if parameter_definition[model_parameter]['Type']=='Aggregate':
                            list_param_definitions_to_evaluate=[self.get_parameter_definition(parameter=i) for i in parameter_definition[model_parameter]['Multiplicative Terms']]
                        else:
                            list_param_definitions_to_evaluate=[parameter_definition]
                        y_values.append(numpy.prod([_auxiliary_functions.evaluate_current_function_value(expression_dictionary=i,growth_rate=mu , medium=medium) for i in list_param_definitions_to_evaluate]))
                    out=pandas.DataFrame()
                    out["growth_rate"]=x_values["growth_rate"]
                    out[model_parameter]=y_values
                    return(out)
                else:
                    if parameter_definition[model_parameter]['Variables'][0] in medium.keys():
                        y_values=[]
                        for i in x_values[parameter_definition[model_parameter]['Variables'][0]]:
                            medium[parameter_definition[model_parameter]['Variables'][0]]=i
                            if parameter_definition[model_parameter]['Type']=='Aggregate':
                                list_param_definitions_to_evaluate=[self.get_parameter_definition(parameter=i) for i in parameter_definition[model_parameter]['Multiplicative Terms']]
                            else:
                                list_param_definitions_to_evaluate=[parameter_definition]
                            y_values.append(numpy.prod([_auxiliary_functions.evaluate_current_function_value(expression_dictionary=i,growth_rate=self.Mu , medium=medium) for i in list_param_definitions_to_evaluate]))
                        out=pandas.DataFrame()
                        out[str(parameter_definition[model_parameter]['Variables'][0])]=x_values[parameter_definition[model_parameter]['Variables'][0]]
                        out[model_parameter]=y_values
                        return(out)
                    else:
                        raise Exception('X variable not found')
            else:
                raise Exception("Variable not found in x values")
        else:
            print("WARNING: Parameter dependent on more than one variable")
            input_variable_independent=list(x_values.keys())[0]
            if input_variable_independent in parameter_definition[model_parameter]['Variables']:
                if input_variable_independent=="growth_rate":
                    y_values=[]
                    for mu in x_values["growth_rate"]:
                        if parameter_definition[model_parameter]['Type']=='Aggregate':
                            list_param_definitions_to_evaluate=[self.get_parameter_definition(parameter=i) for i in parameter_definition[model_parameter]['Multiplicative Terms']]
                        else:
                            list_param_definitions_to_evaluate=[parameter_definition]
                        y_values.append(numpy.prod([_auxiliary_functions.evaluate_current_function_value(expression_dictionary=i,growth_rate=mu , medium=medium) for i in list_param_definitions_to_evaluate]))
                    out=pandas.DataFrame()
                    out["growth_rate"]=x_values["growth_rate"]
                    out[model_parameter]=y_values
                    return(out)
                elif input_variable_independent in medium.keys():
                    y_values=[]
                    for i in x_values[parameter_definition[model_parameter]['Variables'][0]]:
                        medium[parameter_definition[model_parameter]['Variables'][0]]=i
                        if parameter_definition[model_parameter]['Type']=='Aggregate':
                            list_param_definitions_to_evaluate=[self.get_parameter_definition(parameter=i) for i in parameter_definition[model_parameter]['Multiplicative Terms']]
                        else:
                            list_param_definitions_to_evaluate=[parameter_definition]
                        y_values.append(numpy.prod([_auxiliary_functions.evaluate_current_function_value(expression_dictionary=i,growth_rate=self.Mu , medium=medium) for i in list_param_definitions_to_evaluate]))
                    out=pandas.DataFrame()
                    out[input_variable_independent]=x_values[input_variable_independent]
                    out[model_parameter]=y_values
                    return(out)
            else:
                raise Exception('X variable not found')

    def get_model_statistics_information(self) -> dict:
        """
        Get model statistics from ModelStructure

        Returns
        -------
        Dictionary with information
        """
        return(_auxiliary_functions.return_model_statistics_information(model_structure=self.ModelStructure))

    def get_general_model_information(self) -> dict:
        """
        Get general model information from ModelStructure

        Returns
        -------
        Dictionary with information on model name, organism , author...
        """
        return(_auxiliary_functions.return_general_model_information(model_structure=self.ModelStructure))

    def get_compartment_information(self,compartment: str) -> dict:
        """
        Get information on model compartment from ModelStructure

        Parameters
        ----------
        compartment : str
            ID of compartment.

        Returns
        -------
        Dictionary with information:
            'ID' : compartment ID in model (type str)
            'associatedProteins' : proteins localised to compartment (type list)
            'associatedEnzymes' : enzyme located in compartment (type list)
            'associatedReactions' : metabolic reactions located in compartment (type list)
            'associatedMacromolecules' : other macromolecules than proteins (DNA,RNA), located in compartment (type list)
            'Capacity_Constraint' : id of corresponding density constraint.
        """
        return(_auxiliary_functions.return_compartment_information(model_structure=self.ModelStructure,compartment=compartment))

    def get_enzyme_information(self,enzyme: str) -> dict:
        """
        Get information on model enzyme from ModelStructure

        Parameters
        ----------
        enzyme : str
            ID of enzyme.

        Returns
        -------
        Dictionary with information:
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
        return(_auxiliary_functions.return_enzyme_information(model_structure=self.ModelStructure,enzyme=enzyme))

    def get_protein_information(self,protein: str) -> dict:
        """
        Get information on model protein from ModelStructure

        Parameters
        ----------
        protein : str
            ID of protein.

        Returns
        -------
        Dictionary with information:
           'ID' : protein ID in model (type str). May be compartment_isoform.
           'ProtoID' : protein ID in model (type str). Equals ID without compartment-specific ending.
           'ExternalIDs' : identifiers of this protein in other namespaces (Locus-tag, Gene-symbol, ...) (type dict)
           'Function' : associated function, according to Uniprot (type str)
           'Name' : name, according to Uniprot (type list)
           'associatedReactions' : Reactions, this protein is a subunit of (type list)
           'associatedEnzymes' : Enzymes, this protein is a subunit of (type list)
           'Compartment' : Location of the protein (type str)
           'AAcomposition' : Which amino-acids the protein is composed of and what is the stoichiometry(type dict)
           'AAnumber' : Length of protein (type int)
           'Weight' : Weight of protein (type float)
           'ProcessRequirements' : Which processes the protein requires for synthesis and maintenance and how much (type dict)
           'SupportsProcess' : Process-machineries, this protein is a subunit of (type list)
           'AssociatedTarget' : Wheter protein represents a (translation) target.
        """
        return(_auxiliary_functions.return_protein_information(model_structure=self.ModelStructure,protein=protein))

    def get_reaction_information(self,reaction: str) -> dict:
        """
        Get information on model reaction from ModelStructure

        Parameters
        ----------
        reaction : str
            ID of reaction.

        Returns
        -------
        Dictionary with information
           'ID' : reaction ID in model (type str)
           'Compartment_Machinery' : Localisation of enzyme-subunit, catalysing this reaction. (type list)
           'Name' : name according to BiGG (type str)
           'OtherIDs' : Other reaction names (eg. BiGG, KEGG) (type list)
           'Formula' : Reaction formula as string (type str)
           'Reactants' : Which metabolites does this reaction consume of and many (type dict)
           'Products' : Which metabolites does this reaction produce of and many (type dict)
           'Reversible' :  Wheter reaction is reversible (type bool)
           'Type' : Type of reaction ('normal' or 'transport') (type str)
           'Compartment_Species' : Location of the metabolites involved with this reaction (type list)
           'Enzyme' : Enzyme catalysing this reaction (type str)
           'Twins' : Isoreactions of this reactions (catalysed by iso-enzymes) (type list)
           'AssociatedTarget' :  Wheter reaction has a flux target.
        """
        return(_auxiliary_functions.return_reaction_information(model_structure=self.ModelStructure,reaction=reaction))

    def get_metabolite_information(self,metabolite: str) -> dict:
        """
        Get information on model metabolite from ModelStructure

        Parameters
        ----------
        metabolite : str
            ID of metabolite.

        Returns
        -------
        Dictionary with information:
          'ID' : meatbolite ID in model (type str)
          'OtherIDs' : identifiers of this metabolite in other namespaces (BiGG, KEGG ...) (type dict)
          'Name' : Name according to BiGG (type str)
          'ReactionsInvolvedWith' : Reactions which produce or consume this metabolite (type list)
          'boundary' : Boundary metabolite (type boolean)
          'Type' :  Type of metabolite (internal exernal or biomass-precursor) (type str)
          'Compartment' : Location of meatbolite (type str)
          'AssociatedTarget' :  Wheter metabolite represents a target
          'MassBalance_Constraint' : Id of associated mass-balance constraint.
        """
        return(_auxiliary_functions.return_metabolite_information(model_structure=self.ModelStructure,metabolite=metabolite))

    def get_process_information(self,process: str) -> dict:
        """
        Get information on model process from ModelStructure

        Parameters
        ----------
        process : str
            Name of process.

        Returns
        -------
        Dictionary with information:
            'ID' : process ID in model (type str)
            'Name' : name of process (same as key) (type str)
            'Initiation' : reaction-string of initiation (type list)
            'Composition' : Which proteins the process-machinery is composed of and how many (type dict)
            'Components' : Substrates to process (type dict)
            'Capacity_Constraint' :  ID of associated capacity constraint
        """
        return(_auxiliary_functions.return_process_information(model_structure=self.ModelStructure,process=process))

    def get_macro_molecule_information(self,macro_molecule: str) -> dict:
        """
        Get information on model macromolecule from ModelStructure

        Parameters
        ----------
        macro_molecule : str
            ID of macro molecule.

        Returns
        -------
        Dictionary with information:
           'ID' : macromolecule ID in model (type str).
           'ProtoID' : location independent macromolecule ID in model (type str). Equals ID without compartment-specific ending.
           'Type' : DNA or RNA (type str)
           'Compartment' : Location of the macromolecule (type str)
           'Composition' : Which building-blocs the macromolecule is composed of and what is the stoichiometry (type dict)
           'ProcessRequirements' : Which processes the macromolecule requires for synthesis and maintenance and how much (type dict)
           'SupportsProcess' : Process-machineries, this macromolecule is a subunit of (type list)
           'AssociatedTarget' : Target associated with macromolecule (type str)
        """
        return(_auxiliary_functions.return_macro_molecule_information(model_structure=self.ModelStructure,macro_molecule=macro_molecule))

    def get_target_information(self,target: str) -> dict:
        """
        Get information on model target from ModelStructure

        Parameters
        ----------
        target : str
            ID of target.

        Returns
        -------
        Dictionary with information:
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
        return(_auxiliary_functions.return_target_information(model_structure=self.ModelStructure,target=target))

    def get_module_information(self,module: str) -> dict:
        """
        Get information on model module from ModelStructure

        Parameters
        ----------
        module : str
            ID of module.

        Returns
        -------
        Dictionary with information
        """
        return(_auxiliary_functions.return_module_information(model_structure=self.ModelStructure,module=module))

    def get_density_constraint_information(self,density_constraint: str) -> dict:
        """
        Get information on model density-constraint from ModelStructure

        Parameters
        ----------
        density_constraint : str
            ID of constraint.

        Returns
        -------
        Dictionary with information:
           'ID' : density constraint ID in model (type str)
           'AssociatedCompartment' : ID of compartment this constraint defines the capacity for (type str)
           'Type' : Equality or inequality (type str)
           'CapacityParameterID' :  ID of density parameter (type str)
           'CapacityParameter' : Definition of capacity parameter (type dict)
                See doc string of rbatools.rba_session.SessionRBA.get_parameter_definition method
        """
        return(_auxiliary_functions.return_density_constraint_information(model_structure=self.ModelStructure,density_constraint=density_constraint))

    def get_process_constraint_information(self,process_constraint: str) -> dict:
        """
        Get information on model process-capacity constraint from ModelStructure

        Parameters
        ----------
        process_constraint : str
            ID of constraint.

        Returns
        -------
        Dictionary with information:
            'ID' : process-constraint ID in model (type str)
            'AssociatedProcess' : Name of process this constraint relates to (type str)
            'AssociatedProcessID' : ID of process this constraint relates to (type str)
            'Type' : Equality or inequality (type str)
            'CapacityParameterID' :  ID of capacity parameter (type str)
            'CapacityParameter' : Definition of capacity parameter (type dict)
                See doc string of rbatools.rba_session.SessionRBA.get_parameter_definition method
        """
        return(_auxiliary_functions.return_process_constraint_information(model_structure=self.ModelStructure,process_constraint=process_constraint))

    def get_metabolite_constraint_information(self,metabolite_constraint: str) -> dict:
        """
        Get information on model metabolite mass-balance constraint from ModelStructure

        Parameters
        ----------
        metabolite_constraint : str
            ID of constraint.

        Returns
        -------
        Dictionary with information:
           'ID' : metabolite ID in model (type str)
           'AssociatedMetabolite' : ID of metabolite this constraint defines the mass-balance for (same as key) (type str)
           'Type' : Equality or inequality (type dict)
        """
        return(_auxiliary_functions.return_metabolite_constraint_information(model_structure=self.ModelStructure,metabolite_constraint=metabolite_constraint))

    def get_enzyme_constraint_information(self,enzyme_constraint: str) -> dict:
        """
        Get information on model enzyme-capacity constraint from ModelStructure

        Parameters
        ----------
        enzyme_constraint : str
            ID of constraint.

        Returns
        -------
        Dictionary with information:
           'ID' : enzyme-constraint ID in model (type str)
           'AssociatedEnzyme' : ID of enzyme this constraint relates to (type str)
           'AssociatedReaction' : ID of reaction this constraint relates to (type str)
           'Direction': forward or backward efficiency (type str)
           'Type' : Equality or inequality (type str)
           'CapacityParameterID' :  ID of capacity parameter (type str)
           'CapacityParameter' : Definition of capacity parameter (type dict)
                See doc string of rbatools.rba_session.SessionRBA.get_parameter_definition method
        """
        return(_auxiliary_functions.return_enzyme_constraint_information(model_structure=self.ModelStructure,enzyme_constraint=enzyme_constraint))

    def get_compartments(self) -> list:
        """
        Get all model compartment-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_compartments(model_structure=self.ModelStructure))

    def get_enzymes(self) -> list:
        """
        Get all model enzyme-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_enzymes(model_structure=self.ModelStructure))

    def get_proteins(self) -> list:
        """
        Get all model protein-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_proteins(model_structure=self.ModelStructure))

    def get_reactions(self) -> list:
        """
        Get all model reaction-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_reactions(model_structure=self.ModelStructure))

    def get_metabolites(self) -> list:
        """
        Get all model metabolite-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_metabolites(model_structure=self.ModelStructure))

    def get_processes(self) -> list:
        """
        Get all model process-names

        Returns
        -------
        List of names
        """
        return(_auxiliary_functions.return_processes(model_structure=self.ModelStructure))

    def get_macro_molecules(self) -> list:
        """
        Get all model macromolecule-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_macro_molecules(model_structure=self.ModelStructure))

    def get_targets(self) -> list:
        """
        Get all model target-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_targets(model_structure=self.ModelStructure))

    def get_modules(self) -> list:
        """
        Get all model module-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_modules(model_structure=self.ModelStructure))

    def get_density_constraints(self) -> list:
        """
        Get all model density constraint-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_density_constraints(model_structure=self.ModelStructure))

    def get_process_constraints(self) -> list:
        """
        Get all model process-capacity constraint-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_process_constraints(model_structure=self.ModelStructure))

    def get_metabolite_constraints(self) -> list:
        """
        Get all model metabolite mass-balance constraint-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_metabolite_constraints(model_structure=self.ModelStructure))

    def get_enzyme_constraints(self) -> list:
        """
        Get all model enzyme-capacity constraint-IDs

        Returns
        -------
        List of IDs
        """
        return(_auxiliary_functions.return_enzyme_constraints(model_structure=self.ModelStructure))

    def add_reaction(self,reaction_id:str,reactants:dict={},products:dict={},reversible:bool=False,rebuild_model: bool = True):
        """
        Adds reaction to model.

        Parameters
        ----------
        reaction_id : str
            ID of reaction to add.
        reactants : dict
            Reactants of reaction (Metabolite-IDs as keys and their stoichiometries as values)
            Default: {}
        products : dict
            Products of reaction (Metabolite-IDs as keys and their stoichiometries as values)
            Default: {}
        reversible : bool
            Whether reaction is reversible
            Default: True
        rebuild_model : bool
            Wheter the model should be immediately rebuilt
            Default: True
        """

        rxn_to_add=rba.xml.metabolism.Reaction(id_=reaction_id,reversible=reversible)
        for i in reactants.keys():
            rxn_to_add.reactants.append(rba.xml.common.SpeciesReference(species=i,stoichiometry=reactants[i]))
        for i in products.keys():
            rxn_to_add.products.append(rba.xml.common.SpeciesReference(species=i,stoichiometry=products[i]))
        self.model.metabolism.reactions.append(rxn_to_add)
        if rebuild_model:
            self.rebuild_from_model()

    def derive_current_biomass_function(self,from_rba_solution=False) -> pandas.core.frame.DataFrame:
        """
        Returns a biomass-function, corresponding to the currently imposed biomass composition in the RBA-model.
        (Since biomass composition is growth-rate dependent, it changes when other growth-rates are set.)

        Parameters
        ----------
        from_rba_solution : bool
            If True: The biomass-composition, specific to the current solution 
            (growth-rate and decision-variable values) is determined.
            If False: The biomass-composition, specific to the current growth-rate 
            and the corresponding target-values is defined.
        Returns
        -------
        pandas.core.frame.DataFrame with columns "Metabolite" and "Coefficient", indicating metabolic species 
        and their respective 'stoichiometric' coefficients in generated biomass function. 

        """
        if from_rba_solution:
            if self.Problem.Solved:
                RBAproblem = copy.copy(self.Problem.LP)
                model_reactions=[i.id for i in self.model.metabolism.reactions]
                model_metabolites=[i.id for i in self.model.metabolism.species if not i.boundary_condition]
                Cols2remove =[RBAproblem.col_names.index(i) for i in RBAproblem.col_names if not i in model_reactions]
                Rows2remove = [RBAproblem.row_names.index(i) for i in RBAproblem.row_names if not i in model_metabolites]
                A = RBAproblem.A.toarray()
                S_colsremoved = numpy.delete(A, Cols2remove, axis=1)
                S = numpy.delete(S_colsremoved, Rows2remove, axis=0)
                col_names_S = list(numpy.delete(RBAproblem.col_names, Cols2remove))
                row_names_S = list(numpy.delete(RBAproblem.row_names, Rows2remove))
                reaction_fluxes=numpy.array([self.Problem.SolutionValues[i] for i in col_names_S])
                S_rhs = S.dot(reaction_fluxes)
                Metabolite_net_productions = {metabolite: S_rhs[row_names_S.index(metabolite)] for metabolite in row_names_S}
                out=pandas.DataFrame()
                for metabolite in Metabolite_net_productions.keys():
                    if Metabolite_net_productions[metabolite]!=0:
                        out.loc[metabolite,"Metabolite"] = metabolite
                        out.loc[metabolite,"Coefficient"] = -Metabolite_net_productions[metabolite]
            else:
                raise Exception('No solution to obtain biomass composition for')
        else:
            metabolite_constraints=[self.get_metabolite_constraint_information(metabolite_constraint=i)['AssociatedMetabolite'] for i in self.get_metabolite_constraints()]
            metabolite_rhs=self.Problem.get_right_hand_side(constraints = metabolite_constraints)
            out=pandas.DataFrame()
            for metabolite in metabolite_rhs.keys():
                if metabolite_rhs[metabolite]!=0:
                    out.loc[metabolite,"Metabolite"] = metabolite
                    out.loc[metabolite,"Coefficient"] = -metabolite_rhs[metabolite]
            pg_protein_chosen=False
            demand_for_processes={}
            for target in self.get_targets():
                target_info=self.get_target_information(target=target)
                if (target_info["Group"]=="maintenance_atp_target") and (target_info["Type"]=="reaction_fluxes"):
                    reaction=target_info["TargetEntity"]
                    bound_type=target_info["TargetConstraint"]
                    if bound_type=="UpperBound":
                        flux_value=self.Problem.get_ub(variables = reaction)[reaction]
                    else:
                        flux_value=self.Problem.get_lb(variables = reaction)[reaction]
                    reaction_reactants=self.get_reaction_information(reaction=reaction)["Reactants"]
                    reaction_products=self.get_reaction_information(reaction=reaction)["Products"]
                    for reactant in reaction_reactants.keys():
                        if reactant in out["Metabolite"]:
                            out.loc[reactant,"Coefficient"] -= reaction_reactants[reactant]*flux_value
                        else:
                            out.loc[reactant,"Coefficient"] = reaction_reactants[reactant]*flux_value
                    for product in reaction_products.keys():
                        if product in out["Metabolite"]:
                            out.loc[product,"Coefficient"] += reaction_products[product]*flux_value
                        else:
                            out.loc[product,"Coefficient"] += reaction_products[product]*flux_value
                elif (target_info["Group"]=="translation_targets") and (target_info["Type"]=="concentrations"):
                    for process_requirement in self.get_protein_information(protein=target_info["TargetEntity"])["ProcessRequirements"].keys():
                        if process_requirement in demand_for_processes.keys():
                            demand_for_processes[process_requirement]+=self.get_current_parameter_value(parameter=target_info["TargetParameterID"])*self.get_protein_information(protein=target_info["TargetEntity"])["AAnumber"]
                        else:
                            demand_for_processes[process_requirement]=self.get_current_parameter_value(parameter=target_info["TargetParameterID"])*self.get_protein_information(protein=target_info["TargetEntity"])["AAnumber"]
                    if not pg_protein_chosen:
                        pg_protein_chosen=True
                        general_pg_protein_size=self.get_protein_information(protein=target_info["TargetEntity"])["AAnumber"]
                        total_non_pg_density=abs(sum(list(self.Problem.get_right_hand_side(constraints = self.get_density_constraints()).values())))
                        respective_general_pg_protein_concentration=total_non_pg_density/general_pg_protein_size
                        general_pg_protein_composition=self.get_protein_information(protein=target_info["TargetEntity"])["AAcomposition"]
                        general_pg_protein_process_requirements=self.get_protein_information(protein=target_info["TargetEntity"])["ProcessRequirements"]
                        for process in general_pg_protein_process_requirements.keys():
                            if process_requirement in demand_for_processes.keys():
                                demand_for_processes[process_requirement] += respective_general_pg_protein_concentration*general_pg_protein_size
                            else:
                                demand_for_processes[process_requirement] = respective_general_pg_protein_concentration*general_pg_protein_size
                            process_component_map=self.get_process_information(process=process)["Components"]
                            for component in general_pg_protein_composition.keys():
                                if component in list(process_component_map.keys()):
                                    component_reactants=process_component_map[component]["Reactants"]
                                    component_products=process_component_map[component]["Products"]
                                    for reactant in component_reactants.keys():
                                        if reactant in out["Metabolite"]:
                                            out.loc[reactant,"Coefficient"] -= component_reactants[reactant]*respective_general_pg_protein_concentration
                                        else:
                                            out.loc[reactant,"Coefficient"] = component_reactants[reactant]*respective_general_pg_protein_concentration
                                    for product in component_products.keys():
                                        if product in out["Metabolite"]:
                                            out.loc[product,"Coefficient"] += component_products[product]*respective_general_pg_protein_concentration
                                        else:
                                            out.loc[product,"Coefficient"] += component_products[product]*respective_general_pg_protein_concentration
            for process in demand_for_processes.keys():
                nonprotein_subunits={}
                for component in self.get_process_information(process=process_requirement)["Composition"].keys():
                    if component in self.get_macro_molecules():
                        nonprotein_subunits[component]=self.get_process_information(process=process_requirement)["Composition"][component]
                if len(nonprotein_subunits.keys())!=0:
                    demand_flux=demand_for_processes[process]*self.Mu
                    process_capacity_constraint=self.get_process_information(process=process_requirement)["Capacity_Constraint"]
                    if process_capacity_constraint in self.get_process_constraints():
                        try:
                            process_efficiency=self.get_current_parameter_value(parameter=self.get_process_constraint_information(process_constraint=process_capacity_constraint)["CapacityParameterID"])
                            machinery_demand=demand_flux/process_efficiency
                        except:
                            machinery_demand=0.0
                        if machinery_demand!=0.0:
                            for nonprotein_subunit in nonprotein_subunits.keys():
                                stochiometric_factor=nonprotein_subunits[nonprotein_subunit]
                                demanded_concentration_of_nonprotein_component=machinery_demand*stochiometric_factor
                                respective_composition=self.get_macro_molecule_information(macro_molecule=nonprotein_subunit)["Composition"]
                                required_processes=self.get_macro_molecule_information(macro_molecule=nonprotein_subunit)["ProcessRequirements"]
                                for process_requirement_nonprotein in required_processes:
                                    process_component_map=self.get_process_information(process=process_requirement_nonprotein)["Components"]
                                    for component in respective_composition.keys():
                                        if component in list(process_component_map.keys()):
                                            component_reactants=process_component_map[component]["Reactants"]
                                            component_products=process_component_map[component]["Products"]
                                            for reactant in component_reactants.keys():
                                                if reactant in out["Metabolite"]:
                                                    out.loc[reactant,"Coefficient"] -= component_reactants[reactant]*demanded_concentration_of_nonprotein_component
                                                else:
                                                    out.loc[reactant,"Coefficient"] = component_reactants[reactant]*demanded_concentration_of_nonprotein_component
                                            for product in component_products.keys():
                                                if product in out["Metabolite"]:
                                                    out.loc[product,"Coefficient"] += component_products[product]*demanded_concentration_of_nonprotein_component
                                                else:
                                                    out.loc[product,"Coefficient"] += component_products[product]*demanded_concentration_of_nonprotein_component
        out["Coefficient"]/=self.Mu
        return(out)

    def build_fba_model(self,rba_derived_biomass_function=True,from_rba_solution=True):
        """
        Derives and constructs FBA-problem from the RBA-problem and stores the 
        rbatools.fba_problem.ProblemFBA object as attribute 'FBA'. 
        (By default non-parsimonious FBA-problem)

        Parameters
        ----------
        rba_derived_biomass_function : bool
            If True: A biomass function (named 'R_BIOMASS_targetsRBA'), specific to the current growth rate,
            is added to the problem.
            If False: It is assumed that there exists a biomass function, in the model.
        from_rba_solution : bool
            If True: The biomass-composition, specific to the current solution 
            (growth-rate and decision-variable values) is determined.
            If False: The biomass-composition, specific to the current growth-rate 
            and the corresponding target-values is defined.
         """
        RBAproblem = copy.copy(self.Problem.LP)
        A = RBAproblem.A.toarray()
        non_duplicate_model_reactions=[i.id for i in self.model.metabolism.reactions if not '_duplicate_' in i.id]+list(self.ExchangeReactionMap.values())
        model_metabolites=[i.id for i in self.model.metabolism.species]
        Cols2remove =[RBAproblem.col_names.index(i) for i in RBAproblem.col_names if not i in non_duplicate_model_reactions]
        Rows2remove = [RBAproblem.row_names.index(i) for i in RBAproblem.row_names if not i in model_metabolites]
        for target in self.get_targets():
            target_info=self.get_target_information(target=target)
            if (target_info["Group"]=="maintenance_atp_target") and (target_info["Type"]=="reaction_fluxes"):
                Cols2remove.append(RBAproblem.col_names.index(target_info["TargetEntity"]))
        Anew = numpy.delete(A, Cols2remove, axis=1)
        col_namesNew = list(numpy.delete(RBAproblem.col_names, Cols2remove))
        LBnew = numpy.delete(RBAproblem.LB, Cols2remove)
        UBnew = numpy.delete(RBAproblem.UB, Cols2remove)
        fNew = numpy.delete(RBAproblem.f, Cols2remove)
        Anew2 = numpy.delete(Anew, Rows2remove, axis=0)
        row_namesNew = list(numpy.delete(RBAproblem.row_names, Rows2remove))
        row_signsNew = list(numpy.delete(RBAproblem.row_signs, Rows2remove))
        if rba_derived_biomass_function:
            BMfunction=self.derive_current_biomass_function(from_rba_solution=from_rba_solution)
            col_namesNew.append('R_BIOMASS_targetsRBA')
            LBnew = numpy.append(LBnew, 0.0)
            UBnew = numpy.append(UBnew, 10000.0)
            fNew = numpy.append(fNew, 0.0)
            BMrxnCol = numpy.zeros((len(row_namesNew), 1))
            for metabolite in BMfunction.index:
                if metabolite in row_namesNew:
                    if BMfunction.loc[metabolite,"Coefficient"]!=0:
                        BMrxnCol[row_namesNew.index(metabolite), 0]=numpy.float64(BMfunction.loc[metabolite,"Coefficient"])
            Anew2 = numpy.append(Anew2, BMrxnCol, axis=1)

        Matrix1 = ProblemMatrix()
        Matrix1.A = scipy.sparse.coo_matrix(Anew2)
        Matrix1.LB = LBnew
        Matrix1.UB = UBnew
        Matrix1.row_signs = row_signsNew
        Matrix1.row_names = row_namesNew
        Matrix1.col_names = col_namesNew
        Matrix1.f = fNew
        Matrix1.b = numpy.array([0.0]*Anew2.shape[0])
        LP1 = LinearProblem(lp_solver=self.lp_solver)
        LP1.load_matrix(Matrix1)
        self.FBA = ProblemFBA(LP1)

