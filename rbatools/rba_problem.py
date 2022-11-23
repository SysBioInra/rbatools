# python 2/3 compatibility
from __future__ import division, print_function

import numpy
import copy
import warnings
import rbatools._auxiliary_functions as _auxiliary_functions
from typing import Union
from rbatools.rba_problem_matrix import ProblemMatrix
from rbatools.rba_lp import LinearProblem
from rbatools._warnings_and_errors import *


class ProblemRBA(object):
    """
    Class holding RBA-problem as mathematical manifestation of the RBA-model.

    Attributes
    ----------
    classicRBA : boolean
        Indicates that the problem is a classic RBA-problem (as defined in RBApy)
    ClassicRBAmatrix : rba.ConstraintMatrix.matrix object
    Mu : float
        Current Growth rate as numeric value
    LP : rbatools.rba_lp.LinearProblem object
    Enzyme_FWcapacities : list
        List of constraint IDs, which represent forward efficiencies of enzymes
        Created by method 'extract_constraint_types'
    Enzyme_BWcapacities : list
        List of constraint IDs, which represent backward efficiencies of enzymes
        Created by method 'extract_constraint_types'
    ProcessCapacities : list
        List of constraint IDs, which represent net efficiencies of processes
        Created by method 'extract_constraint_types'
    Metabolites : list
        List of constraint IDs, which represent mass balances of metabolites
        Created by method 'extract_constraint_types'
    CompartmentDensities : list
        List of constraint IDs, which represent compartment capacities
        Created by method 'extract_constraint_types'
    Reactions : list
        List of constraint IDs, which represent metabolic reactions
        Created by method 'extract_variable_types'
    Enzymes : list
        List of constraint IDs, which represent enzymes, associated with metabolic reactions
        Created by method 'extract_variable_types'
    Processes : list
        List of constraint IDs, which represent process machineries
        Created by method 'extract_variable_types'
    MuDepIndices_A : list
        List of tuples holding rows (constraint IDs) and columns (variable IDs) of
        constraint-matrix coefficients (LHS), which depend on the growth rate.
        Created by method 'find_growth_rate_dependencies'
    MuDepIndices_b: list
        List of constraint IDs, whos RHS depend on the growth rate
        Created by method 'find_growth_rate_dependencies'
    MuDepIndices_LB: list
        List of variable IDs, whos lower-bounds depend on the growth rate
        Created by method 'find_growth_rate_dependencies'
    MuDepIndices_UB: list
        List of variable IDs, whos upper-bounds depend on the growth rate
        Created by method 'find_growth_rate_dependencies'
    Solved: bool
        Booelean indicating wheter Problem has been successfully solved,
        according to the acceptable solution statuses provided to the solve_lp-method.
        Created by method 'solve_lp'
    SolutionStatus: int
        Numerical value indicating the solution-status of the problem.
        Consult CPLEX documentation for their meaning.
        Created by method 'solve_lp'
    ObjectiveValue: float
        Numeric value of the objective function after optimisation by CPLEX.
        Created by method 'solve_lp'
    SolutionValues: dict
        Solution vector after optimisation by CPLEX.
        (Dictionary with variable IDs as keys and numeric values as values)
        Created by method 'solve_lp'
    DualValues: dict
        Vector of dual-values after optimisation by CPLEX.
        (Dictionary with constraint IDs as keys and numeric values as values)
        Created by method 'solve_lp'
    """
    def __init__(self,matrix,ModelStructure,lp_solver="cplex"):
        """
        Initiates rbatools.rba_problem.ProblemRBA from supplied rba.ConstraintMatrix object.

        Makes sure the problem is consistent and optimisable.

        Parameters
        ----------
        matrix : rba.ConstraintMatrix
        ModelStructure rbatools.rba_model_structure.ModelStructureRBA
        lp_solver : str
            Specifies which LP-solver should be used ('cplex' or 'swiglpk')
            Default: 'cplex'
        """

        self.classicRBA = True
        self.ClassicRBAmatrix = matrix
        self.Mu = 0
        ## Set Mu of ClassicRBAmatrix to 0 ##
        self.build_classic_matrix(Mu=self.Mu)
        ## Initiate LP-object ##
        self.LP = LinearProblem(lp_solver=lp_solver)
        ## Transfer ClassicRBAmatrix information to LP-object ##
        self.LP.load_matrix(matrix=self.ClassicRBAmatrix)
        # Build solver object
        lp_built=self.LP.build_lp() 
        if not lp_built:
            raise LinearProblemError('Unable to build LP-solver object')

        ## Solve for Mu=0 to check if problem is consistent##
        self.solve_lp()
        if not self.Solved:
            ## Problem is wrongly defined and can not be solved ##
            warnings.warn('Inconsistent RBA-problem at growth rate = 0')

        ## Growth-rate dependent indices ##
        ## and constraint/variable types are extractet ##
        self.find_growth_rate_dependencies(ModelStructure=ModelStructure)
        self.extract_constraint_types(ModelStructure=ModelStructure)
        self.extract_variable_types(ModelStructure=ModelStructure)
        self.SolutionType = ''

    def extract_constraint_types(self,ModelStructure):
        """
        Extracts information on the different constraint types in the standard RBA-matrix.

        AddsAttributes: Enzyme_FWcapacities, Enzyme_BWcapacities, ProcessCapacities,
                        Metabolites and CompartmentDensities
        """
        self.Enzyme_FWcapacities = [i["ID"] for i in ModelStructure.EnzymeConstraintsInfo.Elements.values() if i["Direction"]=="forward"]
        self.Enzyme_BWcapacities = [i["ID"] for i in ModelStructure.EnzymeConstraintsInfo.Elements.values() if i["Direction"]=="backward"]
        self.ProcessCapacities = [i["ID"] for i in ModelStructure.ProcessConstraintsInfo.Elements.values()]
        self.Metabolites = [i["AssociatedMetabolite"] for i in ModelStructure.MetaboliteConstraintsInfo.Elements.values()]
        self.CompartmentDensities = [i["ID"] for i in ModelStructure.DensityConstraintsInfo.Elements.values()]

    def extract_variable_types(self,ModelStructure):
        """
        Extracts information on the different variable types in the standard RBA-matrix.

        AddsAttributes: Enzymes, Reactions and Processes
        """
        self.Reactions = [i["ID"] for i in ModelStructure.ReactionInfo.Elements.values()]
        self.Enzymes = [i["ID"] for i in ModelStructure.EnzymeInfo.Elements.values()]
        self.Processes = [i["ID"]+"_machinery" for i in ModelStructure.ProcessInfo.Elements.values()]

    def find_growth_rate_dependencies(self,ModelStructure):
        """
        Extracts information on the growth-rate dependent LP-coefficients.

        AddsAttributes: MuDepIndices_A, MuDepIndices_b, MuDepIndices_LB
                        and MuDepIndices_UB
        """

        ## Construct 2 rba.solver.matrix objects for growth-rate 0 and 1.1 ##
        MuDepIndices_A=[]
        MuDepIndices_b=list(ModelStructure.MetaboliteConstraintsInfo.Elements.keys())
        MuDepIndices_LB=[]
        MuDepIndices_UB=[]

        M0 = copy.deepcopy(self.ClassicRBAmatrix)
        M11 = copy.deepcopy(self.ClassicRBAmatrix)
        M0.build_matrices(0)
        M11.build_matrices(1.1)
        ## Build arrays from matrices ##
        A_0 = M0.A.toarray()
        A_11 = M11.A.toarray()
        ## Find index pairs at which the two constraint-matrices differ ##
        MuDeps = numpy.where(A_11 != A_0)
        ## Transform list of numeric indices in list of tuples with row- and col names ##
        MuDepIndices_A += list(zip([self.ClassicRBAmatrix.row_names[i] for i in MuDeps[0]], [self.ClassicRBAmatrix.col_names[j] for j in MuDeps[1]]))
        MuDepIndices_A+=_auxiliary_functions.get_medium_dependent_coefficients_in_lhs(ModelStructure)
        ## Find rows at which the two righthandsides differ ##
        MuDepIndices_b += [n for n in self.ClassicRBAmatrix.row_names if M0.b[self.ClassicRBAmatrix.row_names.index(
            n)] != M11.b[self.ClassicRBAmatrix.row_names.index(n)]]
        ## Find columns at which the variable bounds differ ##
        MuDepIndices_LB += [n for n in self.ClassicRBAmatrix.col_names if M0.LB[self.ClassicRBAmatrix.col_names.index(
            n)] != M11.LB[self.ClassicRBAmatrix.col_names.index(n)]]
        MuDepIndices_UB += [n for n in self.ClassicRBAmatrix.col_names if M0.UB[self.ClassicRBAmatrix.col_names.index(
            n)] != M11.UB[self.ClassicRBAmatrix.col_names.index(n)]]
        self.MuDependencies = {'FromMatrix': {'A': list(set(MuDepIndices_A)), 'b': list(set(MuDepIndices_b)), 'LB': list(set(MuDepIndices_LB)), 'UB': list(set(MuDepIndices_UB))},
                               'FromParameters': {'A': {}, 'b': {}, 'LB': {}, 'UB': {}}}

    def build_classic_matrix(self, Mu: float):
        """
        Builds standard RBA-matrix according to growth-rate

        Parameters
        ----------
        Mu : float
            Growth rate
        """

        self.Mu = Mu
        self.ClassicRBAmatrix.build_matrices(Mu)

    def solve_lp(self, feasible_stati: list = ["optimal","feasible"], try_unscaling_if_sol_status_is_feasible_only_before_unscaling: bool =True):
        """
        Solves Linear RBA problem.

        When solution-status is amongst the user-defined feasible statuses;
        boolean 'Solved' is set to True and 'ObjectiveValue', 'SolutionValues' and 'DualValues'
        are stored as attributes.

        Parameters
        ----------
        feasible_stati : list of int
            List with identifiers of acceptable solution statuses.
            (consult ILOG-CPLEX documentation for information on them).
            Default: feasible_stati=["optimal","feasible"]
        try_unscaling_if_sol_status_is_feasible_only_before_unscaling :bool
            Try re-solving problem with less aggressive scaling, when solution status is infeasible after unscaling
            Default: True
        """

        self.Solved = False
        try:
            ## Solve cplex LP ##
            self.LP.solve_lp()
            ## Determin solution-status ##
            self.SolutionStatus = self.LP.return_solution_status()
            ## Check if solution status is amongst acceptable ones ##
            if self.SolutionStatus in feasible_stati:
                ## Extract solution-data ##
                self.Solved = True
                self.ObjectiveValue = self.LP.return_objective_value()
                self.SolutionValues = self.LP.return_primal_values()
                self.DualValues = self.LP.return_dual_values()
            else:
                if try_unscaling_if_sol_status_is_feasible_only_before_unscaling:
                    if self.SolutionStatus in ["feasible_only_before_unscaling"]:
                        self.LP.unscale_lp()
                        self.LP.solve_lp()
                        self.SolutionStatus = self.LP.return_solution_status()
                        if self.SolutionStatus in feasible_stati:
                            ## Extract solution-data ##
                            self.Solved = True
                            self.ObjectiveValue = self.LP.return_objective_value()
                            self.SolutionValues = self.LP.return_primal_values()
                            self.DualValues = self.LP.return_dual_values()
        except:
            self.ObjectiveValue = None
            self.SolutionValues = None
            self.DualValues = None
        self.SolutionType = 'Normal'

    def update_growth_rate(self, Mu: float, keepParameters: dict = {}, ModifiedProblem: bool = False):
        """
        Changes growth-rate of problem and sets all associated coefficients.

        Can be provided with 'keepParameters' argument
        to define problem coefficients which should remain unchanged.

        Parameters
        ----------
        Mu : float
            Growth rate
        keepParameters : dict
            Dictionary indicating which elements of the linear problem should
            not be affected when setting growth-rate. Possible keys of dictionary:
                'LHS': List of index-tuples (constraint-ID,variable-ID),
                       indicating elements of the lefthandside (constraint-matrix).
                'RHS': List of constraint IDs, indicating elements of the righthandside (b-vector).
                'LB': List of variable-IDs, indicating elements of the lower-bound vector.
                'UB': List of variable-IDs, indicating elements of the upper-bound vector.
            Default: keepParameters=None
        ModifiedProblem : bool
            Default: False
        """
        ## Pass new matrix and indices of elements to update to LP.update_matrix-method ##
        self.ClassicRBAmatrix.build_matrices(self.Mu)
        inputMatrix = ProblemMatrix()
        inputMatrix.load_matrix(matrix=self.ClassicRBAmatrix)

        if not list(keepParameters.keys()):
            self.LP.update_matrix(matrix=inputMatrix,
                                 Ainds=self.MuDependencies['FromMatrix']['A'],
                                 Binds=self.MuDependencies['FromMatrix']['b'],
                                 CTinds=[],
                                 LBinds=self.MuDependencies['FromMatrix']['LB'],
                                 UBinds=self.MuDependencies['FromMatrix']['UB'],
                                 ModifiedProblem=False)
        else:
            ## If indices are passed, which define elements to be not changed when setting Mu ##
            ## these indices are removed from the update-from-new-to-old indices ##
            if 'LHS' in list(keepParameters.keys()):
                A_idxs = list(set(self.MuDependencies['FromMatrix']['A'])-set(keepParameters['LHS']))
            if 'RHS' in list(keepParameters.keys()):
                B_idxs = list(set(self.MuDependencies['FromMatrix']['b'])-set(keepParameters['RHS']))
            if 'LB' in list(keepParameters.keys()):
                LB_idxs = list(set(self.MuDependencies['FromMatrix']['LB'])-set(keepParameters['LB']))
            if 'UB' in list(keepParameters.keys()):
                UB_idxs = list(set(self.MuDependencies['FromMatrix']['UB'])-set(keepParameters['UB']))
            self.LP.update_matrix(matrix=inputMatrix,
                                 Ainds=A_idxs,
                                 Binds=B_idxs,
                                 CTinds=[],
                                 LBinds=LB_idxs,
                                 UBinds=UB_idxs,
                                 ModifiedProblem=True)
        if ModifiedProblem:
            for i in list(self.MuDependencies['FromParameters']['b'].keys()):
                newPar = self.evaluate_parameter(self.MuDependencies['FromParameters']['b'][i])
                self.set_right_hand_side({i: newPar})
            for i in list(self.MuDependencies['FromParameters']['A'].keys()):
                newPar = self.evaluate_parameter(self.MuDependencies['FromParameters']['A'][i])
                self.set_problem_coefficients({i: newPar})
            for i in list(self.MuDependencies['FromParameters']['LB'].keys()):
                newPar = self.evaluate_parameter(self.MuDependencies['FromParameters']['LB'][i])
                self.set_lb({i: newPar})
            for i in list(self.MuDependencies['FromParameters']['UB'].keys()):
                newPar = self.evaluate_parameter(self.MuDependencies['FromParameters']['UB'][i])
                self.set_ub({i: newPar})

    def set_growth_rate(self, Mu: float, ModelStructure, keepParameters: dict = {}):
        """
        Changes growth-rate of problem and sets all associated coefficients.

        Can be provided with 'keepParameters' argument
        to define problem coefficients which should remain unchanged.

        Parameters
        ----------
        Mu : float
            Growth rate
        ModelStructure : RBA_ModellStructure object.
        keepParameters : dict
            Dictionary indicating which elements of the linear problem should
            not be affected when setting growth-rate. Possible keys of dictionary:
                'LHS': List of index-tuples (constraint-ID,variable-ID),
                       indicating elements of the lefthandside (constraint-matrix).
                'RHS': List of constraint IDs, indicating elements of the righthandside (b-vector).
                'LB': List of variable-IDs, indicating elements of the lower-bound vector.
                'UB': List of variable-IDs, indicating elements of the upper-bound vector.
            Default: keepParameters=None
        """

        self.Mu = float(Mu)
        NumberParDeps=len(ModelStructure.MuDependencies)
        if self.LP.row_names == self.ClassicRBAmatrix.row_names and self.LP.col_names == self.ClassicRBAmatrix.col_names and self.classicRBA and NumberParDeps == 0:
            self.update_growth_rate(Mu=float(Mu), keepParameters=keepParameters, ModifiedProblem=False)
        else:
            self.update_growth_rate(Mu=float(Mu), keepParameters=keepParameters, ModifiedProblem=True)

    def get_constraint_types(self, constraints: Union[list,str] = []) -> dict:
        """
        Extracts type of constraints.

        Parameters
        ----------
        constraints : str or list of str
            Constraints to retreive the type for.
            Either constraint ID or list of constraint IDs to specify the type
            of which constraint to look up.
            This is an optional input; if not provided all constraint are looked up.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            type identification-characters as values.
        """
        if type(constraints) is str:
            names=[constraints]
        elif type(constraints) is list:
            if len(constraints)>0:
                names = constraints
            else:
                names = self.LP.row_names
        return(self.LP.get_constraint_types(constraints=names))

    def set_constraint_types(self, inputDict: dict):
        """
        Sets type of constraints.

        E: = ; L: <= ; G: >=

        Parameters
        ----------
        inputDict : dict
            Dictionary with constraint-IDs as keys and type identification-character as values.
            ({'col1':'E','col2':'L', ...}).
        """
        self.LP.set_constraint_types(inputDict=inputDict)

    def get_objective(self, variables: Union[list,str] = []) -> dict:
        """
        Returns objective coefficient of variables in problem.

        Parameters
        ----------
        variables : str or list of str
            Variables to retreive the objective coefficients for.
            Either variable ID or list of variable IDs to specify
            the coefficients of which variables to look up.
            This is an optional input; if not provided all variables are looked up.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            objective coefficients as values.
        """

        if type(variables) is str:
            vrs=[variables]
        elif type(variables) is list:
            if len(variables)>0:
                vrs = variables
            else:
                vrs = self.LP.col_names
        return(self.LP.get_objective(variables=vrs))

    def set_objective(self, inputDict: dict):
        """
        Sets objective function coefficients.

        Parameters
        ----------
        inputDict : dict
            Dictionary with variable-IDs as keys and new numeric values as values.
            ({'col1':42,'col2':9000, ...}).
        """
        ##Update in cplex.Cplex LP##
        self.LP.set_objective(inputDict={i:inputDict[i] for i in inputDict.keys() if numpy.isfinite(inputDict[i])})

    def clear_objective(self):
        """
        Sets all coefficients of the objective function to zero.
        """
        self.LP.set_objective(inputDict=dict(zip(self.LP.col_names,[0.0]*len(self.LP.col_names))))

    def invert_objective(self):
        """
        Changes sign (optimisation-sense) of objective function.

        """
        current_objective = self.LP.get_objective(variables=self.LP.col_names)
        self.LP.set_objective(inputDict={v:-current_objective[v] for v in current_objective.keys() if current_objective[v]!=0})

    def get_right_hand_side(self, constraints:Union[list,str] = []) -> dict:
        """
        Extracts coefficients of problem's righthand side (B-vector).

        Parameters
        ----------
        constraints : str or list of str
            Constraints to retreive the objective coefficients for.
            Either constraint ID or list of constraint IDs to specify the RHS
            of which constraint to look up.
            This is an optional input; if not provided all constraint are looked up.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            RHS-values as values.
        """

        if type(constraints) is str:
            names=[constraints]
        elif type(constraints) is list:
            if len(constraints)>0:
                names = constraints
            else:
                names = self.LP.row_names
        return(self.LP.get_right_hand_side(constraints=names))

    def set_right_hand_side(self, inputDict: dict):
        """
        Set coefficients of the problems' RHS (b-vector).
        Parameters
        ----------
        inputDict : dict
            Dictionary with constraint-IDs as keys and new numeric values as values.
            ({'row1':42,'row2':9000, ...}).
        """
        self.LP.set_right_hand_side(inputDict={i:inputDict[i] for i in inputDict.keys() if numpy.isfinite(inputDict[i])})

    def calculate_left_hand_side(self, constraints:Union[list,str] = []) -> dict:
        """
        Calculates value of problem's lefthand side
        (after multiplying with solution-vector).

        Parameters
        ----------
        constraints : str or list of str
            Constraints to retreive the LHS-value for.
            Either constraint ID or list of constraint IDs to specify the LHS-value
            of which constraint to look up.
            This is an optional input; if not provided all constraint are looked up.
        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            LHS-value as values.
        """
        if type(constraints) is str:
            names=[constraints]
        elif type(constraints) is list:
            if len(constraints)>0:
                names = constraints
            else:
                names = self.LP.row_names
        Sol=numpy.array([self.SolutionValues[i] for i in self.LP.col_names])
        Amat = self.LP.A.toarray()
        multRes = Amat.dot(Sol)
        out = {c: multRes[self.LP.row_names.index(c)] for c in names}
        return(out)

    def get_problem_coefficients(self, inputTuples:Union[list,tuple] = []) -> dict:
        """
        Returns coefficients of LHS of problem.

        Parameters
        ----------
        inputTuples : tuple or list of tuples.
            Tuples hold row and column indices.
            [('row1','col1'),('row2','col2'),...] or ('row1','col1').

        Returns
        ----------
        dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """

        if type(inputTuples) is list:
            tuples = inputTuples
        elif type(inputTuples) is tuple:
            tuples = [inputTuples]
        else:
            raise InputError('Input-argument not tuple or list of tuples')
        return(self.LP.get_problem_coefficients(inputTuples=tuples))

    def set_problem_coefficients(self, inputDict: dict):
        """
        Set coefficients of the problems' LHS (constraint matrix).

        Parameters
        ----------
        inputDict : dict
            Dict with index-tuples ('row1','col1') as keys
            and new numeric values as values.
            ({('row1','col1'):42,('row2','col2'):9000, ...}).
        """

        variables = []
        constraints = []
        coefficients = []
        Changes = []
        for const in list(inputDict.keys()):
            for var in list(inputDict[const].keys()):
                constraints.append(const)
                variables.append(var)
                coefficients.append(numpy.float64(inputDict[const][var]))
        self.LP.set_problem_coefficients(input=list(zip(constraints, variables, coefficients)))

    def get_ub(self, variables: Union[list,str] = []) -> dict:
        """
        Returns upper bounds of problem variables.

        Parameters
        ----------
        variables : str list of str
            Variables to retreive the objective coefficients for.
            Either variable ID or list of variable IDs to specify
            the coefficients of which variables to look up.
            This is an optional input; if not provided all variables are looked up.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            upper bounds as values.
        """

        if type(variables) is str:
            names=[variables]
        elif type(variables) is list:
            if len(variables)>0:
                names = variables
            else:
                names = self.LP.col_names
        return(self.LP.get_ub(variables=names))

    def get_lb(self, variables: Union[list,str] = []) -> dict:
        """
        Returns lower bounds of problem variables.

        Parameters
        ----------
        variables : str list of str
            Variables to retreive the objective coefficients for.
            Either variable ID or list of variable IDs to specify
            the coefficients of which variables to look up.
            This is an optional input; if not provided all variables are looked up.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            lower bounds as values.
        """

        if type(variables) is str:
            names=[variables]
        elif type(variables) is list:
            if len(variables)>0:
                names = variables
            else:
                names = self.LP.col_names
        return(self.LP.get_lb(variables=names))

    def set_lb(self, inputDict: dict):
        """
        Set lower-bounds of the problem variables.

        Parameters
        ----------
        inputDict : dict
            Dictionary with variable-IDs as keys and new numeric values as values.
            ({'col1':42,'col2':9000, ...}).
        """

        self.LP.set_lb(inputDict={i:inputDict[i] for i in inputDict.keys() if numpy.isfinite(inputDict[i])})

    def set_ub(self, inputDict: dict):
        """
        Set upper-bounds of the problem variables.

        Parameters
        ----------
        inputDict : dict
            Dictionary with variable-IDs as keys and new numeric values as values.
            ({'col1':42,'col2':9000, ...}).
        """
        ##Update in cplex.Cplex LP##
        self.LP.set_ub(inputDict={i:inputDict[i] for i in inputDict.keys() if numpy.isfinite(inputDict[i])})

    def evaluate_parameter(self, definition):
        if type(definition) == str:
            return(self.ClassicRBAmatrix._blocks.parameters.__getitem__(definition).value)
        elif type(definition) == dict:
            variables = {i: self.ClassicRBAmatrix._blocks.parameters.__getitem__(
                i).value for i in definition['Variables']}
            return(eval(str(definition['Equation']), variables))
