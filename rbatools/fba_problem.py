from __future__ import division, print_function
import numpy
import scipy
import copy
import json
import warnings
from typing import Union
from rbatools.rba_problem_matrix import ProblemMatrix
from rbatools._warnings_and_errors import *

class ProblemFBA(object):
    """
    Class holding RBA-problem as mathematical manifestation of the RBA-model.

    Attributes
    ----------
    LP : rbatools.rba_lp.LinearProblem object
    parsimonious: bool
        Indicator whether FBA problem is parsimonious
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
    def __init__(self, FBA):
        self.LP = FBA
        self.LP.build_lp()
        self.parsimonious = False

    def parsimonise(self,rxns_to_ignore_in_objective=[]):
        """
        Transforms current FBA-problem into a parsimonious FBA problem.
        Reversible reactions are duplicated into a (irreversible) forward and backward-variants.
        Objective function is defined to minimize total sum of (internal reactions).

        Parameters
        ----------
        rxns_to_ignore_in_objective : list of str
        List of reactions to ignore in objective function 
        (not subject to principle of parsimony)
         """

        A = self.LP.A.toarray()
        reversibleRxns = [i for i in self.LP.col_names if not i.startswith(
            'R_EX_') and self.LP.LB[self.LP.col_names.index(i)] < 0]
        Anew = numpy.zeros((len(self.LP.row_names), len(reversibleRxns)))
        fNew = numpy.zeros(len(reversibleRxns))
        LBnew = numpy.zeros(len(reversibleRxns))
        UBnew = numpy.zeros(len(reversibleRxns))
        addIndexCount = 0
        newCols = []
        for i in reversibleRxns:
            colIndex = self.LP.col_names.index(i)
            LBold = self.LP.LB[colIndex]
            rxnCol = A[:, colIndex]
            Anew[:, addIndexCount] = -rxnCol
            fNew[addIndexCount] = 1
            LBnew[addIndexCount] = 0
            UBnew[addIndexCount] = -LBold
            newCols.append(str('Backward_'+i))
            self.LP.LB[colIndex] = 0
            addIndexCount += 1

        ReverseMatrix = ProblemMatrix()
        ReverseMatrix.A = scipy.sparse.coo_matrix(Anew)
        ReverseMatrix.b = self.LP.b
        ReverseMatrix.f = fNew
        ReverseMatrix.LB = LBnew
        ReverseMatrix.UB = UBnew
        ReverseMatrix.row_signs = self.LP.row_signs
        ReverseMatrix.row_names = self.LP.row_names
        ReverseMatrix.col_names = newCols
        ReverseMatrix.map_indices()
        self.LP.add_matrix(matrix=ReverseMatrix)
        internalRxns = [i for i in self.LP.col_names if not i.startswith(
            'R_EX_') and not i.startswith('R_maintenance_') and not i.startswith('R_BIOMASS_')]
        for i in range(len(self.LP.col_names)):
            if self.LP.col_names[i] in internalRxns:
                if self.LP.col_names[i] not in rxns_to_ignore_in_objective:
                    self.LP.f[i] = 1
                else:
                    self.LP.f[i] = 0
            else:
                self.LP.f[i] = 0
        self.LP.build_lp()
        self.parsimonious = True

    def solve_lp(self, feasible_stati: list = ["optimal","feasible"], try_unscaling_if_sol_status_is_feasible_only_before_unscaling: bool =True):
        """
        Solves Linear fBA problem.

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
                            if self.parsimonious:
                                solution_values=self.LP.return_primal_values()
                                self.SolutionValues={}
                                for rxn in solution_values.keys():
                                    if rxn.startswith('Backward_'):
                                        self.SolutionValues[rxn.split('Backward_')[1]]=-solution_values[rxn]
                                    else:
                                        self.SolutionValues[rxn]=solution_values[rxn]
                            else:
                                self.SolutionValues = self.LP.return_primal_values()
                            self.DualValues = self.LP.return_dual_values()
        except:
            self.ObjectiveValue = None
            self.SolutionValues = None
            self.DualValues = None
        self.SolutionType = 'Normal'

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
            #warnings.warn('Input-argument not tuple or list of tuples')
            return
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

    def export_escher_map(self, Filename):
        with open(Filename, 'w') as fout:
            fout.write(json.dumps({i[2:]: self.SolutionValues[i]
                                   for i in self.SolutionValues.keys()}, indent=4))
