# python 2/3 compatibility
from __future__ import division, print_function

import abc
import numpy
import copy
import itertools
import rbatools._auxiliary_functions as _auxiliary_functions
import warnings
from scipy.sparse import SparseEfficiencyWarning , csc_matrix , coo_matrix , lil_matrix , csr_matrix
from rbatools.rba_problem_matrix import ProblemMatrix
from rbatools._warnings_and_errors import *

try:
    from cplex import SparsePair
    from cplex import Cplex
    cplex_available = True
except ModuleNotFoundError:
    cplex_available = False
try:
    from swiglpk import *
    swiglpk_available = True
except ModuleNotFoundError:
    swiglpk_available = False

SOLVER_MAP={"cplex":{},"swiglpk":{}}
if cplex_available:
    SOLVER_MAP.update({"cplex":{"feasible":[],"feasible_only_before_unscaling":[5],"optimal":[1]}})
if swiglpk_available:
    SOLVER_MAP.update({"swiglpk":{"feasible":[swiglpk.GLP_FEAS],"feasible_only_before_unscaling":[],"optimal":[swiglpk.GLP_OPT]}})

warnings.simplefilter('ignore', SparseEfficiencyWarning)


class LinearProblem(ProblemMatrix):
    """

    Attributes
    ----------
    A : scipy.sparse.coo_matrix
        Lefthandside of constraint Matrix (aka Constraint Matrix)
    b : numpy.array
        Righthandside of Constraint Matrix
    row_signs : list
        Type of constraints ('E' for equality 'L' for lower-or-equal inequality) --> Ax=b or Ax<=b
    f : numyp.array
        Objective function linear get_coefficients
    LB : numpy.array
        Lower bounds of decision-variables
    UB : numpy.array
        Upper bounds of decision-variables
    row_names : list
        Names of constraints
    col_names : list
        Names of decision-variables
    rowIndicesMap : dict
        Dictionary mapping constraint names to their numeric index (generated automatically)
    colIndicesMap : dict
        Dictionary mapping variable names to their numeric index (generated automatically)
    lp_solver : str
        Selected linear problem solver.
    _lp_solver : rbatools.rba_lp._Solver (any descendent class) object
    """

    def __init__(self, *matrix:ProblemMatrix,lp_solver="cplex"):
        """
        Initiates LP object.
        If provided with an rbatools.rba_problem_matrix.ProblemMatrix object as argument, assumes its attributes.
        If not provided with an ProblemMatrix object, initiates as empty LP object.
        Parameters
        ----------
        matrix : rbatools.rba_problem_matrix.ProblemMatrix
            Optional
        lp_solver : str
            Specifies which LP-solver should be used ('cplex' or 'swiglpk')
            Default: 'cplex'
        """
        if lp_solver in SOLVER_MAP:
            self.lp_solver=lp_solver
            self.initiate_lp_solver() 
        else:
            raise InputError(str('Solver {} not supported. Available options: {}'.format(lp_solver," , ".join(list(SOLVER_MAP.keys())))))

        if isinstance(matrix, ProblemMatrix):
            self.load_matrix(matrix=matrix)
        else:
            ProblemMatrix.__init__(self)
        self.map_indices()

    def initiate_lp_solver(self):
        if self.lp_solver=="cplex":
            if cplex_available:
                self._lp_solver=_SolverCPLEX()
            else:
                raise DependencyError(str("Solver {} not found, try another one".format(self.lp_solver)))
        elif self.lp_solver=="swiglpk":
            if swiglpk_available:
                self._lp_solver=_SolverGLPK()
            else:
                raise DependencyError(str("Solver {} not found, make sure package swiglpk is working".format(self.lp_solver)))

    def update_matrix(self, matrix, Ainds=None, Binds=None, CTinds=None, LBinds=None, UBinds=None, ModifiedProblem=True):
        """
        Overwrites coefficients with new values from argument 'matrix'.

        Parameters
        ----------
        matrix : rbatools.rba_problem_matrix.ProblemMatrix
            The matrix with elements to be added
        Ainds : list of tuples
            List of index-pair tuples [(row_coeff1,col_coeff1),(row_coeff2,col_coeff2),...],
            specifying which elements are updated.
            Default is None
            (then all elements of argument matrix (present in both matrices) are taken to update)
        Binds : list
            List of constraint-IDs for indices of the RHS, specifying which RHS Values are updated.
            Default is None (then all constraints (present in both matrices) are taken to update)
        CTinds : list
            List of constraint-IDs, specifying which row_signs are updated.
            Default is None (then all constraint types (present in both matrices) are taken to update)
        LBinds : list
            List of variable-IDs, specifying which lower-bounds values are updated.
            Default is None (then all variables (present in both matrices) are taken to update)
        UBinds : list
            List of variable-IDs, specifying which upper-bounds values are updated.
            Default is None (then all variables (present in both matrices) are taken to update)
        ModifiedProblem : boolean
            Whether problem has been modifed. If false own lp components are fully replaced by input, if rows and  col names are identical .
        """
        ## Check if indices to update are provided, if not use all indices ##
        if Ainds is None:
            Ainds = list(itertools.product(matrix.row_names, matrix.col_names))
        if Binds is None:
            Binds = matrix.row_names
        if CTinds is None:
            CTinds = matrix.row_names
        if LBinds is None:
            LBinds = matrix.col_names
        if UBinds is None:
            UBinds = matrix.col_names

        ####### Update constraint-matrix LHS (A) #######
        if matrix.row_names == self.row_names and matrix.col_names == self.col_names and not ModifiedProblem:
            self.A=matrix.A
            self.b=matrix.b
            self.LB=matrix.LB
            self.UB=matrix.UB
            self.row_signs=matrix.row_signs
        else:
            ### If old and new matrix do not have same elements and order of indices ###
            ## Find elements (index pairs) which are in the old, as well in th new matrix. ##
            x1, x2 = zip(*Ainds)
            rowsOld = tuple([self.rowIndicesMap[i] for i in x1])
            colsOld = tuple([self.colIndicesMap[i] for i in x2])
            rowsNew = tuple([matrix.rowIndicesMap[i] for i in x1])
            colsNew = tuple([matrix.colIndicesMap[i] for i in x2])

            newA = csc_matrix(matrix.A)
            oldA = csc_matrix(self.A)
            #Overwrite old elements at indices with corresponding elements from new matrix#
            oldA[tuple(rowsOld), tuple(colsOld)] = newA[tuple(rowsNew), tuple(colsNew)]
            self.A = coo_matrix(oldA)

            ## Update RHS (b)##
            if len(Binds) > 0:
                if matrix.row_names == self.row_names:
                    ## If old and new LPs have same rows and row-order ##
                    #Find numeric indices of rows to update (same for old and new matrix)#
                    rowsNew = [matrix.rowIndicesMap[i] for i in Binds]
                    #Overwrite old elements at row-indices with corresponding new elements#
                    for bind in list(range(len(rowsNew))):
                        self.b[rowsNew[bind]] = matrix.b[rowsNew[bind]]
                else:
                    x = [i for i in Binds if i in matrix.row_names and i in self.row_names]
                    #Find numeric indices of rows to update (for old and new matrix individually)#
                    rowsNew = [matrix.rowIndicesMap[i] for i in x]
                    rowsOld = [self.rowIndicesMap[i] for i in x]
                    #Overwrite old elements at row-indices with corresponding new elements#
                    for bind in list(range(len(rowsNew))):
                        self.b[rowsOld[bind]] = matrix.b[rowsNew[bind]]

            ## Update Constraint type ##
            if len(CTinds) > 0:
                if matrix.row_names == self.row_names:
                    rowsNew = [matrix.rowIndicesMap[i] for i in CTinds]
                    self.row_signs = matrix.row_signs
                else:
                    for row in CTinds:
                        self.row_signs[self.rowIndicesMap[row]] = matrix.row_signs[matrix.rowIndicesMap[row]]

            ## Update LB##
            if len(LBinds) > 0:
                oLB = numpy.array(self.LB)
                nLB = numpy.array(matrix.LB)
                if matrix.col_names == self.col_names:
                    colsNew = [matrix.colIndicesMap[i] for i in LBinds]
                    oLB[colsNew] = nLB[colsNew]
                else:
                    x = [i for i in LBinds if i in matrix.col_names and i in self.col_names]
                    colsOld = [self.colIndicesMap[i] for i in x]
                    colsNew = [matrix.colIndicesMap[i] for i in x]
                    oLB[colsOld] = nLB[colsNew]
                self.LB = oLB

            ## Update UB##
            if len(UBinds) > 0:
                oUB = numpy.array(self.UB)
                nUB = numpy.array(matrix.UB)
                if matrix.col_names == self.col_names:
                    colsNew = [matrix.colIndicesMap[i] for i in UBinds]
                    oUB[colsNew] = nUB[colsNew]
                else:
                    x = [i for i in UBinds if i in matrix.col_names and i in self.col_names]
                    colsOld = [self.colIndicesMap[i] for i in x]
                    colsNew = [matrix.colIndicesMap[i] for i in x]
                    oUB[colsOld] = nUB[colsNew]
                self.UB = oUB

        self._lp_solver.import_lp(input_lp=self)
        self.build_lp()

    def add_matrix(self, matrix:ProblemMatrix):
        """
        Merges the Problem with the one provided as input-argument to this method.

        Matrix elements unique to input-matrix are added.
        Elements occuring in both matrices are overwritten with the value from new matrix.
        Generates CPLEX problem from merged problem.

        Parameters
        ----------
        matrix : rbatools.rba_problem_matrix.ProblemMatrix
            The matrix with elements to be added
        """
        matrix.map_indices()
        oldA = copy.deepcopy(self.A.toarray())
        if type(matrix.A) is numpy.ndarray:
            matrix.A = matrix.A.astype('float64')
        else:
            matrix.A = matrix.A.toarray().astype('float64')
        ## Determine union of old- and new matrix's row-names.##
        ## Make sure the new names, not present in the old matrix are added to the end.##
        ## Same thing also with column-names ##
        ## Initiate compound RBA matrix and adapt elements to the compound LP with the new dimensions ##
        compoundProblem = ProblemMatrix()
        compoundProblem.row_names = list(self.row_names + list(set(matrix.row_names)-set(self.row_names)))
        compoundProblem.col_names = list(self.col_names + list(set(matrix.col_names)-set(self.col_names)))
        compoundProblem.A = numpy.zeros((len(compoundProblem.row_names), len(compoundProblem.col_names)))
        compoundProblem.a_to_lil()
        compoundProblem.b = numpy.zeros(len(compoundProblem.row_names))
        compoundProblem.row_signs = ['E']*len(compoundProblem.row_names)
        compoundProblem.f = numpy.zeros(len(compoundProblem.col_names))
        compoundProblem.LB = numpy.zeros(len(compoundProblem.col_names))
        compoundProblem.UB = numpy.zeros(len(compoundProblem.col_names))
        compoundProblem.map_indices()

        ## Since it has been made sure that the indices present in the original problem, ##
        ## are present in the compound problem at the beginning in the exact same order; ##
        ## One can now just put the information of the original problem at the beginning ##
        compoundProblem.A[0:oldA.shape[0], 0:oldA.shape[1]] = oldA
        compoundProblem.b[0:oldA.shape[0]] = copy.deepcopy(self.b)
        compoundProblem.f[0:oldA.shape[1]] = copy.deepcopy(self.f)
        compoundProblem.LB[0:oldA.shape[1]] = copy.deepcopy(self.LB)
        compoundProblem.UB[0:oldA.shape[1]] = copy.deepcopy(self.UB)
        compoundProblem.row_signs[0:oldA.shape[0]] = copy.deepcopy(self.row_signs)
        for i in matrix.row_names:
            for j in matrix.col_names:
                NewMatrixRowIndex = matrix.rowIndicesMap[i]
                NewMatrixColIndex = matrix.colIndicesMap[j]
                CompoundMatrixRowIndex = compoundProblem.rowIndicesMap[i]
                CompoundMatrixColIndex = compoundProblem.colIndicesMap[j]
                compoundProblem.A[CompoundMatrixRowIndex,
                                  CompoundMatrixColIndex] = matrix.A[NewMatrixRowIndex, NewMatrixColIndex]

        ## Find numeric indices of new-matrix elements in the new matrix ##
        NewMatrixRowIndices = tuple([matrix.rowIndicesMap[i] for i in matrix.row_names])
        NewMatrixColIndices = tuple([matrix.colIndicesMap[i] for i in matrix.col_names])

        ## Find numeric indices of new-matrix elements in the compound matrix ##
        CompoundMatrixRowIndices = tuple([compoundProblem.rowIndicesMap[i]
                                          for i in matrix.row_names])
        CompoundMatrixColIndices = tuple([compoundProblem.colIndicesMap[i]
                                          for i in matrix.col_names])

        ## Transfer new-matrix elements to compound problem ##
        compoundProblem.f[list(CompoundMatrixColIndices)] = matrix.f[list(NewMatrixColIndices)]
        compoundProblem.LB[list(CompoundMatrixColIndices)] = matrix.LB[list(NewMatrixColIndices)]
        compoundProblem.UB[list(CompoundMatrixColIndices)] = matrix.UB[list(NewMatrixColIndices)]
        compoundProblem.b[list(CompoundMatrixRowIndices)] = matrix.b[list(NewMatrixRowIndices)]

        for i in range(len(NewMatrixRowIndices)):
            compoundProblem.row_signs[CompoundMatrixRowIndices[i]] = matrix.row_signs[NewMatrixRowIndices[i]]

        ## Overwrite old matrix with compound problem. And rebuild CPLEX-LP if there was one before##
        self.load_matrix(compoundProblem)
        self._lp_solver.import_lp(input_lp=self)
        self.build_lp()

    def load_matrix(self, matrix):
        """
        Imports information from compatible object.

        Imports information for attributes from input-objects (named matrix)
        with the respectively named fields, required.
        Required fields are:
        'A','b','row_signs','f','LB','UB','row_names' and 'col_names'.

        Parameters
        ----------
        matrix : any LP-object holding the required fields
            The matrix with elements to be added
        """

        if _auxiliary_functions.check_for_attributes(matrix, ['A', 'b', 'f', 'LB', 'UB', 'row_signs', 'row_names', 'col_names']):
            if type(matrix.A) is coo_matrix:
                self.A = matrix.A.astype('float64')
            elif type(matrix.A) is numpy.ndarray:
                self.A = coo_matrix(matrix.A.astype('float64'))
            elif type(matrix.A) is lil_matrix:
                self.A = coo_matrix(matrix.A.astype('float64'))
            elif type(matrix.A) is csc_matrix:
                self.A = coo_matrix(matrix.A.astype('float64'))
            elif type(matrix.A) is csr_matrix:
                self.A = coo_matrix(matrix.A.astype('float64'))
            else:
                raise InputError('Constraint matrix A must be an array or sparse matrix')
                #warnings.warn('Constraint matrix A must be an array or sparse matrix')
                #return
            if type(matrix.b) is numpy.ndarray:
                self.b = matrix.b.astype('float64')
            else:
                raise InputError('Righthand side vector b must be of type numpy array')
                #warnings.warn('Righthand side vector b must be of type numpy array')
                #return
            if type(matrix.f) is numpy.ndarray:
                self.f = matrix.f.astype('float64')
            else:
                raise InputError('Objective function vector f must be of type numpy array')
                #warnings.warn('Objective function vector f must be of type numpy array')
                #return
            if type(matrix.LB) is numpy.ndarray:
                self.LB = matrix.LB.astype('float64')
            else:
                raise InputError('Lower bound vector LB must be of type numpy array')
                #warnings.warn('Lower bound vector LB must be of type numpy array')
                #return
            if type(matrix.UB) is numpy.ndarray:
                self.UB = matrix.UB.astype('float64')
            else:
                raise InputError('Upper bound vector UB must be of type numpy array')
                #warnings.warn('Upper bound vector UB must be of type numpy array')
                #return
            if type(matrix.row_signs) is list:
                self.row_signs = matrix.row_signs
            else:
                raise InputError('Row signs list must be of type list')
                #warnings.warn('Row signs list must be of type list')
                #return
            if type(matrix.row_names) is list:
                self.row_names = matrix.row_names
            else:
                raise InputError('Row names list must be of type list')
                #warnings.warn('Row names list must be of type list')
                #return
            if type(matrix.col_names) is list:
                self.col_names = matrix.col_names
            else:
                raise InputError('Column names list must be of type list')
                #warnings.warn('Column names list must be of type list')
                #return
        else:
            raise InputError('Input does not have all necessary elements')
            #warnings.warn('Input does not have all necessary elements')
            #return
        self.map_indices()
        self._lp_solver.import_lp(input_lp=self)

    def build_lp(self) -> bool:
        """
        Constructs an  LP-object, based on solver attribute

        Returns
        ----------
        bool : Boolean, indicating whether LP-object could be constructed
        """
        self._lp_solver.build_lp()
        return(True)

    def solve_lp(self):
        """
        solves linear problem
        """
        self._lp_solver.solve_lp()

    def return_solution_status(self) -> dict:
        """
        Returns solution status, after solving the problem.

        Returns
        ----------
        str : solution status
            optimal, feasible or infeasible
            solution stati returned by solvers:
                optimal:
                    swiglpk -> GLP_OPT
                    cplex -> 1
                feasible:
                    swiglpk -> GLP_FEAS
                feasible_only_before_unscaling:
                    cplex ->  5
        """
        return(self._lp_solver.return_solution_status())

    def return_primal_values(self) -> dict:
        """
        Returns (primal)solution vector.

        Returns
        ----------
        dict : Problem variables as keys and solution values as values.
        """
        return(self._lp_solver.return_primal_values())

    def return_dual_values(self) -> dict:
        """
        Returns (dual) solution vector.

        Returns
        ----------
        dict : Problem constraints as keys and dual values (shadow prices) as values.
        """
        return(self._lp_solver.return_dual_values())

    def return_objective_value(self) -> float:
        """
        Returns objective value, after solving the problem.

        Returns
        ----------
        float : objective value
        """
        return(self._lp_solver.return_objective_value())

    def unscale_lp(self):
        """
        Unscale LP
        """
        self._lp_solver.unscale_lp()

    def get_constraint_types(self, constraints: list = []) -> dict:
        """
        Returns constraint types for specified rows.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        constraints : list of row names.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            type identification-characters as values.
        """
        return(self._lp_solver.get_constraint_types(constraints=constraints))

    def set_constraint_types(self, inputDict: dict):
        """
        Sets the type of specified constraints.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and constraint type as values
        """
        self._lp_solver.set_constraint_types(inputDict=inputDict)

    def get_objective(self, variables: list = []) -> dict:
        """
        Returns objective coefficients for specified columns.

        Parameters
        ----------
        variables : list of column names.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            objective coefficients as values.
        """
        return({col:self.f[self.col_names.index(col)] for col in variables})

    def set_objective(self, inputDict: dict):
        """
        Sets coefficients in objective function to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and objective coefficient as values
        """
        self._lp_solver.set_objective(inputDict=inputDict)

    def get_right_hand_side(self, constraints:list = []) -> dict:
        """
        Returns RHS values for specified rows.

        Parameters
        ----------
        constraints : list of row names.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            RHS-values as values.
        """
        return({row:self.b[self.row_names.index(row)] for row in constraints})

    def set_right_hand_side(self,inputDict):
        """
        Sets the rhs-value of specified constraints.

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and rhs-values as values.
        """
        self._lp_solver.set_right_hand_side(inputDict=inputDict)

    def get_problem_coefficients(self, inputTuples:list = []) -> dict:
        """
        Returns problem coefficints of constraint matrix.

        Parameters
        ----------
        inputTuples : list of coordinate tuples of coefficients (row,col).

        Returns
        ----------
        dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """
        return({tuple:self.A.toarray()[tuple[0],tuple[1]] for tuple in inputTuples})

    def set_problem_coefficients(self, input: list):
        """
        Sets coefficients in constraint matrix.

        Parameters
        ----------
        inputDict : dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """
        self._lp_solver.set_problem_coefficients(input=input)

    def get_ub(self,variables):
        """
        Returns upper bounds for specified columns.

        Parameters
        ----------
        variables : list of column names.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            upper bounds as values.
        """
        return({v:self.UB[self.col_names.index(v)] for v in variables})

    def get_lb(self,variables):
        """
        Returns lower bounds for specified columns.

        Parameters
        ----------
        variables : list of column names.

        Returns
        ----------
        dict
            Dictionary with variable-IDs as keys and
            lower bounds as values.
        """
        return({v:self.LB[self.col_names.index(v)] for v in variables})

    def set_lb(self, inputDict: dict):
        """
        Sets variable lower bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and lower bounds as values
        """
        self._lp_solver.set_lb(inputDict=inputDict)

    def set_ub(self, inputDict: dict):
        """
        Sets variable upper bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and upper bounds as values
        """
        self._lp_solver.set_ub(inputDict=inputDict)

class _Solver(abc.ABC):

    def __init__(self):
        self.A = coo_matrix(numpy.array([]))
        self.b = numpy.array([])
        self.f = numpy.array([])
        self.LB = numpy.array([])
        self.UB = numpy.array([])
        self.row_signs = []
        self.row_names = []
        self.col_names = []

    @property
    @abc.abstractmethod
    def name(self):
        pass

    def import_lp(self,input_lp):
        self.A=input_lp.A
        self.b=input_lp.b
        self.row_signs=input_lp.row_signs
        self.LB=input_lp.LB
        self.UB=input_lp.UB
        self.f=input_lp.f
        self.row_names=list(input_lp.row_names)
        self.col_names=list(input_lp.col_names)


    @abc.abstractmethod
    def build_lp(self):
        pass

    @abc.abstractmethod
    def solve_lp(self):
        """
        solves linear problem
        """
        pass

    @abc.abstractmethod
    def return_solution_status(self) -> dict:
        """
        Returns solution status, after solving the problem.

        Returns
        ----------
        str : solution status
            optimal, feasible or infeasible
            solution stati returned by solvers:
                optimal:
                    swiglpk -> GLP_OPT
                    cplex -> 1
                feasible:
                    swiglpk -> GLP_FEAS
                feasible_only_before_unscaling:
                    cplex ->  5
        """
        pass

    @abc.abstractmethod
    def return_primal_values(self) -> dict:
        """
        Returns (primal)solution vector.

        Returns
        ----------
        dict : Problem variables as keys and solution values as values.
        """
        pass

    @abc.abstractmethod
    def return_dual_values(self) -> dict:
        """
        Returns (dual) solution vector.

        Returns
        ----------
        dict : Problem constraints as keys and dual values (shadow prices) as values.
        """
        pass

    @abc.abstractmethod
    def return_objective_value(self) -> float:
        """
        Returns objective value, after solving the problem.

        Returns
        ----------
        float : objective value
        """
        pass

    @abc.abstractmethod
    def unscale_lp():
        """
        Unscale LP
        """
        pass

    @abc.abstractmethod
    def get_constraint_types(self, constraints: list = []) -> dict:
        """
        Returns constraint types for specified rows.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        constraints : list of row names.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            type identification-characters as values.
        """
        pass

    @abc.abstractmethod
    def set_constraint_types(self, inputDict: dict):
        """
        Sets the type of specified constraints.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and constraint type as values
        """
        pass

    @abc.abstractmethod
    def set_objective(self, inputDict: dict):
        """
        Sets coefficients in objective function to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and objective coefficient as values
        """
        pass

    @abc.abstractmethod
    def set_right_hand_side(self,inputDict):
        """
        Sets the rhs-value of specified constraints.

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and rhs-values as values.
        """
        pass

    @abc.abstractmethod
    def set_problem_coefficients(self, input: list):
        """
        Sets coefficients in constraint matrix.

        Parameters
        ----------
        inputDict : dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """
        pass

    @abc.abstractmethod
    def set_lb(self, inputDict: dict):
        """
        Sets variable lower bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and lower bounds as values
        """
        pass

    @abc.abstractmethod
    def set_ub(self, inputDict: dict):
        """
        Sets variable upper bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and upper bounds as values
        """
        pass

class _SolverGLPK(_Solver):

    @property
    def name(self):
        return("swiglpk")

    def build_lp(self):
        """
        Constructs a GLPK-object.
        """

        glp_term_out(GLP_OFF)

        self.glpkLP = glp_create_prob()
        glp_create_index(self.glpkLP)
        glp_scale_prob(self.glpkLP, GLP_SF_AUTO)

        glp_set_obj_dir(self.glpkLP, GLP_MIN)
        glp_add_rows(self.glpkLP, len(self.row_names)) #  3=nrs number of constraints
        glp_add_cols(self.glpkLP, len(self.col_names)) #  3=nrs number of variables

        row_sign_mapping={'E':GLP_FX,'L':GLP_UP,'G':GLP_LO}
        for row_index in range(len(self.row_names)):
            glp_set_row_name(self.glpkLP, row_index+1, self.row_names[row_index]) #sets name of first row to p -> note that indexing starts at 1 not 0
            glp_set_row_bnds(self.glpkLP, row_index+1, row_sign_mapping[self.row_signs[row_index]], self.b[row_index], self.b[row_index]) #bounds the first row between lb 0 and ub 100

        bound_type_map={False:GLP_DB,True:GLP_FX}
        for col_index in range(len(self.col_names)):
            glp_set_col_name(self.glpkLP, col_index+1, self.col_names[col_index])
            glp_set_obj_coef(self.glpkLP, col_index+1, float(self.f[col_index])) #sets name of first variables objective coefficient to to 10 -> note that indexing starts at 1 not 0
            lb=self.LB[col_index]
            ub=self.UB[col_index]
            glp_set_col_bnds(self.glpkLP, col_index+1, bound_type_map[lb==ub], float(lb), float(ub))

        nonzero_row_inds=list(self.A.row)
        nonzero_col_inds=list(self.A.col)
        nonzero_data=list(self.A.data)

        nonzero_row_indices = intArray(len(nonzero_data)+1)
        nonzero_col_indices = intArray(len(nonzero_data)+1)
        nonzero_coefficients = doubleArray(len(nonzero_data)+1)

        for i in range(len(nonzero_data)):
            nonzero_row_indices[i+1]=int(nonzero_row_inds[i]+1)
            nonzero_col_indices[i+1]=int(nonzero_col_inds[i]+1)
            nonzero_coefficients[i+1]=nonzero_data[i]

        glp_load_matrix(self.glpkLP, len(nonzero_data), nonzero_row_indices, nonzero_col_indices, nonzero_coefficients)

        glp_scale_prob(self.glpkLP,GLP_SF_EQ) #GLP_SF_EQ
        #GLP_SF_GM perform geometric mean scaling;
        #GLP_SF_EQ perform equilibration scaling;
        #GLP_SF_2N round scale factors to nearest power of two;
        #GLP_SF_SKIP skip scaling, if the problem is well scaled.

        self.glpk_simplex_params=glp_smcp()
        setattr(self.glpk_simplex_params, "tol_bnd", 1e-9)
        setattr(self.glpk_simplex_params, "tol_dj", 1e-9)
        setattr(self.glpk_simplex_params, "tol_piv", 1e-10)
        glp_init_smcp(self.glpk_simplex_params)

    def solve_lp(self):
        """
        solves linear problem with swiglpk
        """
        message=glp_simplex(self.glpkLP, self.glpk_simplex_params)

    def return_solution_status(self) -> dict:
        """
        Returns solution status, after solving the problem.

        Returns
        ----------
        str : solution status
            optimal, feasible or infeasible
            solution stati returned by solvers:
                optimal:
                    swiglpk -> GLP_OPT
                    cplex -> 1
                feasible:
                    swiglpk -> GLP_FEAS
                feasible_only_before_unscaling:
                    cplex ->  5
        """
        status=glp_get_status(self.glpkLP)
        if status in SOLVER_MAP["swiglpk"]["optimal"]:
            return("optimal")
        elif status in SOLVER_MAP["swiglpk"]["feasible"]:
            return("feasible")
        elif status in SOLVER_MAP["swiglpk"]["feasible_only_before_unscaling"]:
            return("feasible_only_before_unscaling")
        else:
            return("infeasible")

    def return_primal_values(self) -> dict:
        """
        Returns (primal)solution vector.

        Returns
        ----------
        dict : Problem variables as keys and solution values as values.
        """
        return({self.col_names[i]: glp_get_col_prim(self.glpkLP, i+1) for i in range(len(self.col_names))})

    def return_dual_values(self) -> dict:
        """
        Returns (dual) solution vector.

        Returns
        ----------
        dict : Problem constraints as keys and dual values (shadow prices) as values.
        """
        return({self.row_names[i]: glp_get_row_dual(self.glpkLP, i+1) for i in range(len(self.row_names))})

    def return_objective_value(self) -> float:
        """
        Returns objective value, after solving the problem.

        Returns
        ----------
        float : objective value
        """
        return(glp_get_obj_val(self.glpkLP))

    def unscale_lp():
        """
        Unscale LP
        """
        pass

    def get_constraint_types(self, constraints: list = []) -> dict:
        """
        Returns constraint types for specified rows.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        constraints : list of row names.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            type identification-characters as values.
        """
        return({row:self.row_signs[self.row_names.index(row)] for row in constraints})

    def set_constraint_types(self, inputDict: dict):
        """
        Sets the type of specified constraints.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and constraint type as values
        """
        row_sign_mapping={'E':GLP_FX,'L':GLP_UP,'G':GLP_LO}
        for row in list(inputDict.keys()):
            row_index=self.row_names.index(row)
            glp_set_row_bnds(self.glpkLP, row_index+1, row_sign_mapping[inputDict[row]], self.b[row_index], self.b[row_index]) #bounds the first row between lb 0 and ub 100
            self.row_signs[row_index]=inputDict[row]

    def set_objective(self, inputDict: dict):
        """
        Sets coefficients in objective function to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and objective coefficient as values
        """
        for col in list(inputDict.keys()):
            glp_set_obj_coef(self.glpkLP, self.col_names.index(col)+1, float(inputDict[col])) #sets name of first variables objective coefficient to to 10 -> note that indexing starts at 1 not 0
            self.f[self.col_names.index(col)]=float(inputDict[col])

    def set_right_hand_side(self,inputDict):
        """
        Sets the rhs-value of specified constraints.

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and rhs-values as values.
        """
        row_sign_mapping={'E':GLP_FX,'L':GLP_UP,'G':GLP_LO}
        for row in list(inputDict.keys()):
            glp_set_row_bnds(self.glpkLP, self.row_names.index(row)+1, row_sign_mapping[self.row_signs[self.row_names.index(row)]], float(inputDict[row]), float(inputDict[row])) #bounds the first row between lb 0 and ub 100
            self.b[self.row_names.index(row)]=float(inputDict[row])

    def set_problem_coefficients(self, input: list):
        """
        Sets coefficients in constraint matrix.

        Parameters
        ----------
        inputDict : dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """
        A_array=self.A.toarray()
        for tuple in input:
            A_array[tuple[0],tuple[1]]=tuple[2]
        self.A=coo_matrix(A_array)
        self.build_glpk_lp()

    def set_lb(self, inputDict: dict):
        """
        Sets variable lower bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and lower bounds as values
        """
        for col in list(inputDict.keys()):
            if inputDict[col]!=self.UB[self.col_names.index(col)]:
                glp_set_col_bnds(self.glpkLP, self.col_names.index(col)+1, GLP_DB, float(inputDict[col]), float(self.UB[self.col_names.index(col)]))
            else:
                glp_set_col_bnds(self.glpkLP, self.col_names.index(col)+1, GLP_FX, float(inputDict[col]), float(inputDict[col]))
            self.LB[self.col_names.index(col)]=float(inputDict[col])

    def set_ub(self, inputDict: dict):
        """
        Sets variable upper bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and upper bounds as values
        """
        for col in list(inputDict.keys()):
            if inputDict[col]!=self.LB[self.col_names.index(col)]:
                glp_set_col_bnds(self.glpkLP, self.col_names.index(col)+1, GLP_DB, float(self.LB[self.col_names.index(col)]),float(inputDict[col]))
            else:
                glp_set_col_bnds(self.glpkLP, self.col_names.index(col)+1, GLP_FX, float(inputDict[col]), float(inputDict[col]))
            self.UB[self.col_names.index(col)]=float(inputDict[col])

class _SolverCPLEX(_Solver):

    @property
    def name(self):
        return("cplex")

    def build_lp(self):
        """
        Constructs a CPLEX-object.
        """
        lhs = self.A.tolil()
        rows = [SparsePair(nz_ind, data) for nz_ind, data in zip(lhs.rows, lhs.data)]
        # define problem
        cpxLP = Cplex()
        cpxLP.variables.add(obj=list(self.f), ub=list(self.UB),
                            lb=list(self.LB), names=list(self.col_names))
        cpxLP.linear_constraints.add(lin_expr=rows,
                                     rhs=self.b,
                                     senses=self.row_signs,
                                     names=self.row_names)
        cpxLP.objective.set_sense(cpxLP.objective.sense.minimize)
        cpxLP.parameters.feasopt.tolerance.set(1e-9)
        cpxLP.parameters.simplex.tolerances.feasibility.set(1e-9)
        cpxLP.parameters.simplex.tolerances.optimality.set(1e-9)
        cpxLP.parameters.simplex.tolerances.markowitz.set(0.1)
        cpxLP.parameters.barrier.convergetol.set(1e-9)
        cpxLP.parameters.read.scale.set(1)
        cpxLP.set_results_stream(None)
        cpxLP.set_log_stream(None)
        cpxLP.set_warning_stream(None)
        self.cplexLP = cpxLP

    def solve_lp(self):
        """
        solves linear problem with CPLEX
        """
        self.cplexLP.solve()

    def return_solution_status(self) -> dict:
        """
        Returns solution status, after solving the problem.

        Returns
        ----------
        str : solution status
            optimal, feasible or infeasible
            solution stati returned by solvers:
                optimal:
                    cplex -> 1
                feasible_only_before_unscaling:
                    cplex ->  5
        """
        status=self.cplexLP.solution.get_status()
        if status in SOLVER_MAP["cplex"]["optimal"]:
            return("optimal")
        elif status in SOLVER_MAP["cplex"]["feasible"]:
            return("feasible")
        elif status in SOLVER_MAP["cplex"]["feasible_only_before_unscaling"]:
            return("feasible_only_before_unscaling")
        else:
            return("infeasible")

    def return_primal_values(self) -> dict:
        """
        Returns (primal)solution vector.

        Returns
        ----------
        dict : Problem variables as keys and solution values as values.
        """
        return(dict(zip(self.col_names, self.cplexLP.solution.get_values())))

    def return_dual_values(self) -> dict:
        """
        Returns (dual) solution vector.

        Returns
        ----------
        dict : Problem constraints as keys and dual values (shadow prices) as values.
        """
        return(dict(zip(self.row_names, self.cplexLP.solution.get_dual_values())))

    def return_objective_value(self) -> float:
        """
        Returns objective value, after solving the problem.

        Returns
        ----------
        float : objective value
        """
        return(self.cplexLP.solution.get_objective_value())

    def unscale_lp(self):
        """
        Unscale LP
        """
        self.cplexLP.parameters.read.scale.set(-1)

    def get_constraint_types(self, constraints: list = []) -> dict:
        """
        Returns constraint types for specified rows.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        constraints : list of row names.

        Returns
        ----------
        dict
            Dictionary with constraint-IDs as keys and
            type identification-characters as values.
        """
        return(dict(zip(constraints, [self.cplexLP.linear_constraints.get_senses(self.rowIndicesMap[c]) for c in constraints])))

    def set_constraint_types(self, inputDict: dict):
        """
        Sets the type of specified constraints.
        (E: = ; L: <= ; G:>=)

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and constraint type as values
        """
        ##Update in cplex.Cplex LP##
        self.cplexLP.linear_constraints.set_senses(list(zip(inputDict.keys(), inputDict.values())))
        ##Transfer changes to rbatools.LP object##
        self.row_signs = self.cplexLP.linear_constraints.get_senses()

    def set_objective(self, inputDict: dict):
        """
        Sets coefficients in objective function to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and objective coefficient as values
        """
        if list(inputDict.keys()):
            self.cplexLP.objective.set_linear(zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())]))
            ##Transfer changes to rbatools.LP object##
            self.f = self.cplexLP.objective.get_linear()

    def set_right_hand_side(self,inputDict):
        """
        Sets the rhs-value of specified constraints.

        Parameters
        ----------
        inputDict : dict
            Constraint IDs as keys and rhs-values as values.
        """
        if list(inputDict.keys()):
            ##Update in cplex.Cplex LP##
            self.cplexLP.linear_constraints.set_rhs(list(zip(list(inputDict.keys()), [float(i) for i in list(inputDict.values())])))
            ##Transfer changes to rbatools.LP object##
            self.b = self.cplexLP.linear_constraints.get_rhs()

    def set_problem_coefficients(self, input: list):
        """
        Sets coefficients in constraint matrix.

        Parameters
        ----------
        inputDict : dict
            Dictionary with index tuples as keys and
            matrix coefficients as values.
        """
        self.cplexLP.linear_constraints.set_coefficients(input)
        self.A = _auxiliary_functions.convert_cplex_matrix_to_sparse(self)

    def set_lb(self, inputDict: dict):
        """
        Sets variable lower bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and lower bounds as values
        """
        #if list(inputDict.keys()):
            #self.cplexLP.variables.set_lower_bounds([(i,inputDict[i]) for i in inputDict.keys()])
            #self.UB = self.cplexLP.variables.get_lower_bounds()
        for i in inputDict.keys():
            self.cplexLP.variables.set_lower_bounds(i,float(inputDict[i]))
            self.LB[self.col_names.index(i)]=float(inputDict[i])

    def set_ub(self, inputDict: dict):
        """
        Sets variable upper bounds to specified values.

        Parameters
        ----------
        inputDict : dict
            Variable IDs as keys and upper bounds as values
        """
        #if list(inputDict.keys()):
            #self.cplexLP.variables.set_upper_bounds([(i,inputDict[i]) for i in inputDict.keys()])
            #self.UB = self.cplexLP.variables.get_upper_bounds()
        for i in inputDict.keys():
            self.cplexLP.variables.set_upper_bounds(i,float(inputDict[i]))
            self.UB[self.col_names.index(i)]=float(inputDict[i])
