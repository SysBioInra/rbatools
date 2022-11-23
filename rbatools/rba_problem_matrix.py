# python 2/3 compatibility
from __future__ import division, print_function


import numpy
import scipy
import warnings
from scipy.sparse import (SparseEfficiencyWarning)
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix
import rbatools._auxiliary_functions as _auxiliary_functions
from rbatools._warnings_and_errors import *

warnings.simplefilter('ignore', SparseEfficiencyWarning)


class ProblemMatrix(object):
    """
    Class holding RBA-Linear Problem.
    Required Object type for linear problems (matrices), to be used in this RBA-API.
    Should include all fields, required for a linear problem in the COBRA format, as attributes.

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
    """

    def __init__(self):
        """
        Initiates with empty Problem. Attributes can be changed manually.
        """

        self.A = coo_matrix(numpy.array([]))
        self.b = numpy.array([])
        self.f = numpy.array([])
        self.LB = numpy.array([])
        self.UB = numpy.array([])
        self.row_signs = []
        self.row_names = []
        self.col_names = []
        self.map_indices()

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
            elif type(matrix.A) is scipy.sparse.lil_matrix:
                self.A = coo_matrix(matrix.A.astype('float64'))
            elif type(matrix.A) is scipy.sparse.csc_matrix:
                self.A = coo_matrix(matrix.A.astype('float64'))
            elif type(matrix.A) is scipy.sparse.csr_matrix:
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

    def scale_lhs(self, factor):
        self.A = self.A*factor

    def map_indices(self):
        """
        Generates rowIndicesMap and colIndicesMap from attributes.

        rowIndicesMap = Dictionary: {'constraint_name':index,...}
        colIndicesMap = Dictionary: {'variable_name':index,...}
        """
        self.rowIndicesMap = dict(zip(self.row_names, list(range(len(self.row_names)))))
        self.colIndicesMap = dict(zip(self.col_names, list(range(len(self.col_names)))))

    def a_to_lil(self):
        """
        Transforms constraint matrix (LHS) to scipy.sparse.lil_matrix format.
        """
        self.A = scipy.sparse.lil_matrix(self.A)

    def a_to_coo(self):
        """
        Transforms constraint matrix (LHS) to scipy.sparse.coo_matrix format.
        """
        self.A = coo_matrix(self.A)

    def to_float64(self):
        """
        Transforms all numerical attributes to number type numpy.float64.
        """
        self.A = self.A.astype('float64')
        self.b = self.b.astype('float64')
        self.f = self.f.astype('float64')
        self.LB = self.LB.astype('float64')
        self.UB = self.UB.astype('float64')
