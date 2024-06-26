from .constants cimport *
from .cppast cimport *
from libcpp.set cimport set as cpp_set

cdef extern from "../utils/cppoptimization.h"nogil:
    
    cdef cppclass IdTable:
        unordered_map[string,AST*] id_dict
        vector[string]   id_table
        size_t num_id
        IdTable()except+

    cdef struct Variables:
        unordered_map[string,AST*]        independent_vars
        unordered_map[string,AST*]          dependent_vars
        vector[string]  ordered_independ_var_names
        vector[string]  ordered_depend_var_names
        size_t  num_independent_vars
        size_t  num_dependent_vars
        Variables()except+
        void append(string name0, AST* var,bint independent)except+



    cdef struct Constraints:
        vector[AST*]   values
        vector[string] ordered_constraint_names
        unordered_map[size_t,bint] time_consuming_terms
        size_t num_constraint
        Constraints()except+
        void append(string name0, AST* constraint)except+



    cdef struct Optima:
        vector[AST*] values
        vector[string] ordered_optima_names
        unordered_map[size_t,bint] time_consuming_terms
        vector[double] minormax
        size_t num_optima
        Optima()except+
        void append(string name0, AST* optima,bint ismin)except+
    
    
    cdef struct ChromCorrector:
        bint iscorr1
        bint iscorr2
        string corrector1
        double aim_dQx
        vector[string] dependent_corr1
        vector[size_t] position1
        vector[size_t] unique_position1
        vector[double] coeff1
        vector[double] unique_coeff1
        string corrector2
        double aim_dQy
        vector[string] dependent_corr2
        vector[size_t] position2
        vector[size_t] unique_position2
        vector[double] coeff2
        vector[double] unique_coeff2
        ChromCorrector()except+