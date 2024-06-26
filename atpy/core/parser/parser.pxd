
from  .lexer cimport Token,Lexer
from  ..interface.cppast cimport *
from  ..interface.cppelements cimport CppElement
from  ..interface.cppcomponent cimport CppComponent
from  ..interface.cppbeamline cimport CppBeamLine
from libcpp.string cimport string
from libcpp.pair cimport pair
cimport libcpp.vector as cpp_vector
cimport libcpp.algorithm as cpp_algorithm
# from  interface.constants cimport*

ctypedef pair[string,AST* ] newpair

cdef dict value2enum

cdef class Parser:  # 定义语法分析器的类
    cdef:
        CppBeamLine* line
        readonly Lexer       lexer
        readonly Token       current_token
        readonly  list elems_name
        readonly  dict eids
        readonly  list id_table
        readonly  bint isdatabase
        readonly  list vars_name
        readonly  list vars_elem_name
        readonly  list tokens
        readonly  list codes
        readonly  list term_names
        readonly  list ordered_var_positions
        
    
    cdef void _set_database(self,CppBeamLine* line)except*
        

    cdef void eat(self, str kind)except*

    cdef int position(self)except*

    cdef tuple parameter(self)except*

    cdef tuple get_bounds(self)except*

    cdef int get_eid(self)except*

    cdef void set_Identity(self)except*

    cdef void set_variable(self)except*

    cdef void set_constraint(self)except*

    cdef void set_optimize(self)except*
    
    cdef void set_chromaticity(self)except*



    cdef AST* property(self)except*
    
    cdef AST* function(self, str func)except*
        
    cdef AST* factor(self)except*
    
    cdef AST* term(self)except*
    
    cdef AST* expr(self)except*

    cdef double eval(self, str code)except*

    #cdef double parse(self)except*

    