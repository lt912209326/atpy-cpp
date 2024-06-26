from ..interface.constants cimport*
from ..interface.cppast cimport*

cdef dict token_enum
cdef dict enum_token

cdef class Token:
    cdef:
        readonly str kind #TWS,KWD, LOC, GLB,FUN,ID,ASSIGN,DELAY,
        readonly str value
        readonly int line
        readonly int column

cdef class Lexer:
    cdef:
        size_t count
        size_t token_num
        size_t line_num
        size_t column
        readonly list tokens
        dict elems_index
        list codes
        
    
    
    cdef Token get_current_token(self)
    
    cdef get_next_token(self)
    
    cdef str check_next_token(self)
    
    cdef int tokenize(self, str code)
