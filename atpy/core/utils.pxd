
from .interface.cppstructures cimport Status as CppStatus, STAT0
from .interface.constants cimport*

cdef list STAT_INDEPEND_FIELDS

cdef dict PARTICLES_MAP

cdef dict DEFAULT_STATUS

cdef class Status:
    cdef:
        readonly dict _dict
        readonly list bool_type_fields
        readonly list int_type_fields
        readonly list float_type_fields
        readonly list readonly_fields
        CppStatus stat
