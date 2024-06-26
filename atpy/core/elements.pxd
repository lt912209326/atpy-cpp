from libcpp.string cimport string
from .interface.cppelements cimport*
from .interface.constants cimport*


cdef class Line:
    cdef:
        readonly str  name
        readonly dict elems
        readonly list line
        readonly list expand
        readonly list reverse
        readonly list ordered_lines



cdef class Element:
    cdef:
        CppElement         *elem
        readonly str        name
        readonly str        kind
        readonly int        kind_code
        readonly list       keywords
        readonly list       eids

cdef class Marker(Element):
    pass


cdef class Drift(Element):
    pass

cdef class ExactDrift(Element):
    pass

cdef class Dipole(Element):
    pass


cdef class Quadrupole(Element):
    pass


cdef class Sextupole(Element):
    pass


cdef class Octupole(Element):
    pass


cdef class Tuning(Element):
    pass

cdef class Girder(Element):
    pass


cdef class RFCavity(Element):
    pass


cdef class Wiggler(Line):
    pass



