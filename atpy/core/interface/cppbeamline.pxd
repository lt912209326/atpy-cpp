from  .cppelements cimport *
from .cppcomponent cimport *
from .cppoptimization cimport *
from .cppstructures cimport RDTsCache,MultipoleData,STAT0
from libcpp.vector cimport vector

# from libcpp.unordered_map cimport unordered_map
# from libcpp.string cimport string
# from .cppstructures cimport Status as CppStatus


cdef extern from "../physics/beamline/cppbeamline.h"nogil:

    cdef cppclass CppBeamLine:
        CppElementFactory factory
        unordered_map[string, vector[size_t]] elems_position
        # unordered_map[string, AST*] name_table
        # vector[string] ordered_name_list
        IdTable  id_table;
        Variables vars 
        Constraints constraints 
        Optima optima
        ChromCorrector chrom_corrector
        # unordered_map[string, Assign*] var_table
        # vector[string] ordered_var_list

        cpp_set[size_t] ordered_changed_position
        unordered_map[string,vector[CppElement*] ] elems
        unordered_map[string,vector[CppComponent*] ] components
        # vector[CppElement*] elems
        vector[CppComponent*] line
        RDTsCache rdtcache
        # vector[MultipoleData] multdata
        # vector[size_t] multipole_position
        string  name
        size_t length,nelems                              #number of components in the beam line
        CppStatus stat                                #the calculation status of a beam line like radiation,quatumn excitation,driving term,slice,polarization,independent component,on-momentum
        double globals[GLB_NUM]          #global parameters like natural emittance,natural chromaticity...
        # TrackData  trackdata             #trackdata to record the position of all particles/turn
        CppBeamLine()except+
        CppBeamLine(const string name0, const size_t particle, const double energy,const size_t length0, CppStatus stat0)except+
        CppBeamLine(const CppBeamLine&)except+
        # CppBeamLine(CppBeamLine&)         #global parameters like natural emittance,natural chromaticity...
        #~CppBeamLine()
        int save(CppElement* )except-1
        int display(const int,const bint)except+
        int displayTransferMatrix()except+
        int insert(CppElement*,const int pos,const bint reverse)except+
        int append(CppElement*,const bint reverse)except+
        int get_position_at_s(double coordinate)except+
        int get_optics_at_s(double coordinate, double* twsout)except+
        int TwissPropagate()except+
        int calculate()except+                            #calculate beam line according to above status
        int correctChrom(double dQx1, double dQy1)except +
        int compute_off_momentum_local_twiss(const double* tws)except+
        int compute_off_momentum_twiss(double* tws, double dp, bint local_twiss)except+
        int compute_large_off_momentum_tunes(double dp)except+
        int compute_off_momentum_sum_square(double dp)except+
        int compute_theory_touscheklifetime_part()except+
        int computeSecondOrderChromaticities(const double)except+
        int initialoptics()except+
        int linearoptics()except+
        int computeGlobalProperties()except+
        int computeRDTs()except+
        int findClosedOrbit(double* pop)except+
        int recover_twiss()except+
        int compute_sliced_off_momentum_twiss()except+ # sextupole and octupole only
        int compute_off_momentum_RDTs()except+
        int register_RDTs(double dpi,bint valid)except+
        int recover_RDTs(double dp0)except+
        int track(double* beams, size_t nbegin, size_t nend, size_t nturn_begin, size_t nturn_end )except+
        double get_DA_area(double* area_datas )except+
        int paralleltrack(int np, int nturn, double* beams)except+
    
    void deallocate(CppBeamLine* line)
