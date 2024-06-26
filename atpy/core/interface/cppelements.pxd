# cython: language_level=3
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from .cppstructures cimport Status as CppStatus

from .constants cimport*


cdef extern from "../physics/elements/cppelement.h"nogil:
    cdef cppclass CppElement:
        int kind
        size_t nslice
        string name
        vector[size_t] position
        vector[int] keywords
        unordered_map[int,double] values
        double M66[6][6]
        double T66[6][6]
        CppElement()except+
        CppElement(const string)except+
        CppElement(const unordered_map[int,double] arg)except+
        CppElement(const CppElement&)except+
        CppElement& operator=(const CppElement& elem)except+
        inline int set_keywords(int kwd,double value0)except+
        inline int display()except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int compute_TransferMatrix(const double* R66in,const bint reverse,double* R66out)except+
        #~CppElement()


cdef extern from "../physics/elements/cppmarker.h"nogil:
    cdef cppclass CppMarker(CppElement):
        CppMarker()except+
        CppMarker(const string name)except+
        CppMarker(CppMarker&)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+
        #~CppMarker()


cdef extern from "../physics/elements/cppdrift.h"nogil:
    cdef cppclass CppDrift(CppElement):
        CppDrift()except+
        CppDrift(const string name0, double len, size_t nslice0)except+
        CppDrift(const CppDrift& elem0)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
    
    cdef cppclass CppExactDrift(CppDrift):
        CppExactDrift(const string name0, double len, size_t nslice0)except+


cdef extern from "../physics/elements/cpptuning.h"nogil:
    cdef cppclass CppTuning(CppElement):
        CppTuning()except+
        CppTuning(const string name0, 
              double dnux, double dnuy, double betax1, double betay1, double alphax1, double alphay1, double etax1, double etapx1
              ,double betax2, double betay2, double alphax2, double alphay2, double etax2, double etapx2)except+
        CppTuning(const CppTuning& elem0)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int compute_TransferMatrix(const double* R66in,const bint reverse,double* R66out)except+



cdef extern from "../physics/elements/cppdipole.h"nogil:
    cdef cppclass CppDipole(CppElement):
        bint reverse_mat
        double* M66INV
        CppDipole()except+
        CppDipole(string name,  double l,  double angle, double k1,double e1,double e2, int nslice)except+
        CppDipole(CppDipole& elem0)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int compute_TransferMatrix(const double* R66in,const bint reverse,double* R66out)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+
        #~CppDipole()


cdef extern from "../physics/elements/cppquadrupole.h"nogil:
    cdef cppclass CppQuadrupole(CppElement):
        CppQuadrupole()except+
        CppQuadrupole(string name,double l, double k1,int nslice0)except+
        CppQuadrupole(CppQuadrupole& elem0)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+


cdef extern from "../physics/elements/cppsextupole.h"nogil:
    cdef cppclass CppSextupole(CppElement):
        double* DynamicM66
        CppSextupole()except+
        CppSextupole(const string name,double l, double k2, int nslice)except+
        CppSextupole(CppSextupole& elem0)except+
        #~CppSextupole()
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int compute_TransferMatrix(const double* R66in,const bint reverse,double* R66out)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+
        int thin(double* rin, const double rate,const CppStatus* stat)except+
        int drift(double* rin, const double rate,const CppStatus* stat)except+
        int update_DriftMatrix(double* rin,const double len,const CppStatus* stat)except+
        int update_ThinMatrix(double* rin,const double kick, const CppStatus* stat)except+


cdef extern from "../physics/elements/cppoctupole.h"nogil:
    cdef cppclass CppOctupole(CppElement):
        CppOctupole()except+
        CppOctupole(const string name,double l, double k3, int nslice)except+
        CppOctupole(CppOctupole& elem0)except+
        int update_Matrix(const bint reverse, const double* cod , const CppStatus* stat)except+
        int linearoptics(const double* twsin, const double len_rate, const CppStatus* stat,const bint reverse, double* twsout, double* glbout)except+
        int track(double* rin, const CppStatus* stat, const bint reverse)except+

cdef extern from "../physics/elements/cppelementfactory.h"nogil:
    cdef cppclass CppElementFactory:
        # 根据元件类型创建对应的元件对象
        CppElement *CreateElement(CppElement* elem)except+