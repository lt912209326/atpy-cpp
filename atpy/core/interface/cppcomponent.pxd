from .cppelements cimport *
cdef extern from "../physics/beamline/cppcomponent.h"nogil:
    cdef cppclass CppComponent:
        CppElement* elem
        bint reverse
        bint patch
        size_t position
        double* e1
        double* e2
        matrix[double] tws
        double local[LOC_NUM]
        # 组件几何位置Gx,Gy,Gz,theta1,theta2,theta3,误差dx,dy,dz,rot1,rot2,rot3, 真空管孔径Ax,Ay
        # double location[12]
        # double misalign[6]
        CppComponent()except+
        CppComponent(CppElement* elem, const bint reverse, const unsigned int position0)except+
        CppComponent(const CppComponent& comp0)except+
        CppComponent& operator= (const CppComponent& comp0)except+
        #~CppComponent()
        int display(const int,const bint)except+
        int linearoptics(const double* twsin,const CppStatus* stat, double* )except+
        int update_TransferMatrix(double* rin,const CppStatus* stat)except+
        int get_chromaticities(const double* twsin,const CppStatus* stat)except+
        int get_geometry(const double* localin, const CppStatus* stat)except+
        int compute_TransferMatrix(const double* locin, const CppStatus* stat)except+
        int track(double* rin, const CppStatus* stat)except+
        int patch1(double* rin )except+
        int patch2(double* rin )except+
