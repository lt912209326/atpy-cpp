#ifndef _CPPCOMPONENT_H_
#define _CPPCOMPONENT_H_

#include "cppelement.h"
#include "cppelementfactory.h"
class CppComponent
{
    public:
    CppElement* elem;
    bool reverse;
    bool patch;
    size_t position;
    double* e1;
    double* e2;
    matrix<double> tws;
    // matrix<double> cache_tws;
    double local[LOC_NUM];
    // 组件几何位置Gx,Gy,Gz,theta1,theta2,theta3,误差dx,dy,dz,rot1,rot2,rot3, 真空管孔径Ax,Ay
    // double location[12];
    // double misalign[6];
    CppComponent();
    CppComponent(CppElement* elem, const bool reverse, const size_t position0);
    CppComponent(const CppComponent& comp0);
    CppComponent& operator= (const CppComponent& comp0);
    ~CppComponent();
    int display(const int,const bool);
    int linearoptics(const double* twsin,const Status* stat, double* );
    int update_TransferMatrix(double* rin,const Status* stat);
    int get_chromaticities(const double* twsin,const Status* stat);
    int get_geometry(const double* localin, const Status* stat);
    int compute_TransferMatrix(const double* locin, const Status* stat);
    int track(double* rin, const Status* stat);
    int patch1(double* rin );
    int patch2(double* rin );

};

#endif