#ifndef _CPPQUADRUPOLE_H_
#define _CPPQUADRUPOLE_H_

#include "cppelement.h"


class CppQuadrupole:public CppElement
{
    public:
    CppQuadrupole();

    CppQuadrupole(string name,double l, double k1,int nslice0=1);

    CppQuadrupole(CppQuadrupole& elem0):CppElement(elem0){}
    
    int update_Matrix(const bool reverse, const double* cod , const Status* stat);
    
    // inline int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);

    int track(double* rin, const Status* stat, const bool reverse);


};



#endif