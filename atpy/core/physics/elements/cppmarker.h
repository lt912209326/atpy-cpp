#ifndef _CPPMARKER_H_
#define _CPPMARKER_H_

#include "cppelement.h"


class CppMarker:public CppElement
{
    public:

    CppMarker();

    CppMarker(const string name);

    CppMarker(CppMarker&);

    int update_Matrix(const bool reverse, const double* cod , const Status* stat);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);

    int track(double* rin, const Status* stat, const bool reverse);

    // ~CppMarker();
};



#endif