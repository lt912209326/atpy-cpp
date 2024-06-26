#ifndef _CPPDRIFT_H_
#define _CPPDRIFT_H_

#include "cppelement.h"


class CppDrift:public CppElement
{
    public:
    CppDrift();

    CppDrift(const string name0, double len, size_t nslice0=1);

    CppDrift(const CppDrift& elem0):CppElement(elem0){}


    int track(double* rin, const Status* stat, const bool reverse);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);

    int update_Matrix(const bool reverse, const double* cod , const Status* stat);

};


class CppExactDrift:public CppDrift
{
    public:

    CppExactDrift(const string name0, double len, size_t nslice0=1);
    CppExactDrift(const CppExactDrift& elem0):CppDrift(elem0){};



};


#endif