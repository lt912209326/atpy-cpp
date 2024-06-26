#ifndef _CPPSEXTUPOLE_H_
#define _CPPSEXTUPOLE_H_

#include "cppelement.h"

class CppSextupole:public CppElement
{
    public:
    // double* DynamicM66;
    CppSextupole();

    CppSextupole(const string name,double l, double k2, int nslice=4);

    CppSextupole(CppSextupole& elem0):CppElement(elem0){
        DynamicM66=elem0.DynamicM66;
    }

    ~CppSextupole()
    {
        if(DynamicM66){
            free(DynamicM66);
            DynamicM66=nullptr;
        }
    }


    int update_Matrix(const bool reverse, const double* cod , const Status* stat);
    
    virtual int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);

    int track(double* rin, const Status* stat, const bool reverse);

    // int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);

    int thin(double* rin, const double rate,const Status* stat);

    int drift(double* rin, const double rate,const Status* stat);

    int update_DriftMatrix(double* rin,const double len,const Status* stat);
    
    int update_ThinMatrix(double* rin,const double kick, const Status* stat);


};


#endif