#ifndef _CPPOCTUPOLE_H_
#define _CPPOCTUPOLE_H_

#include "cppelement.h"

class CppOctupole:public CppElement
{
public:
    public:
    // double* DynamicM66;

    CppOctupole();

    CppOctupole(const string name,double l, double k3,size_t nslice0=4);

    CppOctupole(CppOctupole& elem0):CppElement(elem0){
        DynamicM66=elem0.DynamicM66;
    }
    ~CppOctupole()
    {
        if(DynamicM66){
            free(DynamicM66);
            DynamicM66=nullptr;
        }
    }

    int update_Matrix(const bool reverse, const double* cod , const Status* stat);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);

    int track(double* rin, const Status* stat, const bool reverse);

private:

    int _thin_track(double* rin, const double kick,const Status* stat);

    int _drift_track(double* rin, const double len,const Status* stat);
    
    int _update_thin_matrix(double* rin,const double kick, const Status* stat);

    int _update_drift_matrix(double* rin,const double len,const Status* stat);


    // int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);




};



#endif