#ifndef _CPPDIPOLE_H_
#define _CPPDIPOLE_H_

#include "cppelement.h"
#include <exception>
#include <new>



class CppDipole:public CppElement{
    public:
    bool reverse_mat;
    double* M66INV;
    CppDipole();
    
    CppDipole(string name,  double l,  double angle, double k1=0,double e1=0.5,double e2=0.5, int nslice=1);

    CppDipole(CppDipole& elem0):CppElement(elem0){
        reverse_mat=elem0.reverse_mat;
        if(reverse_mat){
            // reverse_mat=true;
            M66INV =(double*)calloc(36,__SIZEOF_DOUBLE__);
            if(!M66INV) throw std::runtime_error("M66INV allocate faile in CppDipole::CppDipole(CppDipole& elem0)") ;
            memcpy(M66INV,elem0.M66INV,36*__SIZEOF_DOUBLE__);
        }
        else{
            M66INV=nullptr;
        }
    }

    int update_Matrix(const bool reverse, const double* cod , const Status* stat);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);

    int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);

    int track(double* rin, const Status* stat, const bool reverse);
    
    ~CppDipole();

};

#endif