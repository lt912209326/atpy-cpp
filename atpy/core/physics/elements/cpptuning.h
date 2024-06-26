#ifndef _CPPTUNING_H_
#define _CPPTUNING_H_

#include "cppelement.h"


class CppTuning:public CppElement
{
    public:
    bool reverse_mat;
    double* M66INV;
    CppTuning();

    CppTuning(const string name0, 
            double dnux, double dnuy, double betax1, double betay1, double alphax1, double alphay1, double etax1, double etapx1
            , double betax2, double betay2, double alphax2, double alphay2, double etax2, double etapx2);

    CppTuning(CppTuning& elem0):CppElement(elem0){
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

    ~CppTuning();


    int track(double* rin, const Status* stat, const bool reverse);

    int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout);
    int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);
    int update_Matrix(const bool reverse, const double* cod , const Status* stat);



};


#endif