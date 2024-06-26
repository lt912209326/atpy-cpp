#ifndef _CPPELEMENT_CPP_
#define _CPPELEMENT_CPP_

#include "cppelement.h"


CppElement::CppElement():nslice(1){
    memcpy(M66,EYE66,sizeof(EYE66) );
    memset(T66,0,sizeof(T66));
    DynamicM66=nullptr;
    // memcpy(T66,)
}

CppElement::CppElement(const string name0):nslice(1),name(name0){
    memcpy(M66,EYE66,sizeof(EYE66) );
    memset(T66,0,sizeof(T66));
    DynamicM66=nullptr;
    // memcpy(T66,)
}


CppElement::CppElement(const CppElement& elem):kind(elem.kind),nslice(elem.nslice),name(elem.name),keywords(elem.keywords),values(elem.values)
{
    // kind=elem.kind;
    // nslice=elem.nslice;
    // name=elem.name;
    // keywords=elem.keywords;
    // values=elem.values;
    
    // cout<<"CppElement::CppElement,keywords pointer address:"<< &keywords[0]<<endl;
    memcpy(M66,elem.M66,sizeof(double)*36);
    memcpy(T66,elem.T66,sizeof(double)*36);
    DynamicM66=nullptr;
    // if(elem.DynamicM66);
}


CppElement& CppElement::operator=(const CppElement& elem)
{
    kind=elem.kind;
    nslice=elem.nslice;
    name=elem.name;
    keywords=elem.keywords;
    values=elem.values;
    memcpy(M66,elem.M66,sizeof(double)*36);
    memcpy(T66,elem.T66,sizeof(double)*36);
    DynamicM66=nullptr;
    return *this;
}



int CppElement::compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out)
{
    for(int j=0;j<6;j++){
        R66out[0*6+j] = M66[0][0]*R66in[0*6+j]+ M66[0][1]*R66in[1*6+j] ;
        R66out[1*6+j] = M66[1][0]*R66in[0*6+j]+ M66[1][1]*R66in[1*6+j];

        R66out[2*6+j] = M66[2][2]*R66in[2*6+j]+ M66[2][3]*R66in[3*6+j];
        R66out[3*6+j] = M66[3][2]*R66in[2*6+j]+ M66[3][3]*R66in[3*6+j];

        R66out[4*6+j] = M66[4][4]*R66in[4*6+j]+ M66[4][5]*R66in[5*6+j];
        R66out[5*6+j] = M66[5][4]*R66in[4*6+j]+ M66[5][5]*R66in[5*6+j];
    }
    return 0;
}


int CppElement::track(double* rin, const Status* stat, const bool reverse){
    return 0;
}


int matmul66(double* C,double* A,double* B){
    int i,j,k;
    double aik=0 ;
    const int M=6,N=6,P=6;
    for(i=0;i<M;i++){
        for(j=0;j<N;j++) C[i*N+j]=0.0;
        for(k=0;k<P;k++){
            aik=A[P*i+k];
            for(j=0;j<N;j++){
                C[i*N+j] += aik*B[k*N+j];
            }
        }
    }
    return 0;
}

#endif