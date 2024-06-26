#ifndef _CPPMARKER_CPP_
#define _CPPMARKER_CPP_

#include "cppmarker.h"



CppMarker::CppMarker():CppElement(){
    kind=MARKER;
}

CppMarker::CppMarker(const string name):CppElement(name){
    kind=MARKER;
}

CppMarker::CppMarker(CppMarker& elem0):CppElement(elem0){}



int CppMarker::update_Matrix(const bool reverse, const double* cod , const Status* stat){
    memcpy(M66, EYE66,sizeof(EYE66));
    return 0;
}

int CppMarker::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout){  
    memcpy(twsout,twsin,TWS_NUM*sizeof(double));
    twsout[DCHROMX]=0.0;
    twsout[DCHROMY]=0.0;
    twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    return 0;
}


int CppMarker::track(double* rin, const Status* stat, const bool reverse){return 0;}



// CppMarker::~CppMarker(){}

#endif