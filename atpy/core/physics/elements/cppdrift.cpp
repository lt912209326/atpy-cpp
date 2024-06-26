#ifndef _CPPDRIFT_CPP_
#define _CPPDRIFT_CPP_

#include "cppdrift.h"


CppDrift::CppDrift():CppElement()
{
    kind=DRIFT;
    keywords={L,NSLICE};
    values[L]=0;
    // update_Matrix(false,0);
}

CppDrift::CppDrift(const string name0, const double len, size_t nslice0):CppElement(name0)
{
    kind=DRIFT;
    keywords={L,NSLICE};
    values[L]=len;
    values[NSLICE]=nslice0;
    nslice=nslice0;
    // update_Matrix(false,0);
}

int CppDrift::track(double* rin, const Status* stat, const bool reverse){
    double length=values[L];
    // double pnorm=1/(1+rin[5]);
    double Lpnorm=length/sqrt( sqr(1+rin[5]) -sqr(rin[1])-sqr(rin[3]) );
    // double px=rin[1];
    // double py=rin[3];

    rin[0]+=rin[1]*Lpnorm;
    rin[2]+=rin[3]*Lpnorm;
    // rin[0]+=length*rin[1]/sqrt(sqr(1/pnorm)-sqr(px)-sqr(py));
    // rin[2]+=length*rin[3]/sqrt(sqr(1/pnorm)-sqr(px)-sqr(py));
    // rin[4]+=length*pnorm*0.5*(rin[1]*rin[1]+rin[3]*rin[3]);
    rin[4]+=(1+rin[5])*Lpnorm - length;
    return 0;

}

int CppDrift::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout){
    double len=len_rate*values[L];
    // BETA,ALPHA,GAMMA 
    twsout[COORD]=twsin[COORD]+len;
    twsout[BETAX] = twsin[BETAX] -2*len*twsin[ALPHAX]+len*len*twsin[GAMMAX];
    twsout[ALPHAX]=                     twsin[ALPHAX]-len*twsin[GAMMAX];
    twsout[GAMMAX]= twsin[GAMMAX];
    twsout[BETAY] = twsin[BETAY] -2*len*twsin[ALPHAY]+len*len*twsin[GAMMAY];
    twsout[ALPHAY]=                     twsin[ALPHAY]-len*twsin[GAMMAY];
    twsout[GAMMAY]= twsin[GAMMAY];
    // NUX,NUY
    twsout[NUX] = twsin[NUX] + 0.5*INV_PI*(atan(twsin[GAMMAX]*len - twsin[ALPHAX]) + atan(twsin[ALPHAX]));
    twsout[NUY] = twsin[NUY] + 0.5*INV_PI*(atan(twsin[GAMMAY]*len - twsin[ALPHAY]) + atan(twsin[ALPHAY]));
    //ETA.ETAPX
    twsout[ETAX] = twsin[ETAX] + len*twsin[ETAPX];
    twsout[ETAPX]= twsin[ETAPX];
    twsout[DCHROMX]   =  0.0;
    twsout[DCHROMY]   =  0.0;
    twsout[CHROMX]   =  twsin[CHROMX];
    twsout[CHROMY]   =  twsin[CHROMY];
    twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    return 0;



}

int CppDrift::update_Matrix(const bool reverse, const double* cod , const Status* stat){
    double length=values[L];
    double dp0=cod[5];
    double pnorm=1/(1+dp0);
    // double pnorm=1/(1+cod[5]);
    // double px=cod[1];//pnorm;
    // double py=cod[3];//pnorm;
    double px=cod[1];
    double py=cod[3];

    M66[0][1]=length*pnorm;
    // cout<<"CppDrift::update_Matrix:M12: "<<M66[0][1]<<endl;
    M66[2][3]=length*pnorm;
    // M66[0][1]=length/sqrt(sqr(1/pnorm)-sqr(px)-sqr(py));
    // M66[2][3]=length/sqrt(sqr(1/pnorm)-sqr(px)-sqr(py));
    M66[4][5]=length;
    return 0;
}


CppExactDrift::CppExactDrift(const string name0, const double len, size_t nslice0):CppDrift(name0,len,nslice0)
{
    kind=EXACTDRIFT;
    // keywords={L,NSLICE};
    // values[L]=len;
    // values[NSLICE]=nslice0;
    // nslice=nslice0;
    // update_Matrix(false,0);
}



#endif