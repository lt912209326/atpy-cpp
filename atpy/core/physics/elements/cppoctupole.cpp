#ifndef _CPPOCTUPOLE_CPP_
#define _CPPOCTUPOLE_CPP_

#include "cppoctupole.h"


const double C1 =0.5/(2-pow(2.0,1/3)) ;
const double C2 = 0.5*(1-pow(2.0,1/3))/(2-pow(2.0,1/3));
const double D1 = 1/(2-pow(2.0,1/3));
const double D2 = -pow(2.0,1/3)/(2-pow(2.0,1/3) );

CppOctupole::CppOctupole():CppElement(){
    kind=OCTUPOLE;
    keywords={L,K3,NSLICE};
    values[NSLICE]=nslice;
    DynamicM66=nullptr;
    
    
}

CppOctupole::CppOctupole(const string name,double l, double k3,size_t nslice0):CppElement(name){
    kind=OCTUPOLE;
    keywords={L,K3,NSLICE};
    values[L]=l;
    values[K3]=k3;
    nslice=nslice0;
    values[NSLICE]=nslice;
    DynamicM66=nullptr;
    // update_Matrix(false,0);
}

int CppOctupole::update_Matrix(const bool reverse, const double* cod , const Status* stat){
    double dp0=stat->dp;

    memcpy(M66,EYE66,sizeof(M66));
    if(stat->misaligment || fabs(cod[5])>1e-8 ){
        // cout<<"CppSextupole::update_Matrix"<<"x: "<<cod[0]<<",y: "<<cod[2]<<endl;
        double dp0=cod[5];
        double len=values[L]/nslice, k3=values[K3];//(1+dp0);
        double len1=C1*len;
        double len2=C2*len;
        double kick1=D1*len*k3/6.0;
        double kick2=D2*len*k3/6.0;

        double rin[6]={0};

        if(!DynamicM66){
            DynamicM66=(double*)calloc(nslice*36,sizeof(double));
        }
        memcpy(rin,cod, 6*__SIZEOF_DOUBLE__);
        if(fabs(len)<1e-8 ){
            _thin_track( rin, k3/6.0, stat);
            _update_thin_matrix( rin, k3/6.0, stat);
        }
        else{
            for(int i=0;i<nslice;i++){
                    _update_drift_matrix(rin,len1,stat);
                _drift_track(rin, len1, stat);
                    _update_thin_matrix( rin, kick1, stat);
                _thin_track( rin, kick1, stat);
                    _update_drift_matrix(rin,len2,stat);
                _drift_track(rin, len2, stat);
                    _update_thin_matrix( rin, kick2, stat);
                _thin_track( rin, kick2, stat);
                    _update_drift_matrix(rin,len2,stat);
                _drift_track(rin, len2, stat);
                    _update_thin_matrix( rin, kick1, stat);
                _thin_track( rin, kick1, stat);
                    _update_drift_matrix(rin,len1,stat);
                _drift_track(rin, len1, stat);

                memcpy(DynamicM66+36*i,M66,36*sizeof(double));
            }

        }
    }
    else{
        double len=values[L];
        M66[0][1]=len;
        M66[2][3]=len;
        M66[4][5]=len;
    }
    return 0;
}


int CppOctupole::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout)
{
    double len=len_rate*values[L];
    // BETA,ALPHA,GAMMA 
    if(fabs(len)<1e-8){
        memcpy(twsout,twsin,TWS_NUM*sizeof(double));
        twsout[DCHROMX]=0.0;
        twsout[DCHROMY]=0.0;
        twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    }
    else{
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
        //ETA,ETAPX
        twsout[ETAX] = twsin[ETAX] + len*twsin[ETAPX];
        twsout[ETAPX]= twsin[ETAPX];
        twsout[DCHROMX]   =  0.0;
        twsout[DCHROMY]   =  0.0;
        twsout[CHROMX]   =  twsout[DCHROMX]+twsin[CHROMX];
        twsout[CHROMY]   =  twsout[DCHROMY]+twsin[CHROMY];
        twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    }
    return 0;
}


int CppOctupole::track(double* rin, const Status* stat, const bool reverse){
    double len=values[L]/nslice, k3=values[K3];//(1+stat->dp);,,,
    double len1=C1*len;
    double len2=C2*len;
    double kick1=D1*len*k3/6.0;
    double kick2=D2*len*k3/6.0;
    if(fabs(values[L])<1e-8){
        _thin_track( rin, k3/6.0, stat);
    }
    else{
        for(int i=0;i<nslice;i++){
            // _drift_track(rin, 0.5*len, stat);
            // _thin_track( rin, k3*len, stat);
            // _drift_track(rin, 0.5*len, stat);

            _drift_track(rin, len1, stat);
            _thin_track( rin, kick1, stat);
            _drift_track(rin, len2, stat);
            _thin_track( rin, kick2, stat);
            _drift_track(rin, len2, stat);
            _thin_track( rin, kick1, stat);
            _drift_track(rin, len1, stat);
        }
    }
    return 0;
}




int CppOctupole::_thin_track(double* rin, const double kick,const Status* stat){
    double x0=rin[0], y0=rin[2],dp=rin[5] ;
    double pnorm=1/(1+dp) ;
    rin[1]= pnorm*kick*x0*(x0*x0 - 3.0*y0*y0) + rin[1];
    rin[3]= pnorm*kick*y0*(-3.0*x0*x0 + y0*y0) + rin[3];
        
    return 0;
}

int CppOctupole::_drift_track(double* rin, const double len,const Status* stat){
    double pnorm=1/(1+rin[5]);
    rin[0]=rin[0]+len*rin[1]*pnorm;
    rin[2]=rin[2]+len*rin[3]*pnorm;
    rin[4]-=len*pnorm*0.5*(rin[1]*rin[1]+rin[3]*rin[3]);
    return 0;
}


    
int CppOctupole::_update_thin_matrix(double* rin,const double kick, const Status* stat){

    double dp0=rin[5],pnorm = 1.0/(1+dp0 );
    double x=rin[0], y=rin[2], dkxx= pnorm*kick*x*x, dkxy= -3.0*pnorm*kick*x*y, dkyy=pnorm*kick*y*y;

    M66[1][0]= dkxx*M66[0][0] + M66[1][0] +dkxy*M66[2][0];
    M66[1][1]= dkxx*M66[0][1] + M66[1][1] +dkxy*M66[2][1];
    M66[1][2]= dkxx*M66[0][2] + M66[1][2] +dkxy*M66[2][2];
    M66[1][3]= dkxx*M66[0][3] + M66[1][3] +dkxy*M66[2][3];
    M66[1][4]= dkxx*M66[0][4] + M66[1][4] +dkxy*M66[2][4];
    M66[1][5]= dkxx*M66[0][5] + M66[1][5] +dkxy*M66[2][5];

    M66[3][0]=  dkxy*M66[0][0] +dkyy*M66[2][0] + M66[3][0];
    M66[3][1]=  dkxy*M66[0][1] +dkyy*M66[2][1] + M66[3][1];
    M66[3][2]=  dkxy*M66[0][2] +dkyy*M66[2][2] + M66[3][2];
    M66[3][3]=  dkxy*M66[0][3] +dkyy*M66[2][3] + M66[3][3];
    M66[3][4]=  dkxy*M66[0][4] +dkyy*M66[2][4] + M66[3][4];
    M66[3][5]=  dkxy*M66[0][5] +dkyy*M66[2][5] + M66[3][5];

    return 0;
}


int CppOctupole::_update_drift_matrix(double* rin,const double len,const Status* stat){

    double factor=1.0/(1+rin[5]);
    double normL=factor*len;
    M66[0][0]=M66[0][0]+normL*M66[1][0];
    M66[0][1]=M66[0][1]+normL*M66[1][1];
    M66[0][2]=M66[0][2]+normL*M66[1][2];
    M66[0][3]=M66[0][3]+normL*M66[1][3];
    M66[0][4]=M66[0][4]+normL*M66[1][4];
    M66[0][5]=M66[0][5]+normL*M66[1][5];

    M66[2][0]=M66[2][0]+normL*M66[3][0];
    M66[2][1]=M66[2][1]+normL*M66[3][1];
    M66[2][2]=M66[2][2]+normL*M66[3][2];
    M66[2][3]=M66[2][3]+normL*M66[3][3];
    M66[2][4]=M66[2][4]+normL*M66[3][4];
    M66[2][5]=M66[2][5]+normL*M66[3][5];
    return 0;
}

#endif