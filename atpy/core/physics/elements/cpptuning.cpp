#ifndef _CPPTUNING_CPP_
#define _CPPTUNING_CPP_

#include "cpptuning.h"


CppTuning::CppTuning():CppElement()
{
    kind=TUNING;
    keywords={DNUX,DNUY,BETAX1,ALPHAX1,BETAY1,ALPHAY1,ETAX1,ETAPX1,BETAX2,ALPHAX2,BETAY2,ALPHAY2,ETAX2,ETAPX2 };
    values[DNUX]=values[DNUY]=values[BETAX1]=values[BETAY1]=values[ALPHAX1]=values[ALPHAY1]=values[ETAX1]=values[ETAPX1]=0;
    values[BETAX2]=values[BETAY2]=values[ALPHAX2]=values[ALPHAY2]=values[ETAX2]=values[ETAPX2]=0;
    reverse_mat=false;
    M66INV=nullptr;
    // update_Matrix(false,0);
}

CppTuning::CppTuning(const string name0, 
            double dnux, double dnuy, double betax1, double betay1, double alphax1, double alphay1, double etax1, double etapx1
            , double betax2, double betay2, double alphax2, double alphay2, double etax2, double etapx2):CppElement(name0)
{
    kind=TUNING;
    keywords={DNUX,DNUY,BETAX1,ALPHAX1,BETAY1,ALPHAY1,ETAX1,ETAPX1,BETAX2,ALPHAX2,BETAY2,ALPHAY2,ETAX2,ETAPX2 };
    values[DNUX]=dnux;
    values[DNUY]=dnuy;
    values[BETAX1]=betax1;
    values[BETAY1]=betay1;
    values[ALPHAX1]=alphax1;
    values[ALPHAY1]=alphay1;
    values[ETAX1]=etax1;
    values[ETAPX1]=etapx1;
    values[BETAX2]=betax2;
    values[BETAY2]=betay2;
    values[ALPHAX2]=alphax2;
    values[ALPHAY2]=alphay2;
    values[ETAX2]=etax2;
    values[ETAPX2]=etapx2;
    reverse_mat=false;
    M66INV=nullptr;
    // update_Matrix(false,0);
}


int CppTuning::track(double* rin, const Status* stat, const bool reverse){
    double pnorm=1/(1+rin[5]);
    double dnux=values[DNUX], dnuy=values[DNUY], betax1=values[BETAX1], betay1=values[BETAY1], alphax1=values[ALPHAX1], alphay1=values[ALPHAY1] ;
    double etax1=values[ETAX1], etapx1=values[ETAPX1] ;
    double betax2=values[BETAX2], betay2=values[BETAY2], alphax2=values[ALPHAX2], alphay2=values[ALPHAY2],etax2=values[ETAX2], etapx2=values[ETAPX2] ;
    double x0=rin[0], px0=rin[1], y0=rin[2], py0=rin[3], z0=rin[4], dp0=rin[5] ;
    double x1,px1,y1,py1,z1,dp1;
    double sinx=sin(PIx2*dnux), cosx=cos(PIx2*dnux), siny=sin(PIx2*dnuy), cosy=cos(PIx2*dnuy);
    double M11,M12,M21,M22, M16, M26;

    if(reverse){
        // alphax1=-values[ALPHAX1];
        // alphay1=-values[ALPHAY1] ;
        // etapx1=-etapx1;

        betax1=values[BETAX2];
        betay1=values[BETAY2];
        alphax1=-values[ALPHAX2];
        alphay1=-values[ALPHAY2] ;
        etax1=values[ETAX2];
        etapx1=-values[ETAPX2] ;

        betax2=values[BETAX1];
        betay2=values[BETAY1];
        alphax2=-values[ALPHAX1];
        alphay2=-values[ALPHAY1];
        etax2=values[ETAX1];
        etapx2=-values[ETAPX1] ;
    }
    double xp0,yp0,xp1,xp2,yp1,yp2;
    double sqrt_bx1 = sqrt(betax1 ), sqrt_by1 = sqrt(betay1 ),sqrt_bx2 = sqrt(betax2 ), sqrt_by2 = sqrt(betay2 );
    // M11 = cosx+alphax1*sinx;
    // M12 = betax1*sinx ;
    // M21 = -(1+sqr(alphax1))/betax1*sinx;
    // M22 = (cosx-alphax1*sinx) ;
    // M16 = (1-M11)*etax1 - M12*etapx1;
    // M26 = -M21*etax1 + (1-M22)*etapx1;

    xp0=pnorm*px0;
    yp0=pnorm*py0;
    
    x1=(x0 - etax1*dp0)/sqrt_bx1;
    xp1=(alphax1*x0+betax1*xp0 - (betax1*etapx1+alphax1*etax1 )*dp0 )/sqrt_bx1;
    y1=y0/sqrt_by1;
    yp1=(alphay1*y0+betay1*yp0)/sqrt_by1;
    z1=etapx1*x0-etax1*xp0+z0;
    dp1=dp0;

    rin[0] = sqrt_bx2*(cosx*x1 + sinx*xp1) + etax2*dp1;
    rin[1] = (1+dp1)*( (-(alphax2*cosx + sinx)*x1 + (cosx-alphax2*sinx)*xp1 )/sqrt_bx2 + etapx2*dp1);
    rin[2] = sqrt_by2*(cosy*y1 + siny*yp1);
    rin[3] = (1+dp1)*(( -(alphay2*cosy+siny)*y1 + (cosy-alphay2*siny) *yp1 )/sqrt_by2 ) ;
    rin[4] = ( -( (betax2*etapx2+alphax2*etax2)*cosx+etax2*sinx)*x1 + ( etax2*cosx - (betax2*etapx2+alphax2*etax2)*sinx )*xp1  )/sqrt_bx2 + z1;
    rin[5] = dp1;

    // rin[0]=           (cosx+alphax1*sinx)*x0 +          betax1*sinx*px0 + M16*rin[5];
    // rin[1]= -(1+sqr(alphax1))/betax1*sinx*x0 +  (cosx-alphax1*sinx)*px0 + M26*rin[5];

    // rin[2]=           (cosy+alphay1*siny)*y0 +          betay1*siny*py0;
    // rin[3]= -(1+sqr(alphay1))/betay1*siny*y0 +  (cosy-alphay1*siny)*py0;

    return 0;

}

int CppTuning::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout){
    // BETA,ALPHA,GAMMA 
    if(reverse){

        twsout[COORD]=twsin[COORD];
        twsout[BETAX] = values[BETAX1] ;
        twsout[ALPHAX]= -values[ALPHAX1];
        twsout[GAMMAX]= (1+ twsout[ALPHAX]*twsout[ALPHAX] )/twsout[BETAX];
        twsout[BETAY] = values[BETAY1] ;
        twsout[ALPHAY]= -values[ALPHAY1];
        twsout[GAMMAY]= (1+ twsout[ALPHAY]*twsout[ALPHAY] )/twsout[BETAY];
        // NUX,NUY
        twsout[NUX] = twsin[NUX] + values[DNUX] ;
        twsout[NUY] = twsin[NUY] + values[DNUY] ;
        //ETA.ETAPX
        twsout[ETAX] = values[ETAX1] ;
        twsout[ETAPX]= -values[ETAPX1];
        twsout[DCHROMX]   =  0.0;
        twsout[DCHROMY]   =  0.0;
        twsout[CHROMX]   =  twsin[CHROMX];
        twsout[CHROMY]   =  twsin[CHROMY];
        twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    }
    else{

        twsout[COORD]=twsin[COORD];
        twsout[BETAX] = values[BETAX2] ;
        twsout[ALPHAX]= values[ALPHAX2];
        twsout[GAMMAX]= (1+ twsout[ALPHAX]*twsout[ALPHAX] )/twsout[BETAX];
        twsout[BETAY] = values[BETAY2] ;
        twsout[ALPHAY]= values[ALPHAY2];
        twsout[GAMMAY]= (1+ twsout[ALPHAY]*twsout[ALPHAY] )/twsout[BETAY];
        // NUX,NUY
        twsout[NUX] = twsin[NUX] + values[DNUX] ;
        twsout[NUY] = twsin[NUY] + values[DNUY] ;
        //ETA.ETAPX
        twsout[ETAX] = values[ETAX2] ;
        twsout[ETAPX]= values[ETAPX2];
        twsout[DCHROMX]   =  0.0;
        twsout[DCHROMY]   =  0.0;
        twsout[CHROMX]   =  twsin[CHROMX];
        twsout[CHROMY]   =  twsin[CHROMY];
        twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    }
    return 0;

}


int CppTuning::compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out){
    double* M=nullptr;
    if(!reverse || !reverse_mat){
        M=&M66[0][0]; 
    }
    else{
        M=M66INV;
    }
    // cout<<"here is CppDipole::compute_TransferMatrix 0, reverse?: "<<reverse<<endl;
    for(int j=0;j<6;j++){
        R66out[R11+0*6+j] = M[0*6+0]*R66in[R11+0*6+j]+ M[0*6+1]*R66in[R11+1*6+j] + M[0*6+5]*R66in[R11+5*6+j];
        R66out[R11+1*6+j] = M[1*6+0]*R66in[R11+0*6+j]+ M[1*6+1]*R66in[R11+1*6+j] + M[1*6+5]*R66in[R11+5*6+j];

        R66out[R11+2*6+j] = M[2*6+2]*R66in[R11+2*6+j]+ M[2*6+3]*R66in[R11+3*6+j];
        R66out[R11+3*6+j] = M[3*6+2]*R66in[R11+2*6+j]+ M[3*6+3]*R66in[R11+3*6+j];
    }
    // cout<<"here is CppDipole::compute_TransferMatrix 10 !"<<endl;
    return 0;
}



int CppTuning::update_Matrix(const bool reverse, const double* cod , const Status* stat){

    double dnux=values[DNUX], dnuy=values[DNUY], betax1=values[BETAX1], betay1=values[BETAY1], alphax1=values[ALPHAX1], alphay1=values[ALPHAY1] ;
    double etax1=values[ETAX1], etapx1=values[ETAPX1 ];
    double sinx=sin(PIx2*dnux), cosx=cos(PIx2*dnux), siny=sin(PIx2*dnuy), cosy=cos(PIx2*dnuy);

    if(reverse && (!M66INV) ){
        reverse_mat=true;
        M66INV =(double*)calloc(36,__SIZEOF_DOUBLE__);
        if(!M66INV)throw std::runtime_error("M66INV allocate faile in CppTuning::update_Matrix") ;
        memcpy(M66INV, EYE66, sizeof(EYE66) );
    }
    M66[0][0]=(cosx+alphax1*sinx) ;
    M66[0][1]=betax1*sinx ;
    M66[1][0]=-(1+sqr(alphax1))/betax1*sinx ;
    M66[1][1]=(cosx-alphax1*sinx) ;
    M66[0][5]=(1.0-M66[0][0])*etax1 - M66[0][1]*etapx1;
    M66[1][5]=-M66[1][0]*etax1 + (1- M66[1][1])*etapx1;
    
    M66[2][2]=(cosy+alphay1*siny) ;
    M66[2][3]=betay1*siny;
    M66[3][2]=-(1+sqr(alphay1))/betay1*siny ;
    M66[3][3]=(cosy-alphay1*siny) ;

    if(reverse){
        M66INV[0*6+0]=(cosx-alphax1*sinx) ;
        M66INV[0*6+1]=betax1*sinx ;
        M66INV[1*6+0]=-(1+sqr(alphax1))/betax1*sinx ;
        M66INV[1*6+1]=(cosx+alphax1*sinx) ;
        M66INV[0*6+5]=(1.0-M66INV[0*6+0])*etax1 + M66INV[0*6+1]*etapx1;
        M66INV[1*6+5]=-M66INV[1*6+0]*etax1 - (1-M66INV[1*6+1])*etapx1;
        
        M66INV[2*6+2]=(cosy-alphay1*siny) ;
        M66INV[2*6+3]=betay1*siny;
        M66INV[3*6+2]=-(1+sqr(alphay1))/betay1*siny ;
        M66INV[3*6+3]=(cosy+alphay1*siny) ;

    }
    return 0;
}

CppTuning::~CppTuning(){
    if(M66INV)
        free(M66INV);
        M66INV=nullptr;
}






#endif