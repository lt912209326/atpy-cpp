#ifndef _CPPQUADRUPOLE_CPP_
#define _CPPQUADRUPOLE_CPP_

#include "cppquadrupole.h"



CppQuadrupole::CppQuadrupole():CppElement(){
    kind=QUADRUPOLE;
    keywords={L,K1,NSLICE};
    values[L]=0;
    values[K1]=0;
    values[NSLICE]=nslice;
}


CppQuadrupole::CppQuadrupole(std::string name,double l, double k1, int nslice0):CppElement(name){
    kind=QUADRUPOLE;
    keywords={L,K1,NSLICE};
    // this->name=name;
    values[L]=l;
    values[K1]=k1;
    nslice=nslice0;
    values[NSLICE]=nslice;
    // update_Matrix(false,0);
}



int CppQuadrupole::update_Matrix(const bool reverse, const double* cod , const Status* stat){
    double dp0=cod[5];
    double pnorm=1/(1+dp0);
    double len = values[L];
    double Fx=values[K1]/(1+dp0);
    double Fy=-Fx;
    double kx,ky,Cx,Sx,Cy,Sy;
    kx=ky=Cx=Sx=Cy=Sy=0.0;
    // 水平方向参数
    if(Fx>0.0){
        kx = sqrt(Fx);    
        Cx=cos(kx*len);    
        Sx=sin(kx*len)/kx;  
        ky = sqrt(-Fy);    
        Cy=cosh(ky*len);    
        Sy=sinh(ky*len)/ky;  
    }
    else if(Fx<0.0){
        kx = sqrt(-Fx);    
        Cx=cosh(kx*len);    
        Sx=sinh(kx*len)/kx;    
        ky = sqrt(Fy);    
        Cy=cos(ky*len);    
        Sy=sin(ky*len)/ky; 
    }
    else{
        Cx = 1;
        Sx=len; 
        Cy = 1;
        Sy=len; 
    }
    M66[0][0]=Cx    ;
    M66[0][1]=Sx*pnorm   ;
    M66[1][0]=-Fx*Sx/pnorm    ;
    M66[1][1]=Cx    ;
    
    M66[2][2]= Cy    ;
    M66[2][3]= Sy*pnorm    ;
    M66[3][2]=-Fy*Sy/pnorm    ;
    M66[3][3]= Cy    ;
    return 0;
}
   


int CppQuadrupole::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout){
    double len = len_rate*values[L];
    double Gx=0;
    double dp0=stat->dp0;
    double Fx=values[K1]/(1+dp0);
    double Fy=-Fx;
    double kx,ky,Cx,Sx,Cy,Sy;
    kx=ky=Cx=Sx=Cy=Sy=0.0;

    double betax0,alphax0,gammax0,betay0,alphay0,gammay0,nux0,nuy0,etax0,etapx0;
    twsout[COORD]=twsin[COORD]+len;
    betax0 = twsin[BETAX];
    alphax0= twsin[ALPHAX];
    gammax0= twsin[GAMMAX];
    betay0 = twsin[BETAY];
    alphay0= twsin[ALPHAY];
    gammay0= twsin[GAMMAY];
    nux0=twsin[NUX];
    nuy0=twsin[NUY];
    etax0=twsin[ETAX];
    etapx0=twsin[ETAPX];

    
    // 水平方向参数
    if(Fx>0.0){
        kx = sqrt(Fx);    
        Cx=cos(kx*len);    
        Sx=sin(kx*len)/kx;  
        ky = sqrt(-Fy);    
        Cy=cosh(ky*len);    
        Sy=sinh(ky*len)/ky;  
    }
    else if(Fx<0.0){
        kx = sqrt(-Fx);    
        Cx=cosh(kx*len);    
        Sx=sinh(kx*len)/kx;    
        ky = sqrt(Fy);    
        Cy=cos(ky*len);    
        Sy=sin(ky*len)/ky; 
    }
    else{
        Cx = 1;
        Sx=len; 
        Cy = 1;
        Sy=len; 
    }
    double Cx2=sqr(Cx);
    double Sx2=sqr(Sx);

    // 扇形磁铁区域
    twsout[BETAX] =       Cx*Cx*betax0  -           2.0*Cx*Sx*alphax0 +  Sx*Sx*gammax0;
    twsout[ALPHAX]=    Fx*Cx*Sx*betax0  +  (Cx*Cx - Fx*Sx*Sx)*alphax0 -  Cx*Sx*gammax0;
    twsout[GAMMAX]= Fx*Fx*Sx*Sx*betax0  +          2*Fx*Cx*Sx*alphax0 +  Cx*Cx*gammax0;

    twsout[BETAY] =       Cy*Cy*betay0  -           2.0*Cy*Sy*alphay0 +  Sy*Sy*gammay0;
    twsout[ALPHAY]=    Fy*Cy*Sy*betay0  +  (Cy*Cy - Fy*Sy*Sy)*alphay0 -  Cy*Sy*gammay0;
    twsout[GAMMAY]= Fy*Fy*Sy*Sy*betay0  +          2*Fy*Cy*Sy*alphay0 +  Cy*Cy*gammay0;
        
    twsout[ETAX]=    Cx*etax0 + etapx0*Sx +  Gx*(1-Cx);
    twsout[ETAPX]=  -etax0*Fx*Sx +   etapx0*Cx + Gx*Sx;
    twsout[NUX] = nux0 + (atan(gammax0*Sx/Cx - alphax0) + atan(alphax0))/PIx2;
    twsout[NUY] = nuy0 + (atan(gammay0*Sy/Cy - alphay0) + atan(alphay0))/PIx2;
    
    twsout[DCHROMX]   =  -0.125*INV_PI*( (betax0*Fx -gammax0 )*Sx*Cx + alphax0*(Cx*Cx- Fx*Sx*Sx) + (betax0*Fx + gammax0)*len - alphax0 );
    twsout[DCHROMY]   =  -0.125*INV_PI*( (betay0*Fy -gammay0 )*Sy*Cy + alphay0*(Cy*Cy- Fy*Sy*Sy) + (betay0*Fy + gammay0)*len - alphay0 );
    
    twsout[CHROMX]   =  twsout[DCHROMX]+twsin[CHROMX];
    twsout[CHROMY]   =  twsout[DCHROMY]+twsin[CHROMY];
    twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    // cout<<"H0: "<<twsout[H0]<<endl;
    return 0;
}




int CppQuadrupole::track(double* rin, const Status* stat, const bool reverse){
    
    double x=rin[0],px=rin[1],y=rin[2],py=rin[3],z=rin[4],dp=rin[5];
    double pnorm=1/(1+dp);
    double len = values[L];
    double Fx=values[K1]/(1+dp);
    double Fy=-Fx;
    double kx,ky,Cx,Sx,Cy,Sy;
    kx=ky=Cx=Sx=Cy=Sy=0.0;
    // 水平方向参数
    if(Fx>0.0){
        kx = sqrt(Fx);    
        Cx=cos(kx*len);    
        Sx=sin(kx*len)/kx;  
        ky = sqrt(-Fy);    
        Cy=cosh(ky*len);    
        Sy=sinh(ky*len)/ky;  
    }
    else if(Fx<0.0){
        kx = sqrt(-Fx);    
        Cx=cosh(kx*len);    
        Sx=sinh(kx*len)/kx;    
        ky = sqrt(Fy);    
        Cy=cos(ky*len);    
        Sy=sin(ky*len)/ky; 
    }
    else{
        Cx = 1;
        Sx=len; 
        Cy = 1;
        Sy=len; 
    }
    
    rin[0] =    Cx*x        + pnorm*Sx*px ;
    rin[1] = -Fx*Sx/pnorm*x + Cx*px       ;
    rin[2] =    Cy*y        + pnorm*Sy*py;
    rin[3] = -Fy*Sy/pnorm*y + Cy*py;
    return 0;
}


#endif