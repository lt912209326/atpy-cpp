#ifndef _CPPDIPOLE_CPP_
#define _CPPDIPOLE_CPP_

#include "cppdipole.h"
#include "cppradiationintegrals.h"



CppDipole::CppDipole():CppElement(){
    kind=DIPOLE;
    keywords={L,ANGLE,K1,E1,E2,NSLICE};
    values[L]=values[ANGLE]=values[K1]=values[E1]=values[E2]=0.0;
    values[NSLICE]=nslice;
    reverse_mat=false;
    M66INV=nullptr;
    // update_Matrix(false,0);
        
}


CppDipole::CppDipole(string name,  double l, double angle, double k1,double e1,double e2, int nslice0):CppElement(name){
    kind=DIPOLE;
    keywords={L,ANGLE,K1,E1,E2,NSLICE};
    // this->name=name;
    values[L]=l;
    values[ANGLE]=angle;
    values[K1]=k1;
    values[E1]=e1;
    values[E2]=e2;
    nslice=nslice0;
    values[NSLICE]=nslice;
    reverse_mat=false;
    M66INV=nullptr;
    
    // update_Matrix(false,0);
        
}


int CppDipole::update_Matrix(const bool reverse, const double* cod , const Status* stat){
    
    double dp0=cod[5];
    double len = values[L];
    double pnorm=1/(1+dp0);
    double angle=values[ANGLE];
    double rhoinv=angle/values[L],Gx=rhoinv;
    double Fx=(rhoinv*rhoinv+values[K1])*pnorm;
    double Fy=-values[K1]*pnorm;

    double edge1=angle*values[E1] ;
    double edge2=angle*values[E2]  ;
    double tge1=0,tge2=0;

    if(reverse&& (edge1!=edge2) && (!M66INV) ){
        reverse_mat=true;
        M66INV =(double*)calloc(36,__SIZEOF_DOUBLE__);
        if(!M66INV)throw std::runtime_error("M66INV allocate faile in CppDipole::update_Matrix") ;
        memcpy(M66INV, EYE66, sizeof(EYE66) );
    }
    double kx,ky,Cx,Sx,Cy,Sy,Dx;
    kx=ky=Cx=Sx=Cy=Sy=Dx=0.0;

    
    // 水平方向参数
    if(Fx>0.0){
        kx = sqrt(Fx);    
        Cx=cos(kx*len);    
        Sx=sin(kx*len)/kx;     
        Dx=(1-Cx)/Fx;
    }
    else if(Fx<0.0){
        kx = sqrt(-Fx);    
        Cx=cosh(kx*len);    
        Sx=sinh(kx*len)/kx;    
        Dx=(1-Cx)/Fx;
    }
    else{
        Cx = 1;
        Sx=len; 
        Dx=0.5*len*len;
    }
    // 垂直方向参数
    if(Fy>0.0){
        ky = sqrt(Fy);    
        Cy=cos(ky*len);    
        Sy=sin(ky*len)/ky; 
    }
    else if(Fy<0.0){
        ky = sqrt(-Fy);    
        Cy=cosh(ky*len);    
        Sy=sinh(ky*len)/ky;   
    }
    else{
        Cy = 1;
        Sy=len; 
    }

    // 扇形磁铁区域
    M66[0][0]=Cx    ;
    M66[0][1]=Sx   ;
    M66[1][0]=-Fx*Sx    ;
    M66[1][1]=Cx    ;
    M66[0][5]=pnorm*Gx*Dx;
    M66[1][5]=Gx*Sx;
    
    M66[2][2]= Cy    ;
    M66[2][3]= Sy    ;
    M66[3][2]=-Fy*Sy    ;
    M66[3][3]= Cy    ;
    M66[4][0]=-M66[1][5];
    M66[4][1]=-M66[0][5];
    // 入口、出口边缘不同的DIPOLE反向传输矩阵
    if(M66INV){
        M66INV[0*6+0]=Cx    ;
        M66INV[0*6+1]=Sx   ;
        M66INV[1*6+0]=-Fx*Sx    ;
        M66INV[1*6+1]=Cx    ;
        M66INV[0*6+5]=pnorm*Gx*Dx;//*pnorm;
        M66INV[1*6+5]=Gx*Sx;//*pnorm;
        
        M66INV[2*6+2]= Cy    ;
        M66INV[2*6+3]= Sy    ;
        M66INV[3*6+2]=-Fy*Sy    ;
        M66INV[3*6+3]= Cy    ;
        M66INV[4*6+0]=-M66INV[1*6+5];
        M66INV[4*6+1]=-M66INV[0*6+5];
        
        M66INV[0*6+1]*=pnorm;
        M66INV[1*6+0]/=pnorm;
        M66INV[2*6+3]*=pnorm;
        M66INV[3*6+2]/=pnorm;
    }

    M66[0][1]*=pnorm;
    M66[1][0]/=pnorm;
    M66[2][3]*=pnorm;
    M66[3][2]/=pnorm;



    // 入口边缘场
    
    if(0.0!=edge1) {
        // cout<<"here in CppDipole::update_Matrix: "<<edge1<<endl;
        tge1=tan(edge1); 
        M66[0][0]+= rhoinv*tge1*M66[0][1];    
        M66[1][0]+= rhoinv*tge1*M66[1][1];     
        M66[4][0]+= rhoinv*tge1*M66[4][1];   

        M66[2][2]-= rhoinv*tge1*M66[2][3];    
        M66[3][2]-= rhoinv*tge1*M66[3][3];  }
    if (0.0!=edge2) {
        tge2=tan(edge2);
        M66[1][0]+= rhoinv*tge2*M66[0][0];    
        M66[1][1]+= rhoinv*tge2*M66[0][1];    
        M66[1][5]+= rhoinv*tge2*M66[0][5];    

        M66[3][2]-= rhoinv*tge2*M66[2][2];    
        M66[3][3]-= rhoinv*tge2*M66[2][3];  
    }
    // 入口、出口边缘不同的DIPOLE反向传输矩阵
    if(M66INV){
        if(0!=edge2){
            tge2=tan(edge2); 
            M66INV[0*6+0]+= rhoinv*tge2*M66INV[0*6+1];
            M66INV[1*6+0]+= rhoinv*tge2*M66INV[1*6+1];

            M66INV[2*6+2]-= rhoinv*tge2*M66INV[2*6+3];
            M66INV[3*6+2]-= rhoinv*tge2*M66INV[3*6+3];
        }
        if (0!=edge1) {
            tge1=tan(edge1);
            M66INV[1*6+0]+= rhoinv*tge1*M66INV[0*6+0];
            M66INV[1*6+1]+= rhoinv*tge1*M66INV[0*6+1];
            M66INV[1*6+5]+= rhoinv*tge1*M66INV[0*6+5];

            M66INV[3*6+2]-= rhoinv*tge1*M66INV[2*6+2];
            M66INV[3*6+3]-= rhoinv*tge1*M66INV[2*6+3];
        }

    }

    return 0;

}



int CppDipole::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout){


    double len = len_rate*values[L];
    double angle=values[ANGLE];
    double dp0=stat->dp0;
    double pnorm=1/(1.0+dp0 );
    double rhoinv=angle/values[L],Gx=rhoinv;
    double Fx=(rhoinv*rhoinv + values[K1])*pnorm;
    double Fy=-values[K1]*pnorm;
    double edge1=reverse ? angle*values[E2]: angle*values[E1] ;
    double edge2=reverse ? angle*values[E1] : angle*values[E2]  ;
    double kx,ky,Cx,Sx,Cy,Sy,Dx;
    kx=ky=Cx=Sx=Cy=Sy=Dx=0.0;

    twsout[COORD]=twsin[COORD]+len;
    double betax0,alphax0,gammax0,betay0,alphay0,gammay0,nux0,nuy0,etax0,etapx0;
    double betax1,alphax1,gammax1,betay1,alphay1,gammay1,nux1,nuy1,etax1,etapx1;
    double dchromx1, dchromx2, dchromy1, dchromy2;
    dchromx1 = dchromx2 = dchromy1 = dchromy2 = 0.0;

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

    if(fabs(len_rate-1)<1e-8){
        // cout<<"CppDipole::linearoptics: "<<endl;
        computeRadiationIntegrals(rhoinv,len,values[K1],edge1, edge2, betax0,alphax0,etax0,etapx0,glbout);
    }
    twsout[DCHROMX] = 0.0;
    twsout[DCHROMY] = 0.0;

    // 入口边缘场
    if(0!=edge1){
        double tge1=tan(edge1);
        alphax1= -rhoinv*tge1*betax0+alphax0;
        gammax1=sqr(rhoinv*tge1)*betax0 -2.0*rhoinv*tge1*alphax0 + gammax0;
        alphay1= rhoinv*tge1*betay0+alphay0;
        gammay1=sqr(rhoinv*tge1)*betay0 +2.0*rhoinv*tge1*alphay0 + gammay0;
        etapx1=rhoinv*tge1*etax0+etapx0;

        alphax0=alphax1;
        gammax0=gammax1;
        alphay0=alphay1;
        gammay0=gammay1;
        etapx0=etapx1;

        dchromx1 =  0.25*INV_PI*betax0*tge1*rhoinv ;

        dchromy1 = -0.25*INV_PI*betay0*tge1*rhoinv ;
        twsout[DCHROMX] = dchromx1;
        twsout[DCHROMY] = dchromy1;

    }
    
    // 水平方向参数
    if(Fx>0.0){
        kx = sqrt(Fx);    
        Cx=cos(kx*len);    
        Sx=sin(kx*len)/kx;     
        Dx=(1-Cx)/Fx;

    }
    else if(Fx<0.0){
        kx = sqrt(-Fx);    
        Cx=cosh(kx*len);    
        Sx=sinh(kx*len)/kx;    
        Dx=(1-Cx)/Fx;
    }
    else{
        Cx = 1;
        Sx=len; 
        Dx=0.5*len*len;
    }
    // 垂直方向参数
    if(Fy>0.0){
        ky = sqrt(Fy);    
        Cy=cos(ky*len);    
        Sy=sin(ky*len)/ky; 
    }
    else if(Fy<0.0){
        ky = sqrt(-Fy);    
        Cy=cosh(ky*len);    
        Sy=sinh(ky*len)/ky;   
    }
    else{
        Cy = 1;
        Sy=len; 
    }
    twsout[DCHROMX] += -0.125*INV_PI*( (betax0*Fx - gammax0 )*Sx*Cx + alphax0*(Cx*Cx- Fx*Sx*Sx) + (betax0*Fx + gammax0)*len - alphax0 );
    // else{
    //     twsout[DCHROMX]  +=  -0.125*INV_PI*(-Fy)*( betax0*len  - 2*len*len*alphax0 + gammax0*pow(len,3)/3 );
    // }
    twsout[DCHROMY] += -0.125*INV_PI*( (betay0*Fy - gammay0 )*Sy*Cy + alphay0*(Cy*Cy- Fy*Sy*Sy) + (betay0*Fy + gammay0)*len - alphay0 );

    // 扇形磁铁区域
    betax1 =       Cx*Cx*betax0  -           2.0*Cx*Sx*alphax0 +  Sx*Sx*gammax0;
    alphax1=    Fx*Cx*Sx*betax0  +  (Cx*Cx - Fx*Sx*Sx)*alphax0 -  Cx*Sx*gammax0;
    gammax1= Fx*Fx*Sx*Sx*betax0  +          2*Fx*Cx*Sx*alphax0 +  Cx*Cx*gammax0;

    betay1 =       Cy*Cy*betay0  -           2.0*Cy*Sy*alphay0 +  Sy*Sy*gammay0;
    alphay1=    Fy*Cy*Sy*betay0  +  (Cy*Cy - Fy*Sy*Sy)*alphay0 -  Cy*Sy*gammay0;
    gammay1= Fy*Fy*Sy*Sy*betay0  +          2*Fy*Cy*Sy*alphay0 +  Cy*Cy*gammay0;
        
    etax1=    Cx*etax0 + etapx0*Sx +  Gx*Dx;
    etapx1=  -etax0*Fx*Sx +   etapx0*Cx + Gx*Sx;
    nux1 = nux0 + (atan(gammax0*Sx/Cx - alphax0) + atan(alphax0))/PIx2;
    nuy1 = nuy0 + (atan(gammay0*Sy/Cy - alphay0) + atan(alphay0))/PIx2;


    betax0=betax1;
    alphax0=alphax1;
    gammax0=gammax1;
    betay0=betay1;
    alphay0=alphay1;
    gammay0=gammay1;
    nux0=nux1;
    nuy0=nuy1;
    etax0=etax1;
    etapx0=etapx1;

    if((1.0==len_rate) && (0!=edge2) )
    {
        double tge2=tan(edge2);
        alphax1= -rhoinv*tge2*betax0+alphax0;
        gammax1=sqr(rhoinv*tge2)*betax0 -2.0*rhoinv*tge2*alphax0 + gammax0;
        alphay1= rhoinv*tge2*betay0+alphay0;
        gammay1=sqr(rhoinv*tge2)*betay0 +2.0*rhoinv*tge2*alphay0 + gammay0;

        etapx1=rhoinv*tge2*etax0+etapx0;

        alphax0=alphax1;
        gammax0=gammax1;
        alphay0=alphay1;
        gammay0=gammay1;
        etapx0=etapx1;

        dchromx2 =  0.25*INV_PI*betax0*tge2*rhoinv ;
        dchromy2 = -0.25*INV_PI*betay0*tge2*rhoinv ;
        twsout[DCHROMX] += dchromx2;
        twsout[DCHROMY] += dchromy2;
        // if(stat->printout) cout<<dchromy1<<", "<<dchromy2 <<endl;
    }

    twsout[BETAX] =  betax0 ;
    twsout[ALPHAX]=  alphax0;
    twsout[GAMMAX]=  gammax0;
    twsout[BETAY] =  betay0 ;
    twsout[ALPHAY]=  alphay0;
    twsout[GAMMAY]=  gammay0;
    twsout[NUX]   =  nux0;
    twsout[NUY]   =  nuy0;
    twsout[ETAX]  =  etax0;
    twsout[ETAPX] =  etapx0;
    twsout[CHROMX]   = twsout[DCHROMX] + twsin[CHROMX];
    twsout[CHROMY]   = twsout[DCHROMY] + twsin[CHROMY];
    twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;

    return 0;
}


int CppDipole::compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out){
    double* M=nullptr;
    if(!reverse || !reverse_mat){
        M=&M66[0][0]; 
    }
    else{
        M=M66INV;
    }
    // cout<<"here is CppDipole::compute_TransferMatrix 0, reverse?: "<<reverse<<endl;
    for(int j=0;j<6;j++){
        R66out[R11+0*6+j] = M[0*6+0]*R66in[R11+0*6+j]+ M[0*6+1]*R66in[R11+1*6+j]+ M[0*6+5]*R66in[R11+5*6+j];
        R66out[R11+1*6+j] = M[1*6+0]*R66in[R11+0*6+j]+ M[1*6+1]*R66in[R11+1*6+j]+ M[1*6+5]*R66in[R11+5*6+j];

        R66out[R11+2*6+j] = M[2*6+2]*R66in[R11+2*6+j]+ M[2*6+3]*R66in[R11+3*6+j];
        R66out[R11+3*6+j] = M[3*6+2]*R66in[R11+2*6+j]+ M[3*6+3]*R66in[R11+3*6+j];

        R66out[R11+4*6+j] = M[4*6+0]*R66in[R11+0*6+j]+ M[4*6+1]*R66in[R11+1*6+j]+   M[4*6+4]*R66in[R11+4*6+j]+ M[4*6+5]*R66in[R11+5*6+j];
        R66out[R11+5*6+j] = M[5*6+0]*R66in[R11+0*6+j]+ M[5*6+1]*R66in[R11+1*6+j]+   M[5*6+4]*R66in[R11+4*6+j]+ M[5*6+5]*R66in[R11+5*6+j];
    }
    // cout<<"here is CppDipole::compute_TransferMatrix 10 !"<<endl;
    return 0;
}


int CppDipole::track(double* rin, const Status* stat, const bool reverse){
    double x=rin[0],px=rin[1],y=rin[2],py=rin[3],z=rin[4],dp=rin[5];
    
    double len = values[L];
    double pnorm=1/(1+dp);
    double angle=values[ANGLE];
    double rhoinv=angle/values[L],Gx=rhoinv;
    double Fx=(rhoinv*rhoinv+values[K1])*pnorm;
    double Fy=-values[K1]*pnorm;

    double edge1=reverse ? angle*values[E2]: angle*values[E1] ;
    double edge2=reverse ? angle*values[E1] : angle*values[E2]  ;

    double tge1=0,tge2=0;


    double kx,ky,Cx,Sx,Cy,Sy,Dx;
    kx=ky=Cx=Sx=Cy=Sy=Dx=0.0;

    // 入口边缘场
    if(0.0!=edge1){
        // cout<<"here in CppDipole::update_Matrix: "<<edge1<<endl;
        tge1=tan(edge1); 
        px+=rhoinv*tge1*rin[0];      
        py-= rhoinv*tge1*rin[2];    
    }
    
    // 水平方向参数
    if(Fx>0.0){
        kx = sqrt(Fx);    
        Cx=cos(kx*len);    
        Sx=sin(kx*len)/kx;     
        Dx=(1-Cx)/Fx;
    }
    else if(Fx<0.0){
        kx = sqrt(-Fx);    
        Cx=cosh(kx*len);    
        Sx=sinh(kx*len)/kx;    
        Dx=(1-Cx)/Fx;
    }
    else{
        Cx = 1;
        Sx=len; 
        Dx=0.5*len*len;
    }
    // 垂直方向参数
    if(Fy>0.0){
        ky = sqrt(Fy);    
        Cy=cos(ky*len);    
        Sy=sin(ky*len)/ky; 
    }
    else if(Fy<0.0){
        ky = sqrt(-Fy);    
        Cy=cosh(ky*len);    
        Sy=sinh(ky*len)/ky;   
    }
    else{
        Cy = 1;
        Sy=len;
    }

    // 扇形磁铁区域
    rin[0] =    Cx*x        + pnorm*Sx*px + pnorm*Gx*Dx*dp;
    rin[1] = -Fx*Sx/pnorm*x + Cx*px       + Gx*Sx*dp;
    rin[2] =    Cy*y        + pnorm*Sy*py;
    rin[3] = -Fy*Sy/pnorm*y + Cy*py;
    // rin[4] = ;
    // rin[5] = ;

    if (0.0!=edge2) {
        tge2=tan(edge2);  
        rin[1]+= rhoinv*tge2*rin[0];
        rin[3]-= rhoinv*tge2*rin[2];
    }


    // rin[0]=x; rin[1]=xp; rin[2]=y; rin[3]=yp; rin[4]=z; rin[5]=dp;
    return 0;
}

CppDipole::~CppDipole(){
    if(M66INV)
        free(M66INV);
        M66INV=nullptr;
}

#endif