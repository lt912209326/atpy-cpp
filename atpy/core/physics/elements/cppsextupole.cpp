#ifndef _CPPCOMBINATION24_CPP_
#define _CPPCOMBINATION24_CPP_

#include "cppsextupole.h"

const double C1 =0.5/(2-pow(2.0,1/3)) ;
const double C2 = 0.5*(1-pow(2.0,1/3))/(2-pow(2.0,1/3));
const double D1 = 1/(2-pow(2.0,1/3));
const double D2 = -pow(2.0,1/3)/(2-pow(2.0,1/3) );

CppSextupole::CppSextupole():CppElement()
{
    kind=SEXTUPOLE;
    keywords={L,K2,NSLICE};
    values[NSLICE]=nslice;
    DynamicM66=nullptr;
    
    
}

CppSextupole::CppSextupole(const string name,double l, double k2, int nslice0):CppElement(name)
{
    kind=SEXTUPOLE;
    keywords={L,K2,NSLICE};
    values[L]=l;
    values[K2]=k2;
    nslice=nslice0;
    values[NSLICE]=nslice;
    DynamicM66=nullptr;
    
    // update_Matrix(false,0);
    
}

int CppSextupole::update_Matrix(const bool reverse, const double* cod ,const Status* stat)
{

    memcpy(M66,EYE66,sizeof(M66));
    if(stat->misaligment || fabs(cod[5])>1e-8 ){
        // cout<<"CppSextupole::update_Matrix"<<"x: "<<cod[0]<<",y: "<<cod[2]<<endl;
        double dp0=cod[5];
        double len=values[L]/nslice, k2=values[K2];//(1+dp0);
        double len1=C1*len;
        double len2=C2*len;
        double kick1=D1*len*k2;
        double kick2=D2*len*k2;
        double rin[6]={0};

        // cout<<"CppSextupole::update_Matrix:"<<cod[5]<<endl;

        if(!DynamicM66){
            DynamicM66=(double*)calloc(nslice*36,sizeof(double));
        }

        memcpy(rin,cod, 6*__SIZEOF_DOUBLE__);
        // cout<<"CppSextupole::update_Matrix: "<<rin[5]<<endl;
        for(int i=0;i<nslice;i++){

            // update_DriftMatrix(rin,0.5*len,stat);
            // drift(rin, 0.5*len, stat);
            // update_ThinMatrix( rin,  len*k2, stat);
            // thin( rin, len*k2, stat);
            // update_DriftMatrix(rin,0.5*len,stat);
            // drift(rin, 0.5*len, stat);
            


            update_DriftMatrix(rin,len1,stat);
            drift(rin, len1, stat);
            update_ThinMatrix( rin, kick1, stat);
            thin( rin, kick1, stat);
            update_DriftMatrix(rin,len2,stat);
            drift(rin, len2, stat);
            update_ThinMatrix( rin, kick2, stat);
            thin( rin, kick2, stat);
            update_DriftMatrix(rin,len2,stat);
            drift(rin, len2, stat);
            update_ThinMatrix( rin, kick1, stat);
            thin( rin, kick1, stat);
            update_DriftMatrix(rin,len1,stat);
            drift(rin, len1, stat);

            memcpy(DynamicM66+36*i,M66,36*sizeof(double));
            
        }
        // cout<<"M66[0][1]: "<<M66[0][1]<<endl;
    }
    else{
        double len=values[L];
        M66[0][1]=len;
        M66[2][3]=len;
        M66[4][5]=len;
    }

    // memcpy(M66,EYE66,sizeof(M66));
    return 0;
}


    
int CppSextupole::linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout)
{
    double len=len_rate*values[L],dp0=stat->dp0;
    double Fx=values[K2]/(1+dp0), Fy=-Fx;
    double betax0,alphax0,gammax0,betay0,alphay0,gammay0,nux0,nuy0,etax0,etapx0;
    
    betax0 = twsin[BETAX];
    alphax0= twsin[ALPHAX];
    gammax0= twsin[GAMMAX];
    betay0 = twsin[BETAY];
    alphay0= twsin[ALPHAY];
    gammay0= twsin[GAMMAY];
    etax0=twsin[ETAX];
    etapx0=twsin[ETAPX];
    if((stat->misaligment || fabs(stat->dp0)>1e-12) && DynamicM66 ){
        double nthslice=len_rate*nslice;
        size_t nth=round(nthslice)-1;
        if(fabs(nthslice-round(nthslice))<1e-8){
            
            double *LocM66=DynamicM66+nth*36;
            // cout<<"CppSextupole::linearoptics: "<<nth<<" "<<round(nthslice)<<endl;
            // for(size_t m=0;m<6;m++){
            //     for(size_t n=0;n<6;n++){
            //         cout<<"  "<<LocM66[m*6+n];
            //     }
            //     cout<<endl;
            // }

            twsout[BETAX] =  pow(LocM66[0],2)*betax0 -2*LocM66[0]*LocM66[1]*alphax0 +pow(LocM66[1],2)*gammax0;
        
            // cout<<"CppSextupole::linearoptics: "<<twsout[BETAX]<<endl;

            twsout[ALPHAX] =  -LocM66[0]*LocM66[1*6+0]*betax0+ (1 + 2*LocM66[1]*LocM66[1*6+0])*alphax0 - LocM66[1]*LocM66[1*6+1]*gammax0;
        
            twsout[GAMMAX] =  pow(LocM66[1*6+0],2)*betax0 -2*LocM66[1*6+0]*LocM66[1*6+1]*alphax0 + pow(LocM66[1*6+1],2)*gammax0;

            twsout[BETAY] =  pow(LocM66[2*6+2],2) *betay0 -2*LocM66[2*6+2]*LocM66[2*6+3]*alphay0 +pow(LocM66[2*6+3],2)*gammay0;
        
            twsout[ALPHAY] =  -LocM66[2*6+2]*LocM66[3*6+2]*betay0+ (1 + 2*LocM66[2*6+3]*LocM66[3*6+2])*alphay0 - LocM66[2*6+3]*LocM66[3*6+3]*gammay0;
        
            twsout[GAMMAY] =  pow(LocM66[3*6+2],2)*betay0 -2*LocM66[3*6+2]*LocM66[3*6+2+1]*alphay0 + pow(LocM66[3*6+3],2)*gammay0;
            
            twsout[ETAX]=    LocM66[0]*etax0 + etapx0*LocM66[1];
            twsout[ETAPX]=  etax0*LocM66[1*6+0] +   etapx0*LocM66[1*6+1];
            LocM66=nullptr;
            
        }
    }
    else{

        // BETA,ALPHA,GAMMA 
        twsout[BETAX] = twsin[BETAX] -2*len*twsin[ALPHAX]+len*len*twsin[GAMMAX];
        twsout[ALPHAX]=                     twsin[ALPHAX]-len*twsin[GAMMAX];
        twsout[GAMMAX]= twsin[GAMMAX];
        twsout[BETAY] = twsin[BETAY] -2*len*twsin[ALPHAY]+len*len*twsin[GAMMAY];
        twsout[ALPHAY]=                     twsin[ALPHAY]-len*twsin[GAMMAY];
        twsout[GAMMAY]= twsin[GAMMAY];
        //ETA.ETAPX
        twsout[ETAX] = twsin[ETAX] + len*twsin[ETAPX];
        twsout[ETAPX]= twsin[ETAPX];
        // twsout[DCHROMX]   =  + 0.25*INV_PI*Fx*(0.25*etapx0*gammax0*pow(len,4) + (etax0*gammax0 -2*etapx0*alphax0 )/3*pow(len,3) + 0.5*(etapx0*betax0 -2*etax0*alphax0 )*sqr(len) + etax0*betax0*len );
        // // twsin[CHROMX] -0.0625*INV_PI*( 2.0*(betax0*Fx -gammax0 )*Sx*Cx + 2*alphax0*(Cx*Cx- Fx*Sx*Sx) + 2*(betax0*Fx + gammax0)*len - 2*alphax0 );
        // twsout[DCHROMY]   =  + 0.25*INV_PI*Fy*(0.25*etapx0*gammay0*pow(len,4) + (etax0*gammay0 -2*etapx0*alphay0 )/3*pow(len,3) + 0.5*(etapx0*betay0 -2*etax0*alphay0 )*sqr(len) + etax0*betay0*len );
        // // twsin[CHROMY] -0.0625*INV_PI*( 2.0*(betay0*Fy -gammay0 )*Sy*Cy + 2*alphay0*(Cy*Cy- Fy*Sy*Sy) + 2*(betay0*Fy + gammay0)*len - 2*alphay0 );
    }

    twsout[COORD]=twsin[COORD]+len;
    // NUX,NUY
    twsout[NUX] = twsin[NUX] + 0.5*INV_PI*(atan(twsin[GAMMAX]*len - twsin[ALPHAX]) + atan(twsin[ALPHAX]));
    twsout[NUY] = twsin[NUY] + 0.5*INV_PI*(atan(twsin[GAMMAY]*len - twsin[ALPHAY]) + atan(twsin[ALPHAY]));
    if(fabs(len)<1e-8){
        twsout[DCHROMX]   =  + 0.25*INV_PI*Fx*betax0*etax0;
        twsout[DCHROMY]   =  + 0.25*INV_PI*Fy*betay0*etax0;
    }else{
        twsout[DCHROMX]   =  + 0.25*INV_PI*Fx*(0.25*etapx0*gammax0*pow(len,4) + (etax0*gammax0 -2*etapx0*alphax0 )/3*pow(len,3) + 0.5*(etapx0*betax0 -2*etax0*alphax0 )*sqr(len) + etax0*betax0*len );
        twsout[DCHROMY]   =  + 0.25*INV_PI*Fy*(0.25*etapx0*gammay0*pow(len,4) + (etax0*gammay0 -2*etapx0*alphay0 )/3*pow(len,3) + 0.5*(etapx0*betay0 -2*etax0*alphay0 )*sqr(len) + etax0*betay0*len );
    }
    twsout[CHROMX]   =  twsout[DCHROMX]+twsin[CHROMX];
    twsout[CHROMY]   =  twsout[DCHROMY]+twsin[CHROMY];
    twsout[H0] =  twsout[GAMMAX]*sqr( twsout[ETAX] )+ 2* twsout[ALPHAX]*twsout[ETAX]*twsout[ETAPX] + twsout[BETAX]*sqr(twsout[ETAPX] ) ;
    return 0;
}


int CppSextupole::compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out)
{
    for(int j=0;j<6;j++){
        R66out[0*6+j] = M66[0][0]*R66in[0*6+j]+ M66[0][1]*R66in[1*6+j]+ M66[0][2]*R66in[2*6+j]+ M66[0][3]*R66in[3*6+j]+ M66[0][4]*R66in[4*6+j]+ M66[0][5]*R66in[5*6+j];
        //  +M66[0][4]*R66in[4*6+j]+ M66[0][5]*R66in[5*6+j] ;
        R66out[1*6+j] = M66[1][0]*R66in[0*6+j]+ M66[1][1]*R66in[1*6+j]+ M66[1][2]*R66in[2*6+j]+ M66[1][3]*R66in[3*6+j]+ M66[1][4]*R66in[4*6+j]+ M66[1][5]*R66in[5*6+j];
        //  +M66[1][4]*R66in[4*6+j]+ M66[1][5]*R66in[5*6+j] ;;

        R66out[2*6+j] = M66[2][2]*R66in[2*6+j]+ M66[2][3]*R66in[3*6+j];
        R66out[3*6+j] = M66[3][2]*R66in[2*6+j]+ M66[3][3]*R66in[3*6+j];

        R66out[4*6+j] = M66[4][4]*R66in[4*6+j]+ M66[4][5]*R66in[5*6+j];
        R66out[5*6+j] = M66[5][4]*R66in[4*6+j]+ M66[5][5]*R66in[5*6+j];
    }
    return 0;
}


int CppSextupole::track(double* rin, const Status* stat, const bool reverse){
    double len=values[L]/nslice, k2=values[K2];//(1+stat->dp);,,,
    double len1=C1*len;
    double len2=C2*len;
    double kick1=D1*len*k2;
    double kick2=D2*len*k2;

    if(fabs(values[L])<1e-8){
        thin( rin, k2, stat);
    }
    else{
        for(int i=0;i<nslice;i++){
            drift(rin, len1, stat);
            thin( rin, kick1, stat);
            drift(rin, len2, stat);
            thin( rin, kick2, stat);
            drift(rin, len2, stat);
            thin( rin, kick1, stat);
            drift(rin, len1, stat);
        }

    }
    return 0;
}


// int CppSextupole::compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out){
//     // 
//     memcpy()
//     return 0;
// }


int CppSextupole::thin(double* rin, const double kick,const Status* stat)
{
    double pnorm=1/(1+rin[5]);
    double kick2 = pnorm*kick;
    rin[1]= -0.5*kick2*rin[0]*rin[0] + rin[1] + 0.5*kick2*rin[2]*rin[2];
    rin[3]=kick2*rin[0]*rin[2] + rin[3];
    // rout[4]=rin[4]+length*rin[5];
    // rin[5]=rin[5];
        
    return 0;

}
int CppSextupole::drift(double* rin, const double len,const Status* stat)
{
    double pnorm=1/(1+rin[5]);
    // cout<<" CppSextupole::drift:rin[5]: "<<rin[5]<<endl;
    rin[0]=rin[0]+len*rin[1]*pnorm;
    rin[2]=rin[2]+len*rin[3]*pnorm;
    rin[4]-=len*pnorm*0.5*(rin[1]*rin[1]+rin[3]*rin[3]);
    
    return 0;
}


int CppSextupole::update_DriftMatrix(double* rin, const double len,const Status* stat)
{
    double factor=1.0/(1+rin[5]);
    // double factor=1.0/(1+stat->dp);
    // cout<<" CppSextupole::update_DriftMatrix:dp: "<<rin[5]<<endl;
    M66[0][0]=M66[0][0]+len*factor*M66[1][0];
    M66[0][1]=M66[0][1]+len*factor*M66[1][1];
    M66[0][2]=M66[0][2]+len*factor*M66[1][2];
    M66[0][3]=M66[0][3]+len*factor*M66[1][3];
    M66[0][4]=M66[0][4]+len*factor*M66[1][4];
    M66[0][5]=M66[0][5]+len*factor*M66[1][5];

    M66[2][0]=M66[2][0]+len*factor*M66[3][0];
    M66[2][1]=M66[2][1]+len*factor*M66[3][1];
    M66[2][2]=M66[2][2]+len*factor*M66[3][2];
    M66[2][3]=M66[2][3]+len*factor*M66[3][3];
    M66[2][4]=M66[2][4]+len*factor*M66[3][4];
    M66[2][5]=M66[2][5]+len*factor*M66[3][5];
    return 0;
}


int CppSextupole::update_ThinMatrix(double* rin, const double kick,const Status* stat)
{

    double x=rin[0], y=rin[2], dkx= -0.5*kick*x, dkxy= 0.5*kick*y, dky=kick*x;
    
    M66[1][0]= dkx*M66[0][0] + M66[1][0] +dkxy*M66[2][0];
    M66[1][1]= dkx*M66[0][1] + M66[1][1] +dkxy*M66[2][1];
    M66[1][2]= dkx*M66[0][2] + M66[1][2] +dkxy*M66[2][2];
    M66[1][3]= dkx*M66[0][3] + M66[1][3] +dkxy*M66[2][3];
    M66[1][4]= dkx*M66[0][4] + M66[1][4] +dkxy*M66[2][4];
    M66[1][5]= dkx*M66[0][5] + M66[1][5] +dkxy*M66[2][5];

    M66[3][0]=  dky*M66[2][0] + M66[3][0];
    M66[3][1]=  dky*M66[2][1] + M66[3][1];
    M66[3][2]=  dky*M66[2][2] + M66[3][2];
    M66[3][3]=  dky*M66[2][3] + M66[3][3];
    M66[3][4]=  dky*M66[2][4] + M66[3][4];
    M66[3][5]=  dky*M66[2][5] + M66[3][5];

    return 0;
}

#endif