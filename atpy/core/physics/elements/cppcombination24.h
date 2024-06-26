#ifndef _CPPCOMBINED46_H_
#define _CPPCOMBINED46_H_

#include "cppelement.h"

class CppCombined46:public CppElement
{
    public:
    CppCombined46(double* varin)
    {
        kind=COMBINED;
        keywords={L,ANGLE,K1,E1,E2};
        values[L]=varin[0];
        values[ANGLE]=varin[1];
        values[K1]=varin[2];
        values[E1]=varin[3];
        values[E2]=varin[4];
        
        memcpy(&M66,&EYE66[0],sizeof(EYE66));
        memset(&T66,0,sizeof(T66));
        update_Matrix(false);
        
    }

    inline int update_Matrix(bool reverse)
    {
        double len = values[L];
        double angle=values[ANGLE];
        double rhoinv=angle/len;
        double K=rhoinv*rhoinv + values[K1];
        
        if(K>0.0)
        {
            double k = sqrt(K);
            double C=cos(k*len);
            double S=sin(k*len);
            M66[0][0]=C;
            M66[0][1]=S/k;
            M66[1][0]=-k*S;
            M66[1][1]=C;
            M66[0][5]=rhoinv/K*(1-C);
            M66[1][5]=rhoinv/k*S;
        }
        else if(K<0.0)
        {
            double k = sqrt(-K);
            double C=cosh(k*len);
            double S=sinh(k*len);
            double k1=sqrt(values[K1]);
            M66[0][0]=C;
            M66[0][1]=S/k;
            M66[1][0]=k*S;
            M66[1][1]=C;
            M66[0][5]=-rhoinv/K*(C-1);
            M66[1][5]=rhoinv/k*S;
            // K1<0
        }
        else if(0.0==K)
        {
            M66[0][0]=1.0;
            M66[0][1]=len;
            M66[1][0]=0;
            M66[1][1]=1.0;
            M66[0][5]=0.5*rhoinv*len*len;
            M66[1][5]=rhoinv*len;
        }


        // vertical matrix
        if(values[K1]>0.0)
        {
            // K1>0
            double k1=sqrt(values[K1]);
            double C=cosh(k1*len);
            double S=sinh(k1*len);
            M66[2][2]=C;
            M66[2][3]=S/k1;
            M66[3][2]=k1*S;
            M66[3][3]=C;

            // cout<<k1<<", "<<C<<endl;

            // cout<<M66[2][2]<<", "<<M66[2][3]<<endl;
            // cout<<M66[3][2]<<", "<<M66[3][3]<<endl;
        }
        else if(values[K1]<0.0)
        {
            // K1>0
            double k1=sqrt(-values[K1]);
            double C=cos(k1*len);
            double S=sin(k1*len);
            M66[2][2]=C;
            M66[2][3]=S/k1;
            M66[3][2]=-k1*S;
            M66[3][3]=C;
        }
        else if(0.0==values[K1])
        {
            // K1==0
            M66[2][2]=1.0;
            M66[2][3]=len;
            M66[3][2]=0;
            M66[3][3]=1.0;
        }

        if(0!=angle*values[E1])
        {
            double edge1=angle*values[E1];
            double tge1=tan(edge1);
            M66[0][0]+= rhoinv*tge1*M66[0][1];
            M66[1][0]+= rhoinv*tge1*M66[1][1];

            M66[2][2]-= rhoinv*tge1*M66[2][3];
            M66[3][2]-= rhoinv*tge1*M66[3][3];

        }

        if(0!=angle*values[E2])
        {
            double edge2=angle*values[E2];
            double tge2=tan(edge2);
            M66[1][0]+= rhoinv*tge2*M66[0][0];
            M66[1][1]+= rhoinv*tge2*M66[0][1];
            M66[1][5]+= rhoinv*tge2*M66[0][5];

            M66[3][2]-= rhoinv*tge2*M66[2][2];
            M66[3][3]-= rhoinv*tge2*M66[2][3];

        }
    }

    inline int linearoptics(const double* twsin, const double len_rate, const Status* stat, double* twsout)
    {
        double len = len_rate*values[L];
        double angle=values[ANGLE];
        double rhoinv=angle/values[L];
        double K=rhoinv*rhoinv + values[K1];

        double betax0,alphax0,gammax0,betay0,alphay0,gammay0,nux0,nuy0,etax0,etapx0;
        double betax1,alphax1,gammax1,betay1,alphay1,gammay1,nux1,nuy1,etax1,etapx1;
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

        if(0!=angle*values[E1])
        {
            double edge1=angle*values[E1];
            double tge1=tan(edge1);
            betax1 = betax0;
            alphax1= -rhoinv*tge1*betax0+alphax0;
            gammax1=rhoinv*rhoinv*tge1*tge1*betax0 -2.0*rhoinv*tge1*alphax0 + gammax0;
            betay1 = betay0;
            alphay1= rhoinv*tge1*betay0+alphay0;
            gammay1=rhoinv*rhoinv*tge1*tge1*betay0 +2.0*rhoinv*tge1*alphay0 + gammay0;
            nux1=nux0;
            nuy1=nuy0;
            etax1=etax0;
            etapx1=rhoinv*tge1*etax0+etapx0;

            betax0=betax1;
            alphax0=alphax1;
            gammax0=gammax1;
            betay0=betay1;
            alphay0=alphay1;
            gammay0=gammay1;
            nux0=nux1;
            nuy0=nuy1;
            etax0=etax0;
            etapx0=etapx0;

        }
        

        if(values[K1]>0)
        {
            // K=rhoinv*rhoinv+values[K1]>0
            
            double k = sqrt(K);
            // k in hor.
            double k1 = sqrt(values[K1]);
            // k1 in ver.
            double C=cos(k*len);
            double S=sin(k*len);
            double C2=cos(2*k*len);
            double CH=cosh(k1*len);
            double SH=sinh(k1*len);
            double CH2=cosh(2*k1*len);
            betax1 =   C*C*betax0 -2.0/k*C*S*alphax0+ S*S/(k*k)*gammax0;
            alphax1= k*C*S*betax0 +       C2*alphax0- C*S/k*gammax0;
            gammax1= (k*k)*S*S*betax0 +  2*k*C*S*alphax0+   C*C*gammax0;
            nux1 = nux0 +0.5*INV_PI*(atan(gammax0*tan(k*len)/k - alphax0) + atan(alphax0));
            betay1 =        CH*CH*betay0 -2.0/k1*CH*SH*alphay0+ SH*SH/(k1*k1)*gammay0;
            alphay1=    -k1*CH*SH*betay0 +         CH2*alphay0-      CH*SH/k1*gammay0;
            gammay1=(k1*k1)*SH*SH*betay0 -  2*k1*CH*SH*alphay0+         CH*CH*gammay0;
            nuy1 = nuy0 + 0.5*INV_PI*(atan((k1*betay0*CH + gammay0/k1*SH - alphay0*exp(k1*len))*exp(k1*len))-atan(k1*betay0 - alphay0));
            etax1=    C*etax0 + S/k*etapx0 + rhoinv/k*(1-C);
            etapx1=-k*S*etax0 +   C*etapx0 + rhoinv/k*S;
        }
        
        else if(values[K1]<0.0)
        {
            // k in hor.
            double k1 = sqrt(-values[K1]);
            // k1 in ver.
            double C=cos(k1*len);
            double S=sin(k1*len);
            double C2=cos(2*k1*len);
            betay1 =   C*C*betay0 -2.0/k1*C*S*alphay0+ S*S/(k1*k1)*gammay0;
            alphay1= k1*C*S*betay0 +       C2*alphay0- C*S/k1*gammay0;
            gammay1= k1*k1*S*S*betay0 +  2*k1*C*S*alphay0+   C*C*gammay0;
            nuy1 = nuy0 + 0.5*INV_PI*(atan(gammay0*tan(k1*len)/k1 - alphay0) + atan(alphay0));

            if(K>0)
            {
                // K>0
                double k=sqrt(K);
                double CP=cos(k*len);
                double SP=sin(k*len);
                double CP2=cos(2*k*len);
                betax1 =   CP*CP*betax0 -2.0/k*CP*SP*alphax0+ SP*SP/(k*k)*gammax0;
                alphax1= k*CP*SP*betax0 +       CP2*alphax0- CP*SP/k*gammax0;
                gammax1= (k*k)*SP*SP*betax0 +  2*k*CP*SP*alphax0+   CP*CP*gammax0;
                nux1 = nux0 +0.5*INV_PI*(atan(gammax0*tan(k*len)/k - alphax0) + atan(alphax0));
                etax1=    CP*etax0 + SP/k*etapx0 + rhoinv/k*(1-CP);
                etapx1=-k*SP*etax0 +   CP*etapx0 + rhoinv/k*SP;
            }
            else if(K<0)
            {
                // K<0
                double k=sqrt(-K);
                double CH=cosh(k*len);
                double SH=sinh(k*len);
                double CH2=cosh(2*k*len);
                betax1 =     CH*CH*betax0 -2.0/k*CH*SH*alphax0+ SH*SH/(k*k)*gammax0;
                alphax1=  -k*CH*SH*betax0 +        CH2*alphax0-     CH*SH/k*gammax0;
                gammax1= k*k*SH*SH*betax0 -  2*k*CH*SH*alphax0+       CH*CH*gammax0;
                nux1 = nux0 + 0.5*INV_PI*(atan(( (k*betax0-2*alphax0)*CH + gammax0/k*SH)*exp(k*len)+alphax0)
                                        - atan(k*betax0 - alphax0));
                etax1=    CH*etax0 + SH/k*etapx0 + rhoinv/k*(1-CH);
                etapx1=-k*SH*etax0 +   CH*etapx0 + rhoinv/k*SH;
            }
            else
            {
                // K==0,angle>0
                betax1 = betax0 -2*len*alphax0+ len*len*gammax0;
                alphax1=                    alphax0-len*gammax0;
                gammax1= gammax0;
                nux1 = nux0 + 0.5*INV_PI*(atan(gammax0*len - alphax0) + atan(alphax0));
                etax1 = etax0 + len*etapx0 +0.5*rhoinv*len*len;
                etapx1=             etapx0 +    rhoinv*len;
            }
        }
        
        else
        {
            // values[K1]==0,pure dipole
            if(angle!=0)
            {
                double k = sqrt(K);
                double C=cos(k*len);
                double S=sin(k*len);
                double C2=cos(2*k*len);
                double S2=sin(2*k*len);
                betax1 =   C*C*betax0 -2.0/k*C*S*alphax0+ S*S/(k*k)*gammax0;
                alphax1= k*C*S*betax0 +       C2*alphax0- C*S/k*gammax0;
                gammax1= (k*k)*S*S*betax0 +  2*k*C*S*alphax0+   C*C*gammax0;
                betay1 = betay0 -2*len*alphay0+ len*len*gammay0;
                alphay1=                     alphay0-len*gammay0;
                gammay1= gammay0;
                nux1 = nux0 +0.5*INV_PI*(atan(gammax0*tan(k*len)/k - alphax0) + atan(alphax0));
                nuy1 = nuy0 + 0.5*INV_PI*(atan(gammay0*len - alphay0) + atan(alphay0));
                etax1=    C*etax0 + S/k*etapx0 + rhoinv/k*(1-C);
                etapx1=-k*S*etax0 +   C*etapx0 + rhoinv/k*S;
            }
            else
            {
                betax1 = betax0 -2*len*alphax0+ len*len*gammax0;
                alphax1=                    alphax0-len*gammax0;
                gammax1= gammax0;
                betay1 = betay0 -2*len*alphay0+ len*len*gammay0;
                alphay1=                     alphay0-len*gammay0;
                gammay1= gammay0;
                nux1 = nux0 + 0.5*INV_PI*(atan(gammax0*len - alphax0) + atan(alphax0));
                nuy1 = nuy0 + 0.5*INV_PI*(atan(gammay0*len - alphay0) + atan(alphay0));
                etax1 = etax0 + len*etapx0 ;
                etapx1=             etapx0 ;

            }
        }
        
        betax0=betax1;
        alphax0=alphax1;
        gammax0=gammax1;
        betay0=betay1;
        alphay0=alphay1;
        gammay0=gammay1;
        nux0=nux1;
        nuy0=nuy1;
        etax0=etax0;
        etapx0=etapx0;

        if((1.0==len_rate) && (0!=angle*values[E2]) )
        {
            double edge2=angle*values[E2];
            double tge2=tan(edge2);
            betax0 = twsin[BETAX];
            alphax0= -rhoinv*tge2*twsin[BETAX]+twsin[ALPHAX];
            gammax0=rhoinv*rhoinv*tge2*tge2*twsin[BETAX] -2.0*rhoinv*tge2*twsin[ALPHAX] + twsin[GAMMAX];
            betay0 = twsin[BETAY];
            alphay0= rhoinv*tge2*twsin[BETAY]+twsin[ALPHAY];
            gammay0=rhoinv*rhoinv*tge2*tge2*twsin[BETAY] +2.0*rhoinv*tge2*twsin[ALPHAY] + twsin[GAMMAY];

            nux1=nux0;
            nuy1=nuy0;
            etax1=etax0;
            etapx1=rhoinv*tge2*etax0+etapx0;

            betax0=betax1;
            alphax0=alphax1;
            gammax0=gammax1;
            betay0=betay1;
            alphay0=alphay1;
            gammay0=gammay1;
            nux0=nux1;
            nuy0=nuy1;
            etax0=etax0;
            etapx0=etapx0;
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
    }

};


#endif