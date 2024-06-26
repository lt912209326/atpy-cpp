#ifndef _CPPRADIATIONINTEGRALS_CPP_
#define _CPPRADIATIONINTEGRALS_CPP_
#include "cppradiationintegrals.h"
int computeRadiationIntegrals(double rhoinv,double blen,double k1,double edge1, double edge2, double betxi,double alfxi,double dxi,double dpxi,double* glb){

    // '''
    //     !----------------------------------------------------------------------*
    //     !     Purpose:                                                         *
    //     !     Calculate synchrotron radiation integrals contribution of        *
    //     !     single element with parameters passed as input.                  *
    //     !                                                                      *
    //     !     Input:                                                           *
    //     !     rhoinv (double) inverse radius of curvature                      *
    //     !     blen (double) length of element                                  *
    //     !     k1 (double) gradient of element                                  *
    //     !     e1, e2 (double) pole face rotations at entrance and exit         *
    //     !     betxi, alfxi, dxi, dpxi (double) twiss parameters in x plane     *
    //     !     Output:                                                          *
    //     !     I[8] synchrotron radiation integral                              *
    //     !                                                                      *
    //     !     Author: Ghislain Roy - June 2014                                 *
    //     !----------------------------------------------------------------------*
    // '''

    // # local variables;
    double  dx2, gamx, dispaverage, curlyhaverage, lq;
    double  betx, alfx, dx, dpx, u0x, u1x, u2x;
    double  gammai, betxaverage, k1n;
    double k2, k, klen;
    double zero, one, two;
    double Cx,Cy,Sx,Sy;
    Cx=Cy=Sx=Sy=0;

    betx = betxi;
    dx = dxi;
    // # k1 = 0;
    zero, one, two = 0.0, 1.0, 2.0;

    alfx = alfxi - betxi*rhoinv*tan(edge1);
    dpx = dpxi + dxi*rhoinv*tan(edge1);

    gamx = (1+alfx*alfx)/betx;

    // # effect of poleface rotation;

    // # global gradient combining weak focusing and dipole gradient;
    // # k2 can be positive or negative and k can be real or imaginary;
    // # k2 = rhoinv*rhoinv + k1;
    k2 = rhoinv*rhoinv + k1;
    if(k2>0){

        k = sqrt(k2);
        Cx=cos(k*blen);
        Sx=sin(k*blen)/k;
    }
    else if(k2<0){

        k = sqrt(-k2);
        Cx=cosh(k*blen);
        Sx=sinh(k*blen)/k;
    }
    else{
        k=0;
        Cx=1;
        Sx=blen;
    }
    // # k = sqrt(k2);
    klen = k*blen;

    // # propagation of dispersion at exit;
    dx2 = dx*Cx + dpx*Sx + rhoinv*(1-Cx)/k2;

    dispaverage = dx*Sx/blen + dpx*(1-Cx)/(k2*blen)+ rhoinv*(blen - Sx)/(k2*blen);

    curlyhaverage =  gamx*dx*dx + 2*alfx*dx*dpx + betx*dpx*dpx 
                + 2*rhoinv*blen*( -(gamx*dx + alfx*dpx)*(blen-Sx)/(k2*blen*blen)  + (alfx*dx + betx*dpx)*(1-Cx)/(k2*blen*blen)) 
                + rhoinv*rhoinv/(k2*k2*blen)*( 
                       0.5*gamx*(3*blen - 4*Sx + Sx*Cx) 
                     - alfx*pow(1-Cx,2) 
                     + 0.5*k2*betx*(blen-Cx*Sx));

    if (fabs(rhoinv) > 1e-10){
        glb[RI1 ] += dispaverage * rhoinv * blen;
        glb[RI2 ] += rhoinv*rhoinv * blen;
        glb[RI3 ] += pow(fabs(rhoinv),3) * blen;
        glb[RI3A] += pow(rhoinv,3)*blen;
        glb[RI4 ] += dispaverage*rhoinv*(pow(rhoinv,2) + 2*k1)* blen- rhoinv*rhoinv*(dx*tan(edge1) + dx2*tan(edge2));
        glb[RI5 ] += curlyhaverage * pow(fabs(rhoinv),3) * blen;
        // cout<<"computeRadiationIntegrals: "<<glb[RI5]<<endl;
    }
    return 0;


}


#endif