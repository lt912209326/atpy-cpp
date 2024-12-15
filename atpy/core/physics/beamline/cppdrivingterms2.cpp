#ifndef _CPPDRIVINGTERMS2_CPP_
#define _CPPDRIVINGTERMS2_CPP_


#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <exception>
#include "cppconstants.h"
#include "cppcomponent.h"
#include "cppdrivingterms.h"




// 只适合周期结构中对称右半部分、且仅仅dnux_dJx, dnux_nJy, dnuy_dJx完成，其他的还需继续

void computmdrivingTerms_nonperiod(vector<CppComponent*>& bline, vector<MultipoleData>& md, const vector<size_t>& mult_pos, vector<size_t>& sext_slice_index, const int multipole_slice,int nPeriods, double* GLB, const Status& stat )
{
    /*  Basmd on J. Bengtsson, SLS Note 9/97, March 7, 1997, with corrections per W. Guo (NSLS) */
    /*  Revismd to follow C. X. Wang AOP-TN-2009-020 for second-order terms */
    complex<double>  h20001, h00201, h10002, h10010, h10100;  //h11001, h00111,
    complex<double> h21000, h30000, h10110, h10020, h10200;
    complex<double> h31000, h40000; //, h11002, h00112, h22000, h11110, h00220, 
    complex<double> h20110, h11200, h20020, h20200, h00310, h00400, h11002, h00112;
    double dnux_dJx, dnux_dJy, dnuy_dJy, d2Qx, d2Qy, dQx, dQy ,etaxpp,betaxp,betayp;
    double h11001, h00111;
    complex<double> t1, t2, t3, t4;
    complex<double> ii;
    complex<double> periodicFactor[9][9];

    double tune[2]={bline[bline.size()-1]->tws(-1,NUX), bline[bline.size()-1]->tws(-1,NUY) } ;
    CppElement* elem;

#define PF(i,j) (periodicFactor[4+i][4+j])
    double betax1, betay1, phix1, phiy1, etax1, termSign;
    double b2L, a2L, b3L, b4L, nux, nuy;
    // ELEMENT_LIST *elem;
    double two=2, three=3, four=4;
    // ELEMDATA *md = NULL;
    size_t  iE, jE, i, j, idxi, idxj;
    //double sqrt8, sqrt2;
    double tilt;

    //sqrt8 = sqrt((double)8);
    //sqrt2 = sqrt((double)2);
    ii = complex<double>(0,1);

    /*  accumulate real and imaginary parts */
    h11001 = h00111 =0.0;
    h20001 = h00201 = h10002 = complex<double>(0,0);
    h21000 = h30000 = h10110 = h10020 = h10200 = complex<double>(0,0);
    h31000 = h40000 =  complex<double>(0,0);      //h11002 = h00112=
    h20110 = h11200 = h20020 = h20200 = h00310 = h00400= h11002= h00112 = complex<double>(0,0);
    h10100 = h10010 = complex<double>(0,0);

    if(stat.printout )cout<<nPeriods<<endl;
    dnux_dJx = dnux_dJy = dnuy_dJy = 0;
    if (nPeriods!=1) {
        double a1, a2;
        for (i=0; i<9 ;i++) {
            for (j=0; j<9; j++) {
                a1 = PIx2*(tune[0]*(i-4)+tune[1]*(j-4));
                a2 = a1/nPeriods;
                periodicFactor[i][j] = (exp(ii*a1)-1.0)/(exp(ii*a2)-1.0);
            }
        }
    } else {
        for (i=0 ;i<9 ;i++)
            for (j=0 ;j<9 ;j++)
                periodicFactor[i][j] = 1;
    }


    int comp_pos_index=-1, comp_pos=bline[mult_pos[0] ]->position, slice_pos=-1,nslice=bline[mult_pos[0] ]->elem->nslice, sext_slice_cnt=0 ;
    for(size_t i=0;i<multipole_slice;i++)
    {
        if(!(slice_pos<nslice)|| comp_pos_index<0 ){
            // 切片计数大于等于元件切片数或开始时，组件元素更换
            comp_pos_index+=1;
            comp_pos = bline[mult_pos[comp_pos_index ] ]->position;
            nslice=bline[mult_pos[comp_pos_index ] ]->elem->nslice;
            slice_pos=0;
            if(stat.nonlineartermonly){
                bline[comp_pos-1]->tws(-1,CHROMX)=dQx;
                bline[comp_pos-1]->tws(-1,CHROMY)=dQy;
                bline[comp_pos]->linearoptics(&(bline[comp_pos-1]->tws[-1] ) ,&stat, GLB);
                dQx+=bline[comp_pos]->tws(-1,DCHROMX);
                dQy+=bline[comp_pos]->tws(-1,DCHROMY);
            }

            a2L = b2L = b3L = b4L = 0;
            elem=bline[comp_pos]->elem;
            switch (elem->kind) 
            {
            case SEXTUPOLE:
                b3L = elem->values[K2] * elem->values[L]/elem->nslice/2;
                break;
            case DIPOLE:
                b2L = elem->values[K1] * elem->values[L]/elem->nslice;
                break;
            case QUADRUPOLE:
                b2L = elem->values[K1] * elem->values[L]/elem->nslice;
                break;
            case OCTUPOLE:
                b4L = elem->values[K3] * elem->values[L]/elem->nslice/6;
                break;
            default:
                break;
            }      
        }
        if(SEXTUPOLE==elem->kind){
            sext_slice_index[sext_slice_cnt]=i;
            sext_slice_cnt+=1;
        }
        
        if (a2L || b2L || b3L || b4L) {
            if (i>0 && slice_pos>0) {
                betax1 = (bline[comp_pos]->tws(slice_pos,BETAX)  + bline[comp_pos]->tws(slice_pos-1,BETAX) )/2;
                etax1  = (bline[comp_pos]->tws(slice_pos,ETAX)  + bline[comp_pos]->tws(slice_pos-1,ETAX) )/2;
                phix1  = PIx2*(bline[comp_pos]->tws(slice_pos,NUX)  + bline[comp_pos]->tws(slice_pos-1,NUX) )/2;
                betay1 = (bline[comp_pos]->tws(slice_pos,BETAY)  + bline[comp_pos]->tws(slice_pos-1,BETAY) )/2;
                phiy1  = PIx2*(bline[comp_pos]->tws(slice_pos,NUY)  + bline[comp_pos]->tws(slice_pos-1,NUY) )/2;
            } 
            else if(i>0 && 0==slice_pos){
                betax1 = (bline[comp_pos]->tws(slice_pos,BETAX)  + bline[comp_pos-1]->tws(-1,BETAX) )/2;
                etax1  = (bline[comp_pos]->tws(slice_pos,ETAX)  + bline[comp_pos-1]->tws(-1,ETAX) )/2;
                phix1  = PIx2*(bline[comp_pos]->tws(slice_pos,NUX)  + bline[comp_pos-1]->tws(-1,NUX) )/2;
                betay1 = (bline[comp_pos]->tws(slice_pos,BETAY)  + bline[comp_pos-1]->tws(-1,BETAY) )/2;
                phiy1  = PIx2*(bline[comp_pos]->tws(slice_pos,NUY)  + bline[comp_pos-1]->tws(-1,NUY) )/2;
            }
            else {
                betax1 =  bline[comp_pos]->tws(slice_pos,BETAX);
                etax1  = bline[comp_pos]->tws(slice_pos,ETAX);
                phix1  = PIx2*bline[comp_pos]->tws(slice_pos,NUX);
                betay1 =  bline[comp_pos]->tws(slice_pos,BETAY);
                phiy1  = PIx2*bline[comp_pos]->tws(slice_pos,NUY);
            }

            md[i].s = bline[comp_pos]->tws(slice_pos,COORD);
            md[i].b2L = b2L;
            md[i].b3L = b3L;
            md[i].etax = etax1;
            md[i].betax = betax1;
            md[i].betax2 = sqr(betax1);
            md[i].rbetax = sqrt(betax1);
            md[i].phix = phix1;
            md[i].px[1] = exp(ii*phix1);
            md[i].px[2] = md[i].px[1]*md[i].px[1];
            md[i].px[3] = md[i].px[1]*md[i].px[2];
            md[i].px[4] = md[i].px[1]*md[i].px[3];
            md[i].betay = betay1;
            md[i].betay2 = sqr(betay1);
            md[i].rbetay = sqrt(betay1);
            md[i].phiy = phiy1;
            md[i].py[1] = exp(ii*phiy1);
            md[i].py[2] = md[i].py[1]*md[i].py[1];
            md[i].py[3] = md[i].py[1]*md[i].py[2];
            md[i].py[4] = md[i].py[1]*md[i].py[3];

            if (a2L) {
                /*  linear coupling terms */
                h10010 += (a2L/4)*md[i].rbetax*md[i].rbetay*md[i].px[1]/md[i].py[1]*PF(1, -1);
                h10100 += (a2L/4)*md[i].rbetax*md[i].rbetay*md[i].px[1]*md[i].py[1]*PF(1, 1);
            }
            if (b2L || b3L) {
                /*  first-order chromatic terms */
                /*  h11001 and h00111 */
                h11001 += (b3L*betax1*etax1/2-b2L*betax1/4)*nPeriods;
                h00111 += (b2L*betay1/4-b3L*betay1*etax1/2)*nPeriods;
                
                /*  h20001, h00201 */
                h20001 += (b3L*betax1*etax1/2-b2L*betax1/4)/2*md[i].px[2]*PF(2,0);
                h00201 += (b2L*betay1/4-b3L*betay1*etax1/2)/2*md[i].py[2]*PF(0,2);

                /*  h10002 */
                h10002 += (b3L*md[i].rbetax*sqr(etax1)-b2L*md[i].rbetax*etax1)/2*md[i].px[1]*PF(1,0);
            }
            if (md[i].b3L ){
                /*  first-order geometric terms from sextupoles */
                /*  h21000 */
                h21000 += b3L*md[i].rbetax*betax1/8*md[i].px[1]*PF(1,0);

                    /*  h30000 */
                h30000 += b3L*md[i].rbetax*betax1/24*md[i].px[3]*PF(3,0);
                    
                    /*  h10110 */
                h10110 += -b3L*md[i].rbetax*betay1/4*md[i].px[1]*PF(1,0);

                    /*  h10020 and h10200 */
                h10020 += -b3L*md[i].rbetax*betay1/8*md[i].px[1]*conj(md[i].py[2])*PF(1,-2);

                h10200 += -b3L*md[i].rbetax*betay1/8*md[i].px[1]*md[i].py[2]*PF(1,2);
            }
            if (b4L) {
                /*  second-order terms from leading order effects of octupoles */
            /*  Ignoring a large number of terms that are not also driven by sextupoles */
                // cout<<"computmdrivingTerms: "<<endl;
                dnux_dJx += 3*b4L*md[i].betax2/(8*PI)*nPeriods;
                dnux_dJy -= 3*b4L*betax1*betay1/(4*PI)*nPeriods;
                dnuy_dJy += 3*b4L*md[i].betay2/(8*PI)*nPeriods;
            
                // h22000 += 3*b4L*md[i].betax2/32*nPeriods;
                // h11110 += -3*b4L*betax1*betay1/8*nPeriods;
                // h00220 += 3*b4L*md[i].betay2/32*nPeriods;
                h31000 += b4L*md[i].betax2/16*md[i].px[2]*PF(2,0);
                h40000 += b4L*md[i].betax2/64*md[i].px[4]*PF(4,0);
                h20110 += -3*b4L*betax1*betay1/16*md[i].px[2]*PF(2,0);
                h11200 += -3*b4L*betax1*betay1/16*md[i].py[2]*PF(0,2);
                h20020 += -3*b4L*betax1*betay1/32*md[i].px[2]*conj(md[i].py[2])*PF(2,-2);
                h20200 += -3*b4L*betax1*betay1/32*md[i].px[2]*md[i].py[2]*PF(2,2);
                h00310 += b4L*md[i].betay2/16*md[i].py[2]*PF(0,2);
                h00400 += b4L*md[i].betay2/64*md[i].py[4]*PF(0,4);
            }
            // i++;
        }
        else{
            md[i].s = bline[comp_pos]->tws(slice_pos,COORD);
            md[i].b2L = b2L;
            md[i].b3L = b3L;
        }
        // LOC[i][LH11001]=h11001.real();
        // LOC[i][LH00111]=h00111.real();
        // LOC[i][LH11001]=LOC[i-1][LH11001]+LOC[i][CHROMX];
        // LOC[i][LH00111]=LOC[i-1][LH00111]+LOC[i][CHROMY];
        bline[comp_pos]->tws(slice_pos,LH11001)=  h11001/PI;
        bline[comp_pos]->tws(slice_pos,LH00111)=  h00111/PI;
        bline[comp_pos]->tws(slice_pos,LH20001)=abs(h20001);
        bline[comp_pos]->tws(slice_pos,LH00201)=abs(h00201);
        bline[comp_pos]->tws(slice_pos,LH10110)=abs(h10110);
        bline[comp_pos]->tws(slice_pos,LH21000)=abs(h21000);
        bline[comp_pos]->tws(slice_pos,LH30000)=abs(h30000);
        bline[comp_pos]->tws(slice_pos,LH10200)=abs(h10200);
        bline[comp_pos]->tws(slice_pos,LH10002)=abs(h10002);
        bline[comp_pos]->tws(slice_pos,LH10020)=abs(h10020);
            // elem = elem->succ;
        slice_pos+=1;
    }

    /*  Doi with the leading-order quad and sext terms */
    GLB[H11001] = h11001/PI;
    GLB[H00111] = h00111/PI;
    // if(stat.nonlineartermonly){
    //     GLB[H11001] = dQx*nPeriods;
    //     GLB[H00111] = dQy*nPeriods;
    // }
    GLB[H20001] = abs(h20001);
    GLB[H00201] = abs(h00201);
    GLB[H10002] = abs(h10002);

    GLB[H10100] = abs(h10100);
    GLB[H10010] = abs(h10010);

    GLB[H21000] = abs(h21000);
    GLB[H30000] = abs(h30000);
    GLB[H10110] = abs(h10110);
    GLB[H10020] = abs(h10020);
    GLB[H10200] = abs(h10200);
    
    size_t bline_length=bline.size();
    size_t num_local_RDTs=TWS_NUM-LH11001;
    size_t elem_kind=0;
    for(size_t i=1;i<bline_length;i++){
        elem_kind=bline[i]->elem->kind;
        if(  (stat.combineddipole && DIPOLE==elem_kind) || QUADRUPOLE==elem_kind || SEXTUPOLE==elem_kind || OCTUPOLE==elem_kind) continue;
        memcpy(&bline[i]->tws(-1,LH11001),&bline[i-1]->tws(-1,LH11001), num_local_RDTs*__SIZEOF_DOUBLE__ );
    }

    if(!(stat.leaderordertermonly)) {
        /*  compute sextupole contributions to second-order terms */
        if (nPeriods!=1) {    
            throw std::runtime_error("Computating of higher-order driving terms not available when n_periods!=1");
        }

        nux = tune[0];
        nuy = tune[1];
        // h11002-=0.5/PI*h11001;
        // h00112-=0.5/PI*h00111;

        vector<size_t>::iterator iiter,jiter;
        double nux_x2=2*nux, nuy_x2=2*nuy;
        // for (iE=0; iE<multipole_slice; iE++) 
        for(iiter=sext_slice_index.begin();iiter!=sext_slice_index.end();iiter++)
        {
            iE=(*iiter);
            // cout<<"cppdrivingterms: iE: "<<iE<<endl;
            if (md[iE].b3L)
            {
                for (jiter=sext_slice_index.begin(); jiter!=sext_slice_index.end(); jiter++) 
                {
                    jE=(*jiter);
                    termSign = SIGN(md[iE].s - md[jE].s);

                    if (md[jE].b3L) 
                    {
                        dnux_dJx += 2.0*md[iE].b3L*md[jE].b3L/(-16*PI)*pow(md[iE].rbetax*md[jE].rbetax, 3)*
                            (3*cos(fabs(md[iE].phix+md[jE].phix)-PI*nux_x2)/sin(PI*nux_x2) + cos(fabs(3*(md[iE].phix+md[jE].phix))-3*PI*nux_x2)/sin(3*PI*nux_x2));

                        dnux_dJx += 2.0*md[iE].b3L*md[jE].b3L/(-16*PI)*pow(md[iE].rbetax*md[jE].rbetax, 3)*
                            (3*cos(fabs(md[iE].phix-md[jE].phix)-PI*nux_x2)/sin(PI*nux_x2) + cos(fabs(3*(md[iE].phix-md[jE].phix))-3*PI*nux_x2)/sin(3*PI*nux_x2));


                        dnux_dJy += 2.0*md[iE].b3L*md[jE].b3L/(8*PI)*sqrt(md[iE].betax*md[jE].betax)*md[iE].betay*
                            (2*md[jE].betax*cos(fabs(md[iE].phix+md[jE].phix)-PI*nux_x2)/sin(PI*nux_x2) 
                            - md[jE].betay*cos(fabs(md[iE].phix+md[jE].phix)+2*fabs(md[iE].phiy+md[jE].phiy)-PI*(nux_x2+2*nuy_x2))/sin(PI*(nux_x2+2*nuy_x2))
                            + md[jE].betay*cos(fabs(md[iE].phix+md[jE].phix)-2*fabs(md[iE].phiy+md[jE].phiy)-PI*(nux_x2-2*nuy_x2))/sin(PI*(nux_x2-2*nuy_x2)));
                            
                        dnux_dJy += 2.0*md[iE].b3L*md[jE].b3L/(8*PI)*sqrt(md[iE].betax*md[jE].betax)*md[iE].betay*
                            (2*md[jE].betax*cos(fabs(md[iE].phix-md[jE].phix)-PI*nux_x2)/sin(PI*nux_x2) 
                            - md[jE].betay*cos(fabs(md[iE].phix-md[jE].phix)+2*fabs(md[iE].phiy-md[jE].phiy)-PI*(nux_x2+2*nuy_x2))/sin(PI*(nux_x2+2*nuy_x2))
                            + md[jE].betay*cos(fabs(md[iE].phix-md[jE].phix)-2*fabs(md[iE].phiy-md[jE].phiy)-PI*(nux_x2-2*nuy_x2))/sin(PI*(nux_x2-2*nuy_x2)));


                        dnuy_dJy += 2.0*md[iE].b3L*md[jE].b3L/(-16*PI)*sqrt(md[iE].betax*md[jE].betax)*md[iE].betay*md[jE].betay*
                            (4*cos(fabs(md[iE].phix+md[jE].phix)-PI*nux_x2)/sin(PI*nux_x2) 
                            + cos(fabs(md[iE].phix+md[jE].phix)+2*fabs(md[iE].phiy+md[jE].phiy)-PI*(nux_x2+2*nuy_x2))/sin(PI*(nux_x2+2*nuy_x2)) 
                            + cos(fabs(md[iE].phix+md[jE].phix)-2*fabs(md[iE].phiy+md[jE].phiy)-PI*(nux_x2-2*nuy_x2))/sin(PI*(nux_x2-2*nuy_x2)));

                        dnuy_dJy += 2.0*md[iE].b3L*md[jE].b3L/(-16*PI)*sqrt(md[iE].betax*md[jE].betax)*md[iE].betay*md[jE].betay*
                            (4*cos(fabs(md[iE].phix-md[jE].phix)-PI*nux_x2)/sin(PI*nux_x2) 
                            + cos(fabs(md[iE].phix-md[jE].phix)+2*fabs(md[iE].phiy-md[jE].phiy)-PI*(nux_x2+2*nuy_x2))/sin(PI*(nux_x2+2*nuy_x2)) 
                            + cos(fabs(md[iE].phix-md[jE].phix)-2*fabs(md[iE].phiy-md[jE].phiy)-PI*(nux_x2-2*nuy_x2))/sin(PI*(nux_x2-2*nuy_x2)));
                        // termSign = SIGN(md[iE].s - md[jE].s);
                        if (termSign) 
                        {
                            // h11002+=termSign*ii/8*md[iE].rbetax*md[jE].rbetax*md[jE].betax*md[iE].etax*(md[iE].b3L*md[iE].etax-md[iE].b2L)*md[jE].b3L*
                            // (exp(ii*(md[iE].phix-md[jE].phix)) - exp(-ii*(md[iE].phix-md[jE].phix)) );
                            // h00112-=termSign*ii/8*md[iE].rbetax*md[jE].rbetax*md[jE].betay*md[iE].etax*(md[iE].b3L*md[iE].etax-md[iE].b2L)*md[jE].b3L*(
                            //     exp(ii*(md[iE].phix-md[jE].phix)) - exp(-ii*(md[iE].phix-md[jE].phix)) );
                            /*  geometric terms */
                            // h22000 += (1./64)*termSign*ii*md[iE].b3L*md[jE].b3L*
                            //     md[iE].rbetax*md[jE].rbetax*md[iE].betax*md[jE].betax*
                            //     (md[iE].px[3]*conj(md[jE].px[3]) + three*md[iE].px[1]*conj(md[jE].px[1]));
                            h31000 += (1./32)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betax*md[jE].betax*
                                md[iE].px[3]*conj(md[jE].px[1]);
                            t1 = conj(md[iE].px[1])*md[jE].px[1];
                            t2 = md[iE].px[1]*conj(md[jE].px[1]);
                            // h11110 += (1./16)*termSign*ii*md[iE].b3L*md[jE].b3L*
                            //     md[iE].rbetax*md[jE].rbetax*md[iE].betay*
                            //     (md[jE].betax*(t1 - conj(t1) )
                            //     + md[jE].betay*md[iE].py[2]*conj(md[jE].py[2])*(conj(t1) + t1));
                            t1 = exp(-ii*(md[iE].phix-md[jE].phix));
                            t2 = conj(t1);
                            h11200 += (1./32)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betay*exp(ii*(2*md[iE].phiy))*
                                (md[jE].betax*(t1 - t2) + two*md[jE].betay*(t2 + t1) );
                            h40000 += (1./64)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betax*md[jE].betax*
                                md[iE].px[3]*md[jE].px[1];
                            h20020 += (1./64)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betay*
                                (md[jE].betax*conj(md[iE].px[1]*md[iE].py[2])*md[jE].px[3]
                                -(md[jE].betax+four*md[jE].betay)*md[iE].px[1]*md[jE].px[1]*conj(md[iE].py[2]) );
                            h20110 += (1./32)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betay*
                                (md[jE].betax*( conj(md[iE].px[1])*md[jE].px[3]- md[iE].px[1]*md[jE].px[1] ) 
                                 + two*md[jE].betay*md[iE].px[1]*md[jE].px[1]*md[iE].py[2]*conj(md[jE].py[2]) );
                            h20200 += (1./64)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betay*
                                (md[jE].betax*conj(md[iE].px[1])*md[jE].px[3]*md[iE].py[2]-
                                (md[jE].betax-four*md[jE].betay)*md[iE].px[1]*md[jE].px[1]*md[iE].py[2] );
                            // h00220 += (1./64)*termSign*ii*md[iE].b3L*md[jE].b3L*
                            //     md[iE].rbetax*md[jE].rbetax*md[iE].betay*md[jE].betay*
                            //     (md[iE].px[1]*md[iE].py[2]*conj(md[jE].px[1]*md[jE].py[2])
                            //      + four*md[iE].px[1]*conj(md[jE].px[1])
                            //      - conj(md[iE].px[1]*md[jE].py[2])*md[jE].px[1]*md[iE].py[2]
                            //     );
                            h00310 += (1./32)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                md[iE].rbetax*md[jE].rbetax*md[iE].betay*md[jE].betay*md[iE].py[2]*
                                (md[iE].px[1]*conj(md[jE].px[1]) -conj(md[iE].px[1])*md[jE].px[1] );
                            h00400 += (1./64)*termSign*ii*md[iE].b3L*md[jE].b3L*
                                        md[iE].rbetax*md[jE].rbetax*md[iE].betay*md[jE].betay*
                                        md[iE].px[1]*conj(md[jE].px[1])*md[iE].py[2]*md[jE].py[2];

                        }
                    }
                }
            }
        }

    // GLB[H22000] = abs(h22000);
    // GLB[H11110] = abs(h11110);
    // GLB[H00220] = abs(h00220);
    GLB[H31000] = abs(h31000);
    GLB[H40000] = abs(h40000);
    GLB[H20110] = abs(h20110);
    GLB[H11200] = abs(h11200);
    GLB[H20020] = abs(h20020);
    GLB[H20200] = abs(h20200);
    GLB[H00310] = abs(h00310);
    GLB[H00400] = abs(h00400);
    // GLB[D2QX] = d2Qx*nPeriods;
    // GLB[D2QY] = d2Qy*nPeriods;
    // GLB[D2QX] = h11002.real();
    // GLB[D2QY] = h00112.real();

    GLB[DNUX_DJX]=0.5*dnux_dJx*nPeriods;
    GLB[DNUX_DJY]=0.5*dnux_dJy*nPeriods;
    GLB[DNUY_DJY]=0.5*dnuy_dJy*nPeriods;
    }

    if(stat.printout )cout<<"computmdrivingTerms_nonperiod finished"<<endl;
}

#endif