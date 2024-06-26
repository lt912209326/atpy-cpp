#ifndef _CPPDRIVINGTERMS3_CPP_
#define _CPPDRIVINGTERMS3_CPP_


#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <exception>
#include "cppconstants.h"
#include "cppcomponent.h"
#include "cppdrivingterms.h"



// 只适合周期结构中对称右半部分、且仅仅dnux_dJx, dnux_nJy, dnuy_dJx完成，其他的还需继续

void computmdrivingTerms_fast_period(vector<CppComponent*>& bline, vector<MultipoleData>& md, const vector<size_t>& mult_pos, vector<size_t>& sext_slice_index, const int multipole_slice,int nPeriods, double* GLB, const Status& stat )
{
    /*  Basmd on J. Bengtsson, SLS Note 9/97, March 7, 1997, with corrections per W. Guo (NSLS) */
    /*  Revismd to follow C. X. Wang AOP-TN-2009-020 for second-order terms */
    complex<double>  h20001, h00201, h10010, h10100, h10002;// h11001, h00111,
    complex<double> h21000, h30000, h10110, h10020, h10200;
    // conjugate terms of above
    complex<double>  h02001, h00021, h12000, h03000, h01110, h01200, h01020;
    // second order driving terms
    complex<double> h31000, h40000, h22000, h11110, h00220;
    complex<double> h20110, h11200, h20020, h20200, h00310, h00400, h11002, h00112;
    

    double h11001, h00111;
    // driving term at every slice
    complex<double>  tmph20001, tmph00201, tmph10002, tmph10010, tmph10100;  //h11001, h00111,
    // conjugate terms of above
    complex<double>  tmph02001, tmph00021, tmph12000, tmph03000, tmph01110, tmph01200, tmph01020;
    complex<double> tmph21000, tmph30000, tmph10110, tmph10020, tmph10200;
    // second order driving terms
    complex<double> tmph31000, tmph40000; //, h11002, h00112, h22000, h11110, h00220, 

    double  h20001_fluc, h00201_fluc, h10002_fluc, h10010_fluc, h10100_fluc, h11001_fluc, h00111_fluc;
    // conjugate terms of above
    double  h02001_fluc, h00021_fluc, h12000_fluc, h03000_fluc, h01110_fluc, h01200_fluc, h01020_fluc;
    double  h21000_fluc, h30000_fluc, h10110_fluc, h10020_fluc, h10200_fluc;
    // second order driving terms
    // second order driving terms
    double h31000_fluc, h40000_fluc, h22000_fluc, h11110_fluc, h00220_fluc;
    double h20110_fluc, h11200_fluc, h20020_fluc, h20200_fluc, h00310_fluc, h00400_fluc, h11002_fluc, h00112_fluc;


    double dnux_dJx, dnux_dJy, dnuy_dJy, d2Qx, d2Qy, dQx, dQy ,etaxpp,betaxp,betayp;
    complex<double> t1, t2, t3, t4;
    complex<double> ii;
    complex<double> periodicFactor[9][9];

    double tune[2]={bline[bline.size()-1]->tws(-1,NUX), bline[bline.size()-1]->tws(-1,NUY) } ;
    CppElement* elem;

#define PF(i,j) (periodicFactor[4+i][4+j])
    double betax1, betay1, phix1, phiy1, etax1, termSign,gammax,gammay;
    double b1L, b2L, a2L, b3L, b4L, nux, nuy;
    // ELEMENT_LIST *elem;
    double two=2, three=3, four=4;
    // ELEMDATA *md = NULL;
    size_t  iE, jE, i, j, idxi, idxj;
    //double sqrt8, sqrt2;
    double tilt;

    //sqrt8 = sqrt((double)8);
    //sqrt2 = sqrt((double)2);
    ii = complex<double>(0,1);
    b1L = b2L = a2L = b3L = b4L = 0;
    /*  accumulate real and imaginary parts */
    h10100 = h10010 = complex<double>(0,0);
    h11001 = h00111 =0.0;
    h10002 = h20001 = h00201 = h21000 = h30000 = h10110 = h10020 = h10200 = complex<double>(0,0);
    // complex<double>(0,0);      //h11002 = h00112=
    h31000 = h40000 =  h20110 = h11200 = h20020 = h20200 = h00310 = h00400= h11002= h00112 = complex<double>(0,0);
    // conjugate driving terms
    h02001= h00021= h12000= h03000= h01110= h01200= h01020=complex<double>(0,0);
    tmph02001= tmph00021= tmph12000= tmph03000= tmph01110= tmph01200= tmph01020=complex<double>(0,0);
    
    h10100_fluc= h11001_fluc= h10010_fluc= h00111_fluc=0;
    h20001_fluc= h00201_fluc= h10002_fluc= h21000_fluc= h30000_fluc= h10110_fluc= h10020_fluc= h10200_fluc=0;
    h31000_fluc = h40000_fluc =  h20110_fluc = h11200_fluc = h20020_fluc = h20200_fluc = h00310_fluc = 0;
    h00400_fluc= h22000_fluc= h11110_fluc= h00220_fluc=h11002_fluc= h00112_fluc = 0;
    

    if(stat.printout )cout<<"cppdrivingterms3.cpp: "<<nPeriods<<endl;
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


    size_t comp_pos_index=-1, comp_pos=bline[mult_pos[0] ]->position, slice_pos=-1,nslice=bline[mult_pos[0] ]->elem->nslice, sext_slice_cnt=0 ;
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

            b1L=a2L = b2L = b3L = b4L = 0;
            elem=bline[comp_pos]->elem;
            switch (elem->kind) 
            {
            case SEXTUPOLE:
                if( fabs(elem->values[L])<1e-8 ){
                    b3L = 0.5*elem->values[K2];
                }
                else{
                    b3L = 0.5*elem->values[K2] * elem->values[L]/elem->nslice;
                }
                break;
            case DIPOLE:
                b2L = elem->values[K1] * elem->values[L]/elem->nslice;
                break;
            case QUADRUPOLE:
                b2L = elem->values[K1] * elem->values[L]/elem->nslice;
                break;
            case OCTUPOLE:
                if( fabs(elem->values[L])<1e-8 ){
                    b4L = elem->values[K3]/6.0;
                }
                else{
                    b4L = elem->values[K3] * elem->values[L]/elem->nslice/6;
                }
                break;
            case EXACTDRIFT:
                b1L=elem->values[L]/elem->nslice;
                break;
            default:
                break;
            }      
        }
        if(SEXTUPOLE==elem->kind){
            sext_slice_index[sext_slice_cnt]=i;
            sext_slice_cnt+=1;
        }
        
        if (b1L || a2L || b2L || b3L || b4L) {
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
            gammax = bline[comp_pos]->tws(slice_pos,GAMMAX);
            gammay = bline[comp_pos]->tws(slice_pos,GAMMAY);

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

            // reset local drving terms variations
            tmph20001 = tmph00201 = tmph21000 = tmph30000 = tmph10110 = tmph10020 = tmph10200 = complex<double>(0,0);
            // conjugate terms
            tmph02001= tmph00021= tmph12000= tmph03000= tmph01110= tmph01200= tmph01020=complex<double>(0,0);

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
                tmph20001 = (b3L*betax1*etax1/2-b2L*betax1/4)/2*md[i].px[2]*PF(2,0);
                tmph00201 = (b2L*betay1/4-b3L*betay1*etax1/2)/2*md[i].py[2]*PF(0,2);

                tmph02001 = std::conj(tmph20001);
                tmph00021 = std::conj(tmph00201);

                /*  h10002 */
                tmph10002 = (b3L*md[i].rbetax*sqr(etax1)-b2L*md[i].rbetax*etax1)/2*md[i].px[1]*PF(1,0);
            }
            if (md[i].b3L ){
                /*  first-order geometric terms from sextupoles */
                /*  h21000 */
                tmph21000 = b3L*md[i].rbetax*betax1/8*md[i].px[1]*PF(1,0);

                    /*  h30000 */
                tmph30000 = b3L*md[i].rbetax*betax1/24*md[i].px[3]*PF(3,0);
                    
                    /*  h10110 */
                tmph10110 = -b3L*md[i].rbetax*betay1/4*md[i].px[1]*PF(1,0);

                    /*  h10020 and h10200 */
                tmph10020 = -b3L*md[i].rbetax*betay1/8*md[i].px[1]*conj(md[i].py[2])*PF(1,-2);

                tmph10200 = -b3L*md[i].rbetax*betay1/8*md[i].px[1]*md[i].py[2]*PF(1,2);

                
                tmph12000 = std::conj(tmph21000);
                tmph03000 = std::conj(tmph30000);
                tmph01110 = std::conj(tmph10110);
                tmph01200 = std::conj(tmph10020);
                tmph01020 = std::conj(tmph10200);
            }
            if (b4L) {
                /*  second-order terms from leading order effects of octupoles */
            /*  Ignoring a large number of terms that are not also driven by sextupoles */
                // cout<<"computmdrivingTerms: "<<endl;
                h22000 += 3*b4L*md[i].betax2/32*nPeriods;
                h11110 -= 3*b4L*betax1*betay1/8*nPeriods;
                h00220 += 3*b4L*md[i].betay2/32*nPeriods;
            
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
            if(b1L ){
                h22000 -= 3*b1L*gammax*gammax/64*nPeriods;
                h11110 -= b1L*gammax*gammay/16*nPeriods;
                h00220 -= 3*b1L*gammay*gammay/64*nPeriods;
            }
            // i++;

                
            h22000 +=  ii*( 3.0*(h21000*tmph12000 - h12000*tmph21000) 
                          + 9.0*(h30000*tmph03000 - h03000*tmph30000) );
            h11110 +=  ii*(   2.0*(h21000*tmph01110 - h01110*tmph21000  ) 
                            - 2.0*(h12000*tmph10110 - h10110*tmph12000 ) 
                            - 4.0*(h10020*tmph01200 - h01200*tmph10020 ) 
                            + 4.0*(h10200*tmph01020 - h01020*tmph10200 ) );
            h00220 +=  ii*(   (h10020*tmph01200 - h01200*tmph10020) 
                            + (h10110*tmph01110 - h01110*tmph10110) 
                            + (h10200*tmph01020 - h01020*tmph10200));

            h31000 +=  ii*( 6.0*(h30000*tmph12000 - h12000*tmph30000)      ) ;
            h40000 +=  ii*( 3.0*(h30000*tmph21000 - h21000*tmph30000)     ) ;
            h20110 +=  ii*( 3.0*(h30000*tmph01110 - h01110*tmph30000) 
                              - (h21000*tmph10110 - h10110*tmph21000) 
                              + 4.0*(h10200*tmph10020 - h10020*tmph10200)   ) ;
            h11200 +=  ii*( 2.0*(h10200*tmph12000 - h12000*tmph10200) 
                            + 2.0*(h21000*tmph01200 - h01200*tmph21000) 
                            + 2.0*(h10200*tmph01110 - h01110*tmph10200) 
                            -  2.0*(h10110*tmph01200 - h01200*tmph10110)      ) ;
            h20020 +=  ii*( -(h21000*tmph10020 - h10020*tmph21000) 
                        + 3.0*(h30000*tmph01020 - h01020*tmph30000) 
                        + 2.0*(h10110*tmph10020 - h10020*tmph10110)  ) ;
            h20200 +=  ii*( 3.0*(h30000*tmph01200 - h01200*tmph30000) 
                            + (h10200*tmph21000 - h21000*tmph10200) 
                            - 2.0*(h10110*tmph10200 - h10200*tmph10110)   ) ;
            h00310 +=  ii*( (h10200*tmph01110 - h01110*tmph10200) 
                            + (h10110*tmph01200 - h01200*tmph10110)    ) ;
            h00400 +=  ii*( (h10200*tmph01200 - h01200*tmph10200)     ) ;
            // + 9.0*h30000*tmph03000

            // h22000 +=  -1.0/64*ii*( 3.0*h21000*tmph12000 + h30000*tmph03000 );
            // h11110 +=  1.0/16*ii*( 2.0*h21000*tmph01110 + h10020*tmph01200 +  h10200*tmph01020);
            // h00220 +=  -1.0/64*ii*( 4.0*h10110*tmph01110 + h10020*tmph01200 + h10200*tmph01020 );

            // h31000 +=  1.0/64*ii*( 2.0*h30000*tmph12000   ) ;
            // h40000 +=  1.0/64*ii*( h30000*tmph12000     ) ;
            // h20110 +=  1.0/64*ii*( h30000*tmph01110 + h21000*tmph10110 + 2.0*h10200*tmph12000   ) ;
            // h11200 +=  2.0*ii*( h10200*tmph12000 + h21000*tmph01200 + h10200*tmph01110 +  h10110*tmph01200      ) ;
            // h20020 +=  ii*( h21000*tmph10020 + h30000*tmph01020 + 4.0*h10110*tmph10020  ) ;
            // h20200 +=  ii*( h30000*tmph01200 + h10200*tmph21000 + 4.0*h10110*tmph10200   ) ;
            // h00310 +=  1.0/32*ii*( h10200*tmph01110 + h10110*tmph01200    ) ;
            // h00400 +=  1.0/64*ii*( h10200*tmph01200     ) ;

            // h10100 = h10010 ;
            // h11001 = h00111 =0.0;
            // not used in second order terms
            h10002 += tmph10002;
            h20001 += tmph20001;
            h00201 += tmph00201;
            //used in second order terms
            h21000 += tmph21000;
            h30000 += tmph30000;
            h10110 += tmph10110;
            h10020 += tmph10020;
            h10200 += tmph10200;
            // conjugate driving terms
            // h02001= h00021= 
            h12000 += tmph12000;
            //not used in second order terms
            h03000 += tmph03000;
            h01110 += tmph01110;
            h01200 += tmph01200;
            h01020 += tmph01020;

        }
        else{
            md[i].s = bline[comp_pos]->tws(slice_pos,COORD);
            md[i].b2L = b2L;
            md[i].b3L = b3L;
        }

        
            // h22000 +=  ii*( 3.0*h21000*tmph12000 + 9.0*h30000*tmph03000 );
            // h11110 +=  ii*( 2.0*h21000*tmph01110 - 2.0*h12000*tmph10110 - 4.0*h10020*tmph01200 +  4.0*h10200*tmph01020);
            // h00220 +=  ii*(   h10020*tmph01200 +   h10110*tmph01110);

            // h31000 +=  ii*( 3.0*h30000*tmph12000 + 9.0*h30000*tmph03000     ) ;
            // h40000 +=  ii*( 6.0*h30000*tmph12000     ) ;
            // h20110 +=  ii*( 3.0*h30000*tmph01110 - h21000*tmph10110 + 4.0*h10200*tmph12000   ) ;
            // h11200 +=  ii*( 2.0*h10200*tmph12000 + 2.0*h21000*tmph01200 + 2.0*h10200*tmph01110 -  2.0*h10110*tmph01200      ) ;
            // h20020 +=  ii*( -h21000*tmph10020 + 3.0*h30000*tmph01020 + 2.0*h10110*tmph10020  ) ;
            // h20200 +=  ii*( 3.0*h30000*tmph01200 + h10200*tmph21000 - 2.0*h10110*tmph10200   ) ;
            // h00310 +=  ii*( h10200*tmph01110 + h10110*tmph01200    ) ;
            // h00400 +=  ii*( h10200*tmph01200     ) ;

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
        // 2nd order rdts
        bline[comp_pos]->tws(slice_pos,LH31000) = abs(h31000) ;
        bline[comp_pos]->tws(slice_pos,LH40000) = abs(h40000) ;
        bline[comp_pos]->tws(slice_pos,LH22000) = abs(h22000) ;
        bline[comp_pos]->tws(slice_pos,LH11110) = abs(h11110) ;
        bline[comp_pos]->tws(slice_pos,LH00220) = abs(h00220) ;
        bline[comp_pos]->tws(slice_pos,LH20110) = abs(h20110) ;
        bline[comp_pos]->tws(slice_pos,LH11200) = abs(h11200) ;
        bline[comp_pos]->tws(slice_pos,LH20020) = abs(h20020) ;
        bline[comp_pos]->tws(slice_pos,LH20200) = abs(h20200) ;
        bline[comp_pos]->tws(slice_pos,LH00310) = abs(h00310) ;
        bline[comp_pos]->tws(slice_pos,LH00400) = abs(h00400) ;
        // bline[comp_pos]->tws(slice_pos,LH11002) = abs(h11002) ;
        // bline[comp_pos]->tws(slice_pos,LH00112) = abs(h00112) ;
            
        slice_pos+=1;
        if(slice_pos == nslice){
            // calculate the sum of fluctuation
            h11001_fluc += abs(h11001/PI);
            h00111_fluc += abs(h00111/PI);
            h20001_fluc += abs(h20001);
            h00201_fluc += abs(h00201);
            h10002_fluc += abs(h10002);

            h10100_fluc += abs(h10100);
            h10010_fluc += abs(h10010);

            h21000_fluc += abs(h21000);
            h30000_fluc += abs(h30000);
            h10110_fluc += abs(h10110);
            h10020_fluc += abs(h10020);
            h10200_fluc += abs(h10200);
                    
            h31000_fluc += abs(h31000) ;
            h40000_fluc += abs(h40000) ;
            h22000_fluc += abs(h22000) ;
            h11110_fluc += abs(h11110) ;
            h00220_fluc += abs(h00220) ;
            h20110_fluc += abs(h20110) ;
            h11200_fluc += abs(h11200) ;
            h20020_fluc += abs(h20020) ;
            h20200_fluc += abs(h20200) ;
            h00310_fluc += abs(h00310) ;
            h00400_fluc += abs(h00400) ;
            // h11002_fluc += abs(h11002) ;
            // h00112_fluc += abs(h00112) ;

        }
    }

    if(stat.printout )cout<<"cppdrivingterms3.cpp: "<<2<<endl;
    /*  Doi with the leading-order quad and sext terms */
    if(stat.rdt_fluctuation){
        // calculate the mean of fluctuation
        
        GLB[H11001] = h11001_fluc/comp_pos_index ;
        GLB[H00111] = h00111_fluc/comp_pos_index ;
        GLB[H20001] = h20001_fluc/comp_pos_index ;
        GLB[H00201] = h00201_fluc/comp_pos_index ;
        GLB[H10002] = h10002_fluc/comp_pos_index ;

        GLB[H10100] = h10100_fluc/comp_pos_index ;
        GLB[H10010] = h10010_fluc/comp_pos_index ;

        GLB[H21000] = h21000_fluc/comp_pos_index ;
        GLB[H30000] = h30000_fluc/comp_pos_index ;
        GLB[H10110] = h10110_fluc/comp_pos_index ;
        GLB[H10020] = h10020_fluc/comp_pos_index ;
        GLB[H10200] = h10200_fluc/comp_pos_index ;
        // 2nd order rdts
        GLB[H22000] = h22000_fluc/comp_pos_index ;
        GLB[H11110] = h11110_fluc/comp_pos_index ;
        GLB[H00220] = h00220_fluc/comp_pos_index ;
        GLB[H31000] = h31000_fluc/comp_pos_index ;
        GLB[H40000] = h40000_fluc/comp_pos_index ;
        GLB[H20110] = h20110_fluc/comp_pos_index ;
        GLB[H11200] = h11200_fluc/comp_pos_index ;
        GLB[H20020] = h20020_fluc/comp_pos_index ;
        GLB[H20200] = h20200_fluc/comp_pos_index ;
        GLB[H00310] = h00310_fluc/comp_pos_index ;
        GLB[H00400] = h00400_fluc/comp_pos_index ;
    }
    else{
        
        GLB[H11001] = h11001/PI;
        GLB[H00111] = h00111/PI;
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
        // 2nd order rdts
        GLB[H22000] = abs(h22000);
        GLB[H11110] = abs(h11110);
        GLB[H00220] = abs(h00220);
        GLB[H31000] = abs(h31000);
        GLB[H40000] = abs(h40000);
        GLB[H20110] = abs(h20110);
        GLB[H11200] = abs(h11200);
        GLB[H20020] = abs(h20020);
        GLB[H20200] = abs(h20200);
        GLB[H00310] = abs(h00310);
        GLB[H00400] = abs(h00400);
    }
    
    size_t bline_length=bline.size();
    size_t num_local_RDTs=TWS_NUM-LH11001;
    size_t elem_kind=0;
    for(size_t i=1;i<bline_length;i++){
        elem_kind=bline[i]->elem->kind;
        if( EXACTDRIFT ==elem_kind ||  (stat.combineddipole && DIPOLE==elem_kind)  || QUADRUPOLE==elem_kind || SEXTUPOLE==elem_kind || OCTUPOLE==elem_kind) continue;
        memcpy(&bline[i]->tws(-1,LH11001),&bline[i-1]->tws(-1,LH11001), num_local_RDTs*__SIZEOF_DOUBLE__ );
    }

    // GLB[D2QX] = d2Qx*nPeriods;
    // GLB[D2QY] = d2Qy*nPeriods;
    // GLB[D2QX] = h11002.real();
    // GLB[D2QY] = h00112.real();


    GLB[DNUX_DJX]= -2.0*INV_PI*h22000.real()*nPeriods;
    GLB[DNUX_DJY]=     -INV_PI*h11110.real()*nPeriods;
    GLB[DNUY_DJY]= -2.0*INV_PI*h00220.real()*nPeriods;


}

#endif