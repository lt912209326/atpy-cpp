#ifndef _CPPSTRUCTURES_H_
#define _CPPSTRUCTURES_H_

#include "cppconstants.h"



struct Status
{
    bool misaligment;                   //misalignment and machine error considered
    bool radiation;                     //synchrotron radiation dumping connsidered
    bool fluctuation;                    //quantum excitation considered
    bool fringe;                        //soft frindge considered
    bool edge;                          //dipole edge will be calculated
    bool lossmap;                       //loss data and information will be calculated
    bool fma;                           //frequent map, tracking data will be produced and saved for fma
    bool slice;                         //slice calculated, only exit of components calculated and less compute source needed if false
    bool combineddipole;                //dipole slice counted when calculate drving terms if true
    bool computedrivingterms;           //driving terms calculated
    bool leaderordertermonly;           //only first-order driving terms calculated, less compute source needed
    bool nonlineartermonly;             //only nonlinear terms calculated such as DTs, first-order chromaticity
    bool linear;                        //only linear element calculated if true
    bool period;                        //a periodic solution is needed
    size_t  npara;                      //thread number for parallel computing
    size_t totalslice;                  //total slice number of the whole beamline
    size_t multipoleslice;              //total slice number of the whole nonlinear elements
    double Trx;                         //the trace of horizontal transfer matrix
    double Try;                         //the trace of vertical transfer matrix
    double dp0;                         //actually dp/p0
    double dp;                          //off-momentum range: (dp0-dp,dp0+dp)
    bool expand;                        //whether enable different component parameters for the same element.
    size_t particle;
    double energy;
    bool second_order_chrom;
    bool third_order_chrom;
    size_t nperiods;
    bool  printout;
    bool transfermatrix;
    double mincouple;
    size_t track_lines;
    size_t track_turns;
    double monitor_dp;                          //monitor to overlook the tunes at dp0 ± monitor_dp
    bool  lazy_compute;
    double larger_monitor_dp;
    bool fast_2nd_order_RDTs;
    double max_betax;
    double max_etax;
    int64_t NP;
    double rf_dp;
    bool rdt_fluctuation;
    bool local_twiss;
    bool off_momentum_rdts;
    double off_rdts_observer;
    uint32_t max_da_range ;
    int64_t chrom_refpt;
    

}; // 添加状态时记得修改 cython文件中 对应状态以及__setitem__ 对应设置

const Status STAT0={ 
                    //  misaligment ,radiation ,fluctuation ,fringe ,edge ,
                     false, false, false, false, false, 
                    //  lossmap ,fma ,slice ,combineddipole ,computedrivingterms ,
                     false, false, true, false, false, 
                    //  leaderordertermonly ,nonlineartermonly ,linear ,period ,npara ,
                     false, false, true, true, 1, 
                    //  totalslice ,multipoleslice ,Trx ,Try ,dp0 ,
                     0, 0, 10, 10, 0, 
                    //  dp ,expand ,particle ,energy ,second_order_chrom ,
                     0, false,ELECTRON,2E9,false,
                    //  third_order_chrom ,nperiods ,printout ,transfermatrix ,mincouple ,
                     false,1,false,true,0.005,
                    //  track_lines ,track_turns ,monitor_dp, larger_monitor_dp ,fast_2nd_order_RDTs
                     13,100,0.01,true,1.0, false,
                    //  max_betax, max_etax, NP, rf_dp, rdt_fluctuation, local_twiss, off_momentum_rdts, off_rdts_observer, max_da_range 
                    0,0,50000000000, 0.02, false, false, false, 0.01, 50, -1
};

/*
const Status STAT0={
    .misaligment =  false,
    .radiation =  false,
    .fluctuation =  false,
    .fringe =  false,
    .edge =  false,
    .lossmap =  false,
    .fma =  false,
    .slice =  true,
    .combineddipole =  false,
    .computedrivingterms =  false,
    .leaderordertermonly =  false,
    .nonlineartermonly =  false,
    .linear =  true,
    .period =  true,
    .npara =  1,
    .totalslice =  0,
    .multipoleslice =  0,
    .Trx =  10,
    .Try =  10,
    .dp0 =  0,
    .dp =  0,
    .expand =  false,
    .particle = ELECTRON,
    .energy = 2E9,
    .second_order_chrom = false,
    .third_order_chrom = false,
    .nperiods = 1,
    .printout = false,
    .transfermatrix = true,
    .mincouple = 0.005,
    .track_lines = 13,
    .track_turns = 100
};
*/


struct MultipoleData{
    double betax, betay,etax,betaxp,betayp,etaxpp;
    double rbetax, rbetay;           /* rbetax = sqrt(betax) */
    double betax2, betay2;           /* betax2 = betax^2 */
    double phix, phiy;
    double b2L, b3L, s;
    complex <double> px[5], py[5];    /* px[j]=(exp(i*phix))^j j>0 */
};


class RDTsCache{
    public:
    size_t sext_slice;
    size_t mult_slice;
    vector<size_t> mult_position;
    vector<size_t> sext_position;
    vector<size_t> sext_slice_index;
    vector<MultipoleData> cache;
    RDTsCache(){
        sext_slice=0;
        mult_slice=0;
    }
    RDTsCache(const RDTsCache& cache0){
        sext_slice=cache0.sext_slice;
        mult_slice=cache0.mult_slice;
        mult_position=cache0.mult_position;
        sext_position=cache0.sext_position;
        sext_slice_index=cache0.sext_slice_index;
        cache=cache0.cache;
    }
};




#endif
