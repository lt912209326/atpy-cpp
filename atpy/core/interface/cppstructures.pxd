from .constants cimport *
cimport numpy as np

cdef extern from "../physics/utils/cppstructures.h"nogil:


    cdef struct Status:
        bint misaligment                   #misalignment and machine error considered
        bint radiation                     #synchrotron radiation dumping connsidered
        bint fluctuation                    #quantum excitation considered
        bint fringe                        #soft frindge considered
        bint edge                          #dipole edge will be calculated
        bint lossmap                       #loss data and information will be calculated
        bint fma                           #frequent map, tracking data will be produced and saved for fma
        bint slice                         #slice calculated, only exit of components calculated and less compute source needed if false
        bint combineddipole                #dipole slice counted when calculate drving terms if true
        bint computedrivingterms           #driving terms calculated
        bint leaderordertermonly           #only first-order driving terms calculated, less compute source needed
        bint nonlineartermonly             #only nonlinear terms calculated such as DTs, first-order chromaticity
        bint linear                        #only linear element calculated if true
        bint period                        #a periodic solution is needed
        size_t  npara                      #thread number for parallel computing
        size_t totalslice                  #total slice number of the whole beamline
        size_t multipoleslice              #total slice number of the whole nonlinear elements
        double Trx                         #the trace of horizontal transfer matrix
        double Try                         #the trace of vertical transfer matrix
        double dp0                         #actually dp/p0
        double dp                          #off-momentum range: (dp0-dp,dp0+dp)
        bint expand                        #whether enable different component parameters for the same element.
        size_t particle 
        double energy 
        bint second_order_chrom
        bint third_order_chrom
        size_t nperiods
        bint  printout
        bint transfermatrix
        double mincouple
        size_t track_lines
        size_t track_turns
        double monitor_dp                  #monitor to overlook the tunes at dp0 ± monitor_dp
        bint  lazy_compute
        double larger_monitor_dp
        bint fast_2nd_order_RDTs
        double max_betax
        double max_etax
        np.int64_t  NP
        double rf_dp
        bint rdt_fluctuation
        bint local_twiss
        bint off_momentum_rdts
        double off_rdts_observer
        np.uint32_t max_da_range 
        np.int64_t chrom_refpt

    cdef const Status STAT0

    cdef struct MultipoleData:
        double betax, betay,etax,betaxp,betayp,etaxpp
        double rbetax, rbetay           # rbetax = sqrt(betax) */
        double betax2, betay2           # betax2 = betax^2 */
        double phix, phiy
        double b2L, b3L, s
        complex[double] px[5]
        complex[double] py[5]    # px[j]=(exp(i*phix))^j j]0 */
    


    cdef cppclass RDTsCache:
        size_t sext_slice
        size_t mult_slice
        vector[size_t] mult_position
        vector[size_t] sext_position
        vector[size_t] sext_slice_index
        vector[MultipoleData] cache
        RDTsCache()except+
        RDTsCache(const RDTsCache& cache0)except+


