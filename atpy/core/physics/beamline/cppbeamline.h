#ifndef _CPPBEAMLINE_H_
#define _CPPBEAMLINE_H_
#include "cppelementfactory.h"
#include "cppcomponent.h"
#include "cppstructures.h"

#include "cpptrackdata.h"
#include "cppdrivingterms.h"
#include "cppoptimization.h"


double touscheF(double x);

class CppBeamLine{
    public:
    CppElementFactory factory;
    unordered_map<string, vector<size_t>> elems_position;
    // unordered_map<string, AST*> name_table;
    // unordered_map<string, Assign*> var_table;
    // vector<string> ordered_name_list
    // vector<string> ordered_var_list
    bool need_reinit;
    IdTable  id_table;
    Variables vars ;
    Constraints constraints ;
    Optima optima;
    ChromCorrector chrom_corrector;

    set<size_t> ordered_changed_position;
    unordered_map<string,vector<CppElement*> > elems;
    unordered_map<string,vector<CppComponent*> > components;
    vector<CppComponent*> line;
    RDTsCache rdtcache;
    string  name;
    size_t length,nelems;                              //number of components in the beam line
    Status stat;                                //the calculation status of a beam line; like radiation,quatumn excitation,driving term,slice,polarization,independent component,on-momentum
    double globals[GLB_NUM];          //global parameters like natural emittance,natural chromaticity...
    // TrackData  trackdata;             //trackdata to record the position of all particles/turn
    CppBeamLine();
    CppBeamLine(const string name0, const size_t particle, const double energy,const size_t length0, Status stat0);
    CppBeamLine(const CppBeamLine&);
    // CppBeamLine(CppBeamLine&);         //global parameters like natural emittance,natural chromaticity...
    ~CppBeamLine();
    int save(CppElement* );
    int display(const int,const bool);
    int displayTransferMatrix(){
        for(int k=0;k<=11;k++){
            cout<<"Transfer Matrix at"<<k<<"-th element:"<<endl;
            for(int i=0;i<6;i++){
                for(int j=0;j<6;j++){
                cout<< std::setw(13)<< std::left  << line[k]->local[i*6+j];
                }
                cout<<endl;
            }
        }
        return 0;
    }
    int insert(CppElement*,const int pos,const bool reverse=false);
    int append(CppElement*,const bool reverse);
    int get_position_at_s(double coordinate);
    int get_optics_at_s(double coordinate, double* twsout);
    int TwissPropagate();
    int calculate();                            //calculate beam line according to above status
    int correctChrom(double dQx1, double dQy1);
    int compute_off_momentum_twiss(double* tws, double dp, bool local_twiss);
    int compute_large_off_momentum_tunes(double dp);
    int compute_off_momentum_sum_square(double dp);
    int compute_theory_touscheklifetime_part();
    int computeSecondOrderChromaticities(const double);
    int initialoptics();
    int linearoptics();
    int computeGlobalProperties();
    int computeRDTs();
    int compute_off_momentum_local_twiss(const double* tws0);
    int findClosedOrbit(double* pop);
    int recover_twiss();
    int compute_sliced_off_momentum_twiss(); // sextupole and octupole only
    int compute_off_momentum_RDTs();
    int register_RDTs(double dpi,bool valid);
    int recover_RDTs(double dp0);
    int track(double* beams, size_t nbegin, size_t nend, size_t nturn_begin, size_t nturn_end );
    double get_DA_area(double* area_datas);
    int paralleltrack(int np, int nturn, double* beams);

};


void deallocate(CppBeamLine* line);

#endif