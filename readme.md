# atpy 重构计划
class Component{

public:
 int position;
 int eid;
 int nslices;
 int sliceid;
}

class BeamLine{

    xt::xarray<double,N,TWS_NUM> twiss;
    xt::xarray<double,N,6,7> TM;
    xt::xarray<complex,N,RDTs_NUM> RDTs_cache;
    
    xt::xarray<double,N,TWS_NUM> cache_twiss;
    xt::xarray<double,N,6,7> cache_TM;
    std::unordered_map<string, AST*> globals;

}