#ifndef _CPPBEAMLINE_CPP_
#define _CPPBEAMLINE_CPP_

#include "omp.h"
#include <Eigenvalues>
#include <Dense>
#include <exception>
#include "cppbeamline.h"
#include "twiss.h"
#include <vector>
#include <algorithm>

// using std::min;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::EigenSolver;
using Eigen::Map;
using Eigen::Dynamic;
using Eigen::RowMajor;
typedef Matrix<double,6,6,RowMajor> Matrix6d;

double touschekF(double x){
    // condition x > 0 is not checked in the function, be careful
    double Eular=0.5772;
    double logx, logx2, logx3, logx4, value;

    if(x<0.0013 ){
        return log(Eular/x)-1.5;
    }
    else if(x<10){
        logx=log(x);
        logx2=sqr(logx);
        logx3=logx*logx2;
        logx4=sqr(logx2);

        value=-3.10811-2.19156*logx - 0.615641*logx2 -0.160444*logx3 - 0.0460054*logx4 -0.0105172*logx*logx4 
                - 1.31192e-3*logx2*logx4 -6.3898e-5*logx3*logx4;
        return exp(value );
    }
    else{
        return 1e-16;
    }
}

void deallocate(CppBeamLine* line){
    if(line !=nullptr){
        delete line;
    }
}


CppBeamLine::CppBeamLine(){
    // elems_position["START"]={0};
    // elems.resize(0);
    // elems.emplace_back(new CppMarker("START") );
    
    nelems=0;
    length=0;
    CppElement* tmp_elem=new CppMarker("START");
    this->append(tmp_elem,false);
    delete tmp_elem;
    // elems["START"]->position.emplace_back(0);
    stat.Trx=10;
    stat.Try=10;
    stat.totalslice=0;
    stat.multipoleslice=0;
    fill(globals, globals+GLB_NUM,0.0);

    // rdtcache.mult_slice=0;
    // rdtcache.sext_slice=0;

}

CppBeamLine::CppBeamLine(const string name0, const size_t particle, const double energy,const size_t length0, Status stat0){
    name=name0;
    stat=stat0;
    stat.Trx=10;
    stat.Try=10;
    stat.totalslice=1;
    stat.multipoleslice=0;
    // cout<<"here in CppBeamLine::CppBeamLine(const string name0"<<endl;
    // rdtcache.mult_slice=0;
    // rdtcache.sext_slice=0;

    fill(globals, globals+GLB_NUM,0.0);
    globals[ENERGY]=energy;
    globals[MASS0]= particle==PROTON? 938272029.0:510998.9461;
    globals[GAMMA]=globals[ENERGY]/globals[MASS0] ;
    // elems["START"]=new CppMarker("START");
    nelems=0;
    length=0;
    line.resize(length0+1);
    // cout<<"CppBeamLine::CppBeamLine(: "<<length0<<endl;
    CppElement* tmp_elem=new CppMarker("START");
    this->append(tmp_elem,false);
    delete tmp_elem;
    // this->append(elems["START"],false);
    // elems["START"]->position.emplace_back(0);
}


CppBeamLine::CppBeamLine(const CppBeamLine& bline)
{
    name=bline.name;
    stat=bline.stat;
    factory=bline.factory;
    stat.Trx=10;
    stat.Try=10;
    // globals[ENERGY]=bline.globals[ENERGY];
    // globals[MASS0]= bline.globals[MASS0];
    // globals[GAMMA]= bline.globals[GAMMA];
    // memcpy(globals, bline.globals, GLB_NUM*__SIZEOF_DOUBLE__ );
    length=0;
    nelems=0;
    line.resize(bline.line.size() );
    for(auto iter:bline.line){
        this->append(iter->elem,iter->reverse);
        // cout<<"CppBeamLine::CppBeamLine(&):"<< (iter)->elem<<endl;
    }
    // rdtcache(bline.rdtcache);
    rdtcache=bline.rdtcache;
    // rdtcache.mult_slice=bline.rdtcache.mult_slice;
    // rdtcache.sext_slice=bline.rdtcache.sext_slice;

    // nelems=bline.nelems;
    // length=bline.length;
    memcpy(globals, bline.globals,GLB_NUM*sizeof(double));
}


CppBeamLine::~CppBeamLine(){
    double a=0;
    // cout<<"CppBeamLine not deallocated!"<<endl;
    for(auto &iter:elems){
        // cout<<iter.first<<endl;
        for(auto j :iter.second ){
            // cout<<"CppBeamLine::~CppBeamLine: "<<j<<endl;
            delete j;
        }
    }
    // cout<<"CppBeamLine::~CppBeamLine:0 "<<endl;
    for(auto iter:line){
        // cout<<"CppBeamLine::~CppBeamLine: "<<endl;
        delete iter;
    }
    if(stat.printout) cout<<"CppBeamLine::~CppBeamLine: "<<name<<" deallocated!"<<endl;

}


int CppBeamLine::save(CppElement* elem){
    if(elems.count(elem->name) ){
        *elem=**(elems[elem->name].begin() );
        return 0;
    }
    else{
        return -1;
    }
}



int CppBeamLine::display(const int disflag, const bool detail){
    
    vector<size_t> twiss={ COORD, BETAX, ALPHAX,GAMMAX, BETAY, ALPHAY,GAMMAX, ETAX, ETAPX, NUX, NUY,DCHROMX,DCHROMY};//,CHROMX,CHROMY };
    if(disflag!=DVTs){
        cout<< std::setw(6)<< std::left  <<"No."<< std::setw(14)<< std::left  << "Name";
    }
    cout.precision(8);
    switch (disflag)
    {
    case KWD:
        /* code */
        for(size_t i=0;i<KWD_NUM;i++){
            cout<< std::setw(14)<< std::left  << KEYWORDS_DICT.at(i);
        }
        break;
    case TWS:
        /* code */
        for(size_t i : twiss){
            cout<< std::setw(14)<< std::left  << TWISS_DICT.at(i);
        }
        break;
    case DVTs:
        for(size_t i=DQX;i<GLB_NUM;i++){
            cout<< std::setw(14)<< std::left  << GLOBALS_DICT.at(i);
            cout<< std::setw(14)<< std::right  << globals[i]<<endl;
        }
        break;

    default:
        break;
    }
    cout<<endl;
    if(disflag==DVTs) return 0;
    for(size_t i=0;i<=length;i++){
        line[i]->display(disflag,detail);
    }
    return 0;
}



int CppBeamLine::insert(CppElement* elem0,const int pos,const bool reverse){
    // 需要调整component中position参数
    if(fabs(pos)>length){ 
        throw std::out_of_range("insert position is out of range in function: CppBeamLine::insert!") ;
    }
    length+=1;
    if(0==elems.count(elem0->name)){
        // CppElement* tmp_elem=
        nelems+=1;
        elems[elem0->name ].emplace_back( factory.CreateElement(elem0) );
        elems_position[elem0->name]={length};
        elems[elem0->name ].back()->position={length};
    }
    else if(stat.expand){
        nelems+=1;
        elems[elem0->name].emplace_back( factory.CreateElement(elem0) );
        elems_position[elem0->name].emplace_back(length);
        elems[elem0->name ].back()->position={length};
    }
    else{
        elems_position[elem0->name].emplace_back(length);
        elems[elem0->name ].back()->position.emplace_back(length);
    }
    if(pos>0){
        line.insert(line.begin()+pos,new CppComponent(elems[elem0->name ].back(), reverse,pos) );
    }
    else{
        
        line.insert(line.end()+pos,new CppComponent(elems[elem0->name ].back(), reverse,pos) );
    }
    
    components[ elem0->name ].emplace_back(line.back() );

    stat.totalslice+=elem0->nslice;
    if(elem0->kind==QUADRUPOLE || elem0->kind==SEXTUPOLE || (stat.combineddipole && elem0->kind==DIPOLE)|| 
            elem0->kind==EXACTDRIFT || elem0->kind==OCTUPOLE ){
        stat.multipoleslice+=elem0->nslice;
        rdtcache.mult_slice+=elem0->nslice;
        if(elem0->kind==SEXTUPOLE)rdtcache.sext_slice+=elem0->nslice;
        rdtcache.mult_position.push_back(length);
    }

    return 0;

}


int CppBeamLine::append(CppElement* elem0,const bool reverse){
    // cout<<"here in CppBeamLine:append:"<<elem0->name<<endl;
    // cout<<"reverse: "<<reverse<<endl;
    // cout<<"CppBeamLine::append,elem0 address:"<< elem0<<endl;
    // cout<<"CppBeamLine::append,keywords pointer address:"<< &(elem0->keywords[0])<<endl;

    CppElement* tmp0;
    if(0==nelems){
        length=0;
    }
    else{
        length+=1;
    }
    if(length+1>line.size() ){
        cout<<"CppBeamLine::append: Input components are more than expected!"<<endl;
    }
    // length+=1;
    // cout<<"CppBeamLine::append: "<<0<<" : "<<elem0->name<<endl;
    // cout<<length<<endl;
    if(0==elems.count(elem0->name)){
        // cout<<"here in CppBeamLine:append:"<<0<<"length: "<<length<<endl;
        nelems+=1;
        // elems_position[elem0->name]->position.resize(0);
        tmp0=factory.CreateElement(elem0);
        elems[elem0->name].emplace_back( tmp0 );
        // cout<<"here in CppBeamLine:append:"<<1<<endl;
        // cout<<"CppBeamLine::append,elem0 address:"<< tmp0<<endl;
        line[length]= new CppComponent(elems[ elem0->name ].back(), (bool)reverse, length ) ;
        // cout<<"here in CppBeamLine:append:"<<2<<endl;
        elems_position[elem0->name]={length};
        elems[elem0->name ].back()->position={length};

    }
    else if(stat.expand){
        // cout<<1<<endl;
        
        nelems+=1;
        tmp0=factory.CreateElement(elem0);
        elems[elem0->name].emplace_back( tmp0 );
        // cout<<"CppBeamLine::append,elem0 address:"<< tmp0<<endl;
        line[length]= new CppComponent(elems[ elem0->name ].back(), reverse, length ) ;
        elems_position[elem0->name].emplace_back(length);
        elems[elem0->name ].back()->position={length};

        // elems.emplace_back(factory.CreateElement(elem0) );
        // line.emplace_back( elems.back(),reverse,length );
        
        // elems.back()->position={length};

    }
    else{
        
        // elems[elem0->name].emplace_back( factory.CreateElement(elem0) );
        line[length]=new CppComponent( elems[ elem0->name ].back(), reverse, length ) ;
        elems_position[elem0->name].emplace_back(length);
        elems[elem0->name ].back()->position.emplace_back(length);
    }
    // cout<<"CppBeamLine::append,elem0 address:"<< elems[elem0->name ].back()<<endl;
    components[ elem0->name ].emplace_back(line.back() );
    // line.emplace_back( elems[elem0->name],reverse,length );
    // elems[elem0->name]->position.emplace_back(length);

    // cout<<"CppBeamLine::append: "<<1<<endl;
    stat.totalslice+=elem0->nslice;
    if(elem0->kind==QUADRUPOLE || elem0->kind==SEXTUPOLE || (stat.combineddipole && elem0->kind==DIPOLE) ||
            elem0->kind==EXACTDRIFT || elem0->kind==OCTUPOLE) {
        stat.multipoleslice+=elem0->nslice;
        rdtcache.mult_slice+=elem0->nslice;
        if(elem0->kind==SEXTUPOLE)rdtcache.sext_slice+=elem0->nslice;
        rdtcache.mult_position.push_back(length);
    }
    
    // cout<<"here in CppBeamLine:append:"<<7<<endl;
    // cout<<"CppBeamLine::append: 2"<<line[0]->tws(0,1)<<endl;
    return 0;
}




int CppBeamLine::get_position_at_s(double coordinate){
    size_t pos;
    for(auto &iter:line ){
        pos=iter->position;
        if(iter->local[S] >=coordinate)break;
    }
    return pos;
}



int CppBeamLine::get_optics_at_s(double coordinate, double* twsout){
    size_t pos;
    int i_th=0;
    size_t i;
    double delta_len=0, elem_len=0, len_rate;
    double* twsin=nullptr;
    double tmp_glbs[GLB_NUM ];
    memcpy(tmp_glbs, globals, EMITX*__SIZEOF_DOUBLE__ );
    pos=get_position_at_s(coordinate);
    elem_len = line[pos]->elem->values[L];
    if(0==pos || fabs(elem_len)<1e-15 ){
        memcpy(twsout, &line[pos]->tws[0], TWS_NUM*__SIZEOF_DOUBLE__);
    }
    else{

        delta_len=coordinate-line[pos-1]->local[S];
        len_rate=delta_len/elem_len;
        twsin= &line[pos-1]->tws[-1] ;
        // 设置LHxxxxx离切片中最近的值
        for(i=LH00111;i<TWS_NUM;i++){
            twsout[i]= (1-len_rate)*twsin[i] + len_rate*line[pos]->tws(-1,i) ;
        }

        // i_th=round(delta_len/elem_len*line[pos]->elem->values[NSLICE ] )-1;
        // if(i_th>0){
        //     memcpy(twsout,&line[pos]->tws[i_th ], TWS_NUM*__SIZEOF_DOUBLE__);
        // }
        // else{
        //     memcpy(twsout,twsin, TWS_NUM*__SIZEOF_DOUBLE__);
        // }
        line[pos]->elem->linearoptics(twsin, delta_len/elem_len, &stat,line[pos]->reverse, twsout, tmp_glbs );
        twsin=nullptr;
    }
    return 0;
}


int CppBeamLine::computeGlobalProperties(){
    double gamma0=globals[GAMMA], energy0=globals[ENERGY];
    double *RI= &globals[0],time0, circumference, Jx ;
    double c_light=2.99792458e8 ;
    double C_gamma=8.85e-5;
    // cout<<"CppBeamLine::computeGlobalProperties: RI5"<<RI[RI5]<<" ,RI3"<<RI[RI3]<<endl;


    // #natural emittance
    globals[EMITX]=Cq*sqr(gamma0)*RI[RI5]/(RI[RI2] - RI[RI4]);
    globals[CIRCUMFERENCE]=line[length]->local[S];
    circumference=globals[CIRCUMFERENCE];
    time0=circumference/c_light;
    // #energy loss U0 [keV]
    globals[U0] = 1e6*C_gamma*0.5*INV_PI*pow(energy0*1e-9,4)*RI[RI2];
    // #momentum compaction factor
    globals[ALPHAC] = RI[RI1]/globals[CIRCUMFERENCE ];
    // #damping factor D Jx=1-D,Jy=1,Jz=2+D
    globals[DAMP_FACTOR] = RI[RI4]/RI[RI2];
    globals[JX] =1-globals[DAMP_FACTOR];
    globals[JY] =1;
    globals[JZ] =2+globals[DAMP_FACTOR];
    // #square of Energy dispersion
    globals[ESPREAD] = sqrt(Cq*gamma0*gamma0*RI[RI3]/(2*RI[RI2] + RI[RI4]) );
    if(fabs(RI[RI3])>1e-15){
        // #spin polarization I[6]:I3a in S.Y. Lee's book
        globals[SPIN] = -1.6/sqrt(3)*RI[RI3A]/RI[RI3];
    }
    Jx=globals[JX];
    //Damping time ms (eV->keV s->ms)
    globals[TAUY] = 2*energy0*time0/globals[U0] ;
    globals[TAUX] = globals[TAUY]/Jx ;
    globals[TAUZ] = globals[TAUY]/(3-Jx) ;
    globals[QX]=line[length]->tws(-1,NUX);
    globals[QY]=line[length]->tws(-1,NUY);
    // globals[H11001]=line[length]->tws(-1,CHROMX);
    // globals[H00111]=line[length]->tws(-1,CHROMY);
    RI=nullptr;
    return 0;
}


int CppBeamLine::initialoptics(){
    // 初始化起点光学函数
    memcpy(line[0]->local+R11, EYE66, sizeof(EYE66));
    // fill(&(line[0]->tws[0] ) , &(line[0]->tws[0] )+TWS_NUM, 1E10);

    // if(ordered_changed_position.size()>0){
    //     for(size_t index:ordered_changed_position){
    //         line[index]->update_TransferMatrix(line[index-1]->local ,&stat );
    //     }
    // }
    // else{
    //     for(int i=1;i<=length;i++){
    //         line[i]->update_TransferMatrix(line[i-1]->local ,&stat );
    //     }
    // }
    // if(stat.period || stat.transfermatrix){
    //     for(int i=1;i<length+1;i++){
    //         line[i]->compute_TransferMatrix(line[i-1]->local ,&stat );
    //     }
    // }
    line[0]->local[CODDP]=stat.dp0;
    findClosedOrbit(line[0]->local+CODX);
    if(stat.printout) cout<<"initialoptics 1"<<endl;
    if(stat.period){
        line[0]->tws=0.0;
        int nslice=line[length]->elem->nslice;
        double Trx= line[length]->local[R11]+line[length]->local[R22];
        double Try= line[length]->local[R33]+line[length]->local[R44];
        stat.Trx=Trx;
        stat.Try=Try;
        if(fabs( Trx)<2 && fabs(Try)<2  ){
            double* local=line[length]->local;
            double cosmux = 0.5*Trx;
            double cosmuy = 0.5*Try;
            double cscmux = local[R12]/( fabs(local[R12])*sqrt(1-cosmux*cosmux) ); //1/sin(mux);
            double cscmuy = local[R34]/( fabs(local[R34])*sqrt(1-cosmuy*cosmuy) ); //1/sin(muy);
            // 重设初始点LHxxxxx为0
            fill(&(line[0]->tws[0] ) , &(line[0]->tws[0] )+TWS_NUM, 0.0);

            line[0]->tws(0,BETAX ) = local[R12]*cscmux;
            line[0]->tws(0,ALPHAX) = 0.5*(local[R11]-local[R22])*cscmux;
            line[0]->tws(0,GAMMAX) = -local[R21]*cscmux;
            line[0]->tws(0,BETAY ) = local[R34]*cscmuy;
            line[0]->tws(0,ALPHAY) = 0.5*(local[R33]-local[R44])*cscmuy;
            line[0]->tws(0,GAMMAY) = -local[R43]*cscmuy;
            line[0]->tws(0,ETAX  ) = 0.5*(local[R16]*(1.0-local[R22]) + local[R12]*local[R26])/(1.0-cosmux);
            line[0]->tws(0,ETAPX ) = 0.5*(local[R26]*(1.0-local[R11]) + local[R21]*local[R16])/(1.0-cosmux);
            local=nullptr;
        }
        else{
            fill(&(line[0]->tws[0] ) , &(line[0]->tws[0] )+TWS_NUM, 1E10);
        }
    }
    else{
        line[0]->tws(0,GAMMAX) = (1+sqr(line[0]->tws(0,ALPHAX)  ))/line[0]->tws(0,BETAX)  ;
        line[0]->tws(0,GAMMAY) = (1+sqr(line[0]->tws(0,ALPHAY)  ))/line[0]->tws(0,BETAY)  ;
    }
    
    if(stat.printout) cout<<"initialoptics 2"<<endl;
    return 0;
}


int CppBeamLine::TwissPropagate(){
    int kind;
    if(stat.printout) cout<<"TwissPropagate 2"<<endl;
    // avoid integrate twice
    fill(globals+RI1,globals+RI5+1, 0.0);
    globals[NATURE_CHROMX]=0.0;
    globals[NATURE_CHROMY]=0.0;
    globals[TOTAL_K2L]=0.0;
    for(int i=1;i<=length;i++){
        kind=line[i]->elem->kind;
        line[i]->linearoptics(&(line[i-1]->tws[-1] ) ,&stat, globals);
        if(QUADRUPOLE==kind || DIPOLE==kind ){
            globals[NATURE_CHROMX]+=line[i]->tws(-1,DCHROMX);
            globals[NATURE_CHROMY]+=line[i]->tws(-1,DCHROMY);
        }
        else if(SEXTUPOLE==kind){
            globals[TOTAL_K2L]+=fabs(line[i]->elem->values[K2]*line[i]->elem->values[L]);
        }
        line[i]->get_geometry(line[i-1]->local ,&stat);
    }
    if(stat.printout) cout<<"TwissPropagate 2"<<endl;
    return 0;
}


int CppBeamLine::calculate(){
    double nux=0,nuy=0;
    size_t kind=0;
    // 考虑设置一些状态参数，当某些量改变时，状态值为true即需要重新分配空间，
    // 如insert函数、四六极铁切片数改变等
    //
    //
    //
    // 
    if(stat.computedrivingterms){
        if(rdtcache.cache.size()!=stat.multipoleslice){
            rdtcache.cache.resize(stat.multipoleslice); //分配四极铁六极铁切片内存空间，用于DrivingTerms计算
            cout<<"calculate: cache size      : "<<rdtcache.cache.size()<<endl;
        }
        if(rdtcache.sext_slice_index.size()!=rdtcache.sext_slice){
            rdtcache.sext_slice_index.resize(rdtcache.sext_slice );
            cout<<"calculate: sext cache size : "<<rdtcache.sext_slice<<endl;
        } 
    }
    if(!stat.nonlineartermonly){
        initialoptics();
        // 全局变量重设为0
        fill(globals+EMITX,globals+GLB_NUM, 0.0);
        // 将色品参数及全局驱动项设为1e12,方便优化
        fill(globals+H10010,globals+GLB_NUM, 1e12);
        globals[NATURE_CHROMX]=0;
        globals[NATURE_CHROMY]=0;
        globals[TOTAL_K2L]=0;
        double len;
        if( !stat.period || (stat.period && fabs(stat.Trx)<2.0 && fabs(stat.Trx)<2.0)  ){
            for(int i=1;i<=length;i++){
                kind=line[i]->elem->kind;
                line[i]->linearoptics(&(line[i-1]->tws[-1] ) ,&stat, globals);
                if(QUADRUPOLE==kind || DIPOLE==kind ){
                    globals[NATURE_CHROMX]+=line[i]->tws(-1,DCHROMX);
                    globals[NATURE_CHROMY]+=line[i]->tws(-1,DCHROMY);
                }
                else if(SEXTUPOLE==kind){
                    len=line[i]->elem->values[L];
                    globals[TOTAL_K2L]+=fabs(line[i]->elem->values[K2]*(len>0? len : 1) ) ;
                }
                line[i]->get_geometry(line[i-1]->local ,&stat);
            }
        }
    }
    else{
        // nonlinearonly=True,
        if(fabs(stat.Trx)>2 || fabs(stat.Trx)>2){
            if(stat.printout) cout<<"fabs(stat.Trx)>2 || fabs(stat.Trx)>2 || !stat.computedrivingterms"<<endl;
            // 通常初次计算光学函数
            fill(globals+EMITX,globals+GLB_NUM, 0.0);
            initialoptics();
            globals[NATURE_CHROMX]=0;
            globals[NATURE_CHROMY]=0;
            globals[TOTAL_K2L]=0;
            if( !stat.period || fabs(stat.Trx)>2.0 || fabs(stat.Trx)>2.0  ){
            // 在非周期解或者不稳定时，退出
                cout<<"period is "<<(stat.period?"True":"False")<<", Trx: "<<stat.Trx<<", Try: "<<stat.Try<<endl;
                throw std::logic_error("when nonlinearonly is true,period need be true,|Trx|<2 and |Try|<2 are needed!");
            }
            for(int i=1;i<=length;i++){
                kind=line[i]->elem->kind;
                line[i]->linearoptics(&(line[i-1]->tws[-1] ) ,&stat, globals);
                if(QUADRUPOLE==kind ||  DIPOLE==kind ){
                    globals[NATURE_CHROMX]+=line[i]->tws(-1,DCHROMX);
                    globals[NATURE_CHROMY]+=line[i]->tws(-1,DCHROMY);
                }
                else if(SEXTUPOLE==kind){
                    globals[TOTAL_K2L]+=fabs(line[i]->elem->values[K2]*line[i]->elem->values[L]);
                }
                line[i]->get_geometry(line[i-1]->local ,&stat);
            }
        }
        else{
            //已获得稳定周期解和光学函数，仅计算非线性元件一阶色品
            // for(int i=1;i<=length;i++){
            //     line[i]->linearoptics(&(line[i-1]->tws[-1] ) ,&stat, globals);
            //     line[i]->get_geometry(line[i-1]->local ,&stat);
            // }
            globals[TOTAL_K2L]=0;
            for(int i=1;i<=length;i++){
                if(line[i]->elem->kind==SEXTUPOLE){
                    line[i]->linearoptics(&(line[i-1]->tws[-1] ) ,&stat, globals);
                    globals[TOTAL_K2L]+=fabs(line[i]->elem->values[K2]*line[i]->elem->values[L]);
                }
            }
            for(size_t i=1;i<=length;i++){
                line[i]->get_chromaticities(&(line[i-1]->tws[-1] ) ,&stat);
            }
        }
    }
    computeGlobalProperties();
    compute_theory_touscheklifetime_part();
    // 色品计算，三阶驱动项计算
    if(stat.computedrivingterms && !stat.fast_2nd_order_RDTs && stat.period && fabs(stat.Trx)<2 && fabs(stat.Try)<2  ){
        if(stat.printout) cout<<"calculate: computmdrivingTerms_period: "<<endl;
        computmdrivingTerms_period(line, rdtcache.cache, rdtcache.mult_position, rdtcache.sext_slice_index, stat.multipoleslice,stat.nperiods, globals, stat );
    }
    else if( stat.computedrivingterms && stat.fast_2nd_order_RDTs && stat.period && fabs(stat.Trx)<2 && fabs(stat.Try)<2   ){
        if(stat.printout) cout<<"calculate: computmdrivingTerms_fast_period: "<<endl;
        computmdrivingTerms_fast_period(line, rdtcache.cache, rdtcache.mult_position, rdtcache.sext_slice_index, stat.multipoleslice,stat.nperiods, globals, stat );
    }
    else if( stat.computedrivingterms && !stat.period ){
        if(stat.printout) cout<<"calculate: computmdrivingTerms_nonperiod: "<<endl;
        computmdrivingTerms_nonperiod(line, rdtcache.cache, rdtcache.mult_position, rdtcache.sext_slice_index, stat.multipoleslice,stat.nperiods, globals, stat );
    }
    return 0;
}

int CppBeamLine::computeRDTs(){
    if(stat.computedrivingterms && !stat.fast_2nd_order_RDTs && stat.period && fabs(stat.Trx)<2 && fabs(stat.Try)<2  ){
        if(stat.printout) cout<<"calculate: computmdrivingTerms_period: "<<endl;
        computmdrivingTerms_period(line, rdtcache.cache, rdtcache.mult_position, rdtcache.sext_slice_index, stat.multipoleslice,stat.nperiods, globals, stat );
    }
    else if( stat.computedrivingterms && stat.fast_2nd_order_RDTs && stat.period && fabs(stat.Trx)<2 && fabs(stat.Try)<2   ){
        if(stat.printout) cout<<"calculate: computmdrivingTerms_fast_period: "<<endl;
        computmdrivingTerms_fast_period(line, rdtcache.cache, rdtcache.mult_position, rdtcache.sext_slice_index, stat.multipoleslice,stat.nperiods, globals, stat );
    }
    else if( stat.computedrivingterms && !stat.period ){
        if(stat.printout) cout<<"calculate: computmdrivingTerms_nonperiod: "<<endl;
        computmdrivingTerms_nonperiod(line, rdtcache.cache, rdtcache.mult_position, rdtcache.sext_slice_index, stat.multipoleslice,stat.nperiods, globals, stat );
    }
    return 0;
}

int CppBeamLine::compute_off_momentum_local_twiss(const double* tws0){

    double  precision=1e-10;
    Matrix<double,6,6> beams;
    // double nux=0, nuy=0, tmp_nux=0, tmp_nuy=0;
    double tws1[TWS_NUM]={0};
    // double tws2[TWS_NUM]={0};
    double *tws2;

    double *pbeam0=beams.data();
    double *tmp_local=&line[0]->local[0];

    beams.setIdentity(6,6);
    beams*=precision;
    for(int i=0;i<6;++i){
        for(int j=0;j<6;++j){
            pbeam0[6*i+j]+=line[0]->local[CODX+j];
        }
    }
    
    double tmpnux=0,tmpnuy=0,nux=0,nuy=0,localnux=0,localnuy=0;
    // if(stat.printout) cout<<beams<<endl;
    memcpy(tws1,tws0,__SIZEOF_DOUBLE__*TWS_NUM );
    // 跟踪
    for(size_t j=1;j<length+1;j++){
        for(size_t i=0;i<6;i++){
            line[j]->elem->track(pbeam0+6*i, &stat, line[j]->reverse);
            for(size_t k=0;k<6;k++){
                line[j]->local[R11+6*k+i]=(pbeam0[6*i+k] - line[j]->local[CODX+k] )/precision;
            }
        }
        tws2= &(line[j]->tws[-1]);
        propagate_twiss(tws2, tws1, &line[j]->local[R11]);

        tmpnux=tws2[NUX];
        tmpnuy=tws2[NUY];
        if(isnan(tmpnux) || isnan(tmpnuy) ){
            line[length]->local[LOCAL_NUX]=NAN;
            line[length]->local[LOCAL_NUY]=NAN;
            break;
        }
        if (tmpnux+1e-10 < nux){
            localnux+=(1+tmpnux-nux);
        }
        else{
            localnux+=(tmpnux-nux);
        }
        if (tmpnuy+1e-10 < nuy){
            localnuy+=(1+tmpnuy-nuy);
        }
        else{
            localnuy+=(tmpnuy-nuy);
        }
        if(TUNING == line[j]->elem->kind ){            
            if(line[j]->elem->values[DNUX]<0 || (line[j]->elem->values[DNUX]== 0.0 && tmpnux +1e-10 < nux) ) localnux-=1.0;
            if(line[j]->elem->values[DNUY]<0 || (line[j]->elem->values[DNUY]== 0.0 && tmpnuy +1e-10 < nuy) ) localnuy-=1.0;
        }
        else if( fabs(line[j]->elem->values[L])<1e-8 ){
            if(tmpnux +1e-10< nux ) localnux-=1.0;
            if(tmpnuy +1e-10< nuy ) localnuy-=1.0;
        
        }
        // if(stat.printout) cout<<i<<": "<<"nux: "<<localnux<<" , nuy: "<<localnuy<<endl;
        nux=tmpnux;
        nuy=tmpnuy;
        tws2[NUX]=localnux;
        tws2[NUY]=localnuy;
        line[j]->local[LOCAL_NUX]=localnux;
        line[j]->local[LOCAL_NUY]=localnuy;
    }
    if(stat.printout){
        cout<<"tmpnux: "<<tmpnux<<",tmpnuy: "<<tmpnuy<<endl;
        cout<<"nux: "<<tws2[NUX]<<",nuy: "<<tws2[NUY]<<endl;
    }
    return 0;
}

int CppBeamLine::findClosedOrbit(double* pop){
    double dp=pop[DP];
    double pop0[6]={0};
    Matrix<complex<double>,6,6,RowMajor> eigenvector;
    Matrix<double,Dynamic,Dynamic,RowMajor> solv;

    size_t MaxIter=40;
    double resume=1e10, convergence=1e-30, precision=fmax(1e-13, fmin(1e-12,fabs(dp)*1e-8) );
    double dt=0;
    Matrix<double,6,7> eye67;
    Matrix<double,6,1> beam0;
    Matrix<double,6,7> beams;
    
    double mux;
    double muy;

    Matrix<double,6,6>  lf, eye66;
    Matrix<double,4,4> m44,  eye44;
    Matrix<double,6,6,RowMajor> trsfmat;
    Matrix<double,6,1> xco, dco,d,tol;
    Matrix<double,4,1> b;
    // double *pbeam0=beam0.data();
    double *pbeam0=beams.data();
    beam0.fill(0);
    xco.fill(0);
    b.fill(0);
    dco.fill(0);
    eye67.setIdentity(6,7);
    eye44.setIdentity(4,4);
    // when dp=0,precision==0,then reset precision to a small number

    eye67*=precision;

    trsfmat.fill(0);
    lf.fill(0);
    tol.array()+=1e-30;
    lf(5,4)=1e-12;
    xco(5,0)=dp;
    for(size_t it=0;it<MaxIter;it++){
        beam0=xco;
        for(size_t i=0;i<7;i++){
            beams.col(i)=eye67.col(i)+xco;
        }
        // 跟踪
        // if(it<3) cout<<"before track"<<endl<<beams<<endl;
        for(size_t j=1;j<length+1;j++){
            for(size_t i=0;i<7;i++){
                line[j]->elem->track(pbeam0+6*i, &stat, line[j]->reverse);
            }
            memcpy(line[j]->local+CODX,pbeam0+6*6,6*__SIZEOF_DOUBLE__);
        }
        // if(it<3) cout<<"after track"<<endl<<beams<<endl;
        for(size_t i=0;i<6;i++){
            trsfmat.col(i)=(beams.col(i)-beams.col(6))/precision;
        }
        if( !stat.period){
            resume=0.5*convergence;
            break;
        }

        // if(fabs(trsfmat(5,4))<1e-10 ){
        //     for(size_t i=0;i<6;i++){
        //         trsfmat(4,i)=0;
        //         trsfmat(5,i)=0;
        //     }
        // }

        m44.block<4,4>(0,0)=trsfmat.block<4,4>(0,0);
        b.head(4)=(beams.col(6)-xco).head(4);

        // if( b.array().abs().sum()<convergence ){
        //     if(stat.printout) cout<<"Iter: "<<it<<", xco: "<<endl<<xco.transpose()<<endl;
        //     resume=0.1*convergence;
        //     break;
        // }

        if(isnan(b.dot(b) ) ){
            return -3;
        }
        dco.head(4)=(eye44-m44).fullPivHouseholderQr().solve(b).head(4) ;
        // dco.head(4)=(eye44-m44).bdcSvd().solve(b) ;

        // cout<<"Iter: "<<it<<", dco: "<<dco.transpose()<<endl<<"xco: "<<xco.transpose()<<endl;
        resume= dco.dot(dco);

        if(isnan(resume) ){
            return -2;
        }

        if(resume<convergence ){
            if(stat.printout){
                cout<<"Iter: "<<it<<",closed orbit found: "<<endl;
                cout<<"Delta X: "<<dco.transpose()<<endl;
                cout<<"particle input: "<<xco.transpose()<<endl;
                cout<<"particle out: "<<beams.col(6).transpose()<<endl;
            }
            break;
        }
        xco +=dco;
        for(int icoor=0;icoor<6;icoor++){
           if(fabs(xco(icoor,0) )<1e-15)  xco(icoor,0)= 0.0;
        }
    }
    // cout<<trsfmat<<endl;
    
    if(resume<convergence || !stat.period){
        memcpy(pop,xco.data(), 6*__SIZEOF_DOUBLE__ );
        memcpy(line[0]->local+CODX, xco.data(), 6*__SIZEOF_DOUBLE__ );
        for(size_t i=0;i<7;i++){
            beams.col(i)=eye67.col(i)+xco;
        }
        for(size_t j=1;j<length+1;j++){
            for(size_t i=0;i<7;i++){
                line[j]->elem->track(pbeam0+6*i, &stat, line[j]->reverse);
            }
            for(size_t i=0;i<6;i++){
                trsfmat.col(i)=(beams.col(i)-beams.col(6))/precision;
            }
            memcpy(line[j]->local+R11, trsfmat.data(), 36*__SIZEOF_DOUBLE__);
            memcpy(line[j]->local+CODX, pbeam0+36, 6*sizeof(double));
        }


        memcpy(line[length]->local+R11, trsfmat.data(), 36*__SIZEOF_DOUBLE__);
    }
    else{
        if(stat.printout)cout<<"findcod: Closed orbit not found!"<<endl<<"final xco: "<<endl<<xco.transpose()<<endl<<"dco: "<<endl<<dco.transpose()<<endl;
        return -1;
    }

    return 0;
}

int CppBeamLine::recover_twiss(){
    initialoptics();
    TwissPropagate();
    return 0;
}

int CppBeamLine::compute_sliced_off_momentum_twiss(){
    int iloop=0,islice=0;
    int kind;
    bool reverse;
    double *tws0,*tws1, *DynamicM66,*cod;
    double dp0=stat.dp0;
    stat.dp0=line[0]->local[CODDP];
    // if(stat.printout)cout<<"CppBeamLine::compute_sliced_off_momentum_twiss:1"<<endl;
    for(iloop=1;iloop<length+1;iloop++){
        kind=line[iloop]->elem->kind;
        if( SEXTUPOLE != kind && OCTUPOLE != kind && QUADRUPOLE != kind  ){
            continue;
        }
        else{
            reverse=line[iloop]->reverse;
            cod=&(line[iloop-1]->local[0])+CODX;
            line[iloop]->elem->update_Matrix(reverse, cod, &stat) ;
            DynamicM66 = line[iloop]->elem->DynamicM66;
            // if(stat.printout)cout<<"CppBeamLine::compute_sliced_off_momentum_twiss:2: "<<iloop<<endl;
            if(DynamicM66){
                tws0=&(line[iloop-1]->tws[-1]);
                for(islice=0;islice<line[iloop]->elem->nslice-1;islice++){
                    tws1=&(line[iloop]->tws[islice]);
                    // if(stat.printout){
                    //     cout<<"CppBeamLine::compute_sliced_off_momentum_twiss:2: "<<iloop<<": "<<islice<<endl;
                    //     cout<<tws0[BETAX]<<", "<<tws0[BETAY]<<", "<<tws0[NUX]<<", "<<tws0[NUY]<<", "<<tws0[ALPHAX]<<", "<<tws0[ALPHAY]<<endl;
                        
                    //     for(int i=0;i<TWS_NUM;i++){
                    //         cout<<tws0[i]<<" ";
                    //     }cout<<kind<<"M66:"<<endl;
                    //     for(int i=0;i<6;i++){
                    //         for(int j=0;j<6;j++){
                    //             cout<<DynamicM66[islice*36+i*6+j]<<" ";
                    //         }
                    //         cout<<endl;
                    //     }
                    //     cout<<endl;
                    // }
                    propagate_twiss( tws1, tws0,DynamicM66+islice*36 );
                    
                    // if(stat.printout)cout<<"CppBeamLine::compute_sliced_off_momentum_twiss:3: "<<iloop<<": "<<islice<<endl;
                    tws1[NUX] += tws0[NUX] ;
                    tws1[NUY] += tws0[NUY] ;
                }
            }
            else{
                line[iloop]->linearoptics(&(line[iloop-1]->tws[-1] ) ,&stat, globals);
            }
        }
    }
    stat.dp0=dp0;
    // if(stat.printout)cout<<"CppBeamLine::compute_sliced_off_momentum_twiss:4"<<endl;
    return 0;
}

int CppBeamLine::compute_off_momentum_RDTs(){
    double dpi=0.001, dp0=stat.dp0;
    double tws1[TWS_NUM]={0};
    int res1;
    string name;
    
    if(line[length]->tws(-1,NUX) <1e6 && line[length]->tws(-1,NUX) <1e6 ){
        register_RDTs(dp0,true);
    }
    else{
        register_RDTs(dp0,false);
    }
    for(int i=-1;i<2;i+=2){
        dpi=stat.dp0+i*stat.off_rdts_observer;
        res1=compute_off_momentum_twiss(tws1,dpi,true);
        if(res1<0){
            register_RDTs(dpi,false);
        }
        else{
            compute_sliced_off_momentum_twiss();
            computeRDTs();
            register_RDTs(dpi,true);
        }
    }
    recover_RDTs(stat.dp0);
    return 0;
}

int CppBeamLine::register_RDTs(double dpi, bool valid){
    int iglb;
    double  value;
    double tws1[TWS_NUM]={0};
    string num2str, tmp;
    string name;
    int effect_number,first_zero_pos;
    tmp=to_string( (int)floor(fabs(dpi)/1e-4 ) );
    first_zero_pos= tmp.find("0")==string::npos? tmp.length():tmp.find("0") ;
    // tmp="0100", effect_number=4-2=2 | tmp="0110", effect_number=4-1=3
    effect_number=4-(tmp.length() - first_zero_pos ) ;
    num2str=to_string(fabs(dpi) );
    if(dpi>0.0){
        num2str=string("@p") + num2str.substr(num2str.find(".")+1,effect_number);
    }
    else if(dpi<0.0){
        num2str=string("@m") + num2str.substr(num2str.find(".")+1,effect_number);
    }
    else{
        num2str=string("@0");
    }
    for(iglb=H10010;iglb<H00400+1;iglb++){
        name=GLOBALS_DICT.at(iglb)+num2str;
        value= valid? globals[iglb] : 1e10;
        if( std::find(id_table.id_table.begin(),id_table.id_table.end(), name) != id_table.id_table.end() ){
            ((Identity*)(id_table.id_dict[name]))->expr->value=value;
            // fresh Identity value, since AST-Refer type just get value of Identity, don't fresh Identity
            ((Identity*)(id_table.id_dict[name]))->calc();
        }
        else{
            id_table.id_table.push_back( name );
            value = value;
            id_table.id_dict[name]=new Identity(name,true, new Number(value) )   ;
        }
        if(stat.printout){
            cout<<name<<" : "<<id_table.id_dict[name]->value<<endl;
        }
    }
    return 0;
}


int CppBeamLine::recover_RDTs(double dp0){
    
    int iglb;
    double  value;
    double tws1[TWS_NUM]={0};
    string num2str, tmp;
    string name;
    int effect_number,first_zero_pos;
    tmp=to_string( (int)floor(fabs(0)/1e-4 ) );
    first_zero_pos= tmp.find("0")==string::npos? tmp.length():tmp.find("0") ;
    // tmp="0100", effect_number=4-2=2 | tmp="0110", effect_number=4-1=3
    effect_number=4-(tmp.length() - first_zero_pos ) ;
    num2str=to_string(fabs(dp0) );
    if(dp0>0.0){
        num2str=string("@p") + num2str.substr(num2str.find(".")+1,effect_number);
    }
    else if(dp0<0.0){
        num2str=string("@m") + num2str.substr(num2str.find(".")+1,effect_number);
    }
    else{
        num2str=string("@0");
    }
    for(iglb=H10010;iglb<H00400+1;iglb++){
        name=GLOBALS_DICT.at(iglb)+num2str;
        globals[iglb]= ((Identity*)(id_table.id_dict[name]))->expr->value;
        if(stat.printout){
            cout<<name<<" : "<<id_table.id_dict[name]->value<<endl;
        }
    }
    return 0;
}

int CppBeamLine::correctChrom(double dQx1, double dQy1){
    // chrom_corrector.aim_dQx
    double coeffx1=0, coeffx2=0, coeffy1=0, coeffy2=0; 
    size_t pos=0;
    double betax,betay,alphax,alphay,gammax,gammay,eta,etap;
    double len,len_2, len_3, len_4;
    double dqx=0,dqy=0;
    double dQx0, dQy0;
    double d0,d1,d2;
    double k1=0,k2=0;

    dQx0=line[length]->tws(-1,CHROMX) ;
    dQy0=line[length]->tws(-1,CHROMY);
    if(chrom_corrector.iscorr1){
        for(size_t it=0;it<chrom_corrector.position1.size();it++ ){
            pos=chrom_corrector.position1[it];
            len=line[pos ]->elem->values[L];
            len_2=len*len;
            len_3=len*len_2;
            len_4=len*len_3;
            // beta=0.5*(line[pos-1 ]->tws(-1,BETAX) + line[pos ]->tws(-1,BETAX) );
            // eta =0.5*(line[pos-1 ]->tws(-1,ETAX ) + line[pos ]->tws(-1,ETAX) );
            betax= line[pos-1 ]->tws(-1,BETAX) ;
            eta = line[pos-1 ]->tws(-1,ETAX ) ;
            etap = line[pos-1 ]->tws(-1,ETAPX ) ;
            betay=line[pos-1 ]->tws(-1,BETAY);
            // coeffx1+= chrom_corrector.coeff1[it]*eta*beta*len;

            // if(fabs(line[pos]->elem->values[K2])<1e-10 ){
                
            alphax= line[pos-1 ]->tws(-1,ALPHAX) ;
            alphay= line[pos-1 ]->tws(-1,ALPHAY) ;
            gammax= line[pos-1 ]->tws(-1,GAMMAX) ;
            gammay= line[pos-1 ]->tws(-1,GAMMAY) ;
            if(fabs(len)<1e-8){
                coeffx1 += chrom_corrector.coeff1[it]*eta*betax ;
                coeffy1 += chrom_corrector.coeff1[it]*eta*betay ;
            }else{
                coeffx1 += chrom_corrector.coeff1[it]*(0.25*etap*gammax*len_4 + (eta*gammax -2*etap*alphax )/3*len_3 + 
                            0.5*(etap*betax -2*eta*alphax )*len_2 + eta*betax*len );
                coeffy1 += chrom_corrector.coeff1[it]*(0.25*etap*gammay*len_4 + (eta*gammay -2*etap*alphay )/3*len_3 + 
                            0.5*(etap*betay -2*eta*alphay )*len_2 + eta*betay*len );
            }
            
        }
        dqx=4*PI*(dQx1-dQx0 );

    }
    if(chrom_corrector.iscorr2){
        for(size_t it=0;it<chrom_corrector.position2.size();it++ ){
            pos=chrom_corrector.position2[it];
            len=line[pos ]->elem->values[L];
            len_2=len*len;
            len_3=len*len_2;
            len_4=len*len_3;
            betax= line[pos-1 ]->tws(-1,BETAX) ;
            eta = line[pos-1 ]->tws(-1,ETAX ) ;
            etap = line[pos-1 ]->tws(-1,ETAPX ) ;
            betay=line[pos-1 ]->tws(-1,BETAY);
            alphax= line[pos-1 ]->tws(-1,ALPHAX) ;
            alphay= line[pos-1 ]->tws(-1,ALPHAY) ;
            gammax= line[pos-1 ]->tws(-1,GAMMAX) ;
            gammay= line[pos-1 ]->tws(-1,GAMMAY) ;
            if(fabs(len)<1e-8 ){
            coeffx2 += chrom_corrector.coeff2[it]*eta*betax;
            coeffy2 += chrom_corrector.coeff2[it]*eta*betay;
            }else{
                coeffx2 += chrom_corrector.coeff2[it]*(0.25*etap*gammax*len_4 + (eta*gammax -2*etap*alphax )/3*len_3 + 
                            0.5*(etap*betax -2*eta*alphax )*len_2 + eta*betax*len );
                coeffy2 += chrom_corrector.coeff2[it]*(0.25*etap*gammay*len_4 + (eta*gammay -2*etap*alphay )/3*len_3 + 
                            0.5*(etap*betay -2*eta*alphay )*len_2 + eta*betay*len );
            }
           
        }
        dqy=4*PI*(dQy1-dQy0 );

    }
    if(chrom_corrector.iscorr1 && chrom_corrector.iscorr2){
        d0= -coeffx1*coeffy2+coeffx2*coeffy1;
        d1= -(dqx*coeffy2+dqy*coeffx2);
        d2= (coeffx1*dqy+coeffy1*dqx);
        if (abs(d0)<1e-10 ){
            k1=1e6;
            k2=1e6;
        }
        else{
            k1=d1/d0;
            k2=d2/d0;
        }
        for(size_t it=0;it<chrom_corrector.unique_position1.size();it++ ){
            pos=chrom_corrector.unique_position1[it];
            line[pos ]->elem->values[K2]+=k1*chrom_corrector.unique_coeff1[it] ;
        }
        for(size_t it=0;it<chrom_corrector.unique_position2.size();it++ ){
            pos=chrom_corrector.unique_position2[it];
            line[pos ]->elem->values[K2]+=k2*chrom_corrector.unique_coeff2[it] ;
        }
    }
    else if(chrom_corrector.iscorr1 && !chrom_corrector.iscorr2){
        if(abs(coeffx1)<1e-30 ) coeffx1=1e-10;
        k1=dqx/coeffx1;
        for(size_t it=0;it<chrom_corrector.unique_position1.size();it++ ){
            pos=chrom_corrector.unique_position1[it];
            line[pos ]->elem->values[K2]+=k1*chrom_corrector.unique_coeff1[it] ;
        }
    }
    else if(!chrom_corrector.iscorr1 && chrom_corrector.iscorr2){
        if(abs(coeffy2)<1e-30 ) coeffy2=1e-10;
        k2= -dqy/coeffy2;
        for(size_t it=0;it<chrom_corrector.unique_position2.size();it++ ){
            pos=chrom_corrector.unique_position2[it];
            line[pos ]->elem->values[K2]+=k2*chrom_corrector.unique_coeff2[it] ;
        }
    }
    return 0;

}

int CppBeamLine::compute_off_momentum_twiss(double* tws, double dp, bool local_twiss=true){
    // Status stat0=stat;
    // double dp0=stat.dp0;
    // double Trx,Try;
    int result, result2;
    // double OneTurnTransferMatrixCache[6][6];
    // double tws0[TWS_NUM]={0.0};
    double tws1[TWS_NUM]={0.0},tws2[TWS_NUM]={0.0};
    double cosmux=0, sinmux=0, cosmuy=0, sinmuy=0;
    double *tmpM= &line[length]->local[0];


    fill(tws,tws+TWS_NUM,0.0);

    // memcpy(OneTurnTransferMatrixCache, &line[length]->local[R11], 36*__SIZEOF_DOUBLE__ );
    memcpy(tws1, &line[0]->tws[0], TWS_NUM*__SIZEOF_DOUBLE__ );
    tws1[TWSMODE]=1;
    // memcpy(tws0, &line[length]->tws[0], TWS_NUM*__SIZEOF_DOUBLE__ );


    // if(stat.period){
    //     tws0[NUX]=(line[length]->local[R12]>0? acos(0.5*stat.Trx):PIx2-acos(0.5*stat.Trx) )/PIx2;
    //     tws0[NUY]=(line[length]->local[R34]>0? acos(0.5*stat.Try):PIx2-acos(0.5*stat.Try) )/PIx2;
    // }
    // else{
    //     tws0[NUX]= fmod( line[length]->tws(-1,NUX), 1.0);
    //     tws0[NUY]= fmod( line[length]->tws(-1,NUY), 1.0);
    // }
    // must recovery at the end
    fill(line[0]->local+CODX,line[0]->local+CODDP+1,0);

    line[0]->local[CODDP]=dp;
    result=findClosedOrbit(line[0]->local+CODX );
    
    if(!stat.period ){
        // non-period
        propagate_twiss(&tws2[0], tws1, &tmpM[R11]);
        memcpy(tws,tws2,__SIZEOF_DOUBLE__*TWS_NUM );

        if(!(isnan(tws[NUX] ) || isnan(tws[NUY] )) ){
            if(local_twiss){
                compute_off_momentum_local_twiss(tws1);
                tws[NUX]=line[length]->local[LOCAL_NUX];
                tws[NUY]=line[length]->local[LOCAL_NUY];
            }
        }
        else{
            tws[NUX]=1e10;
            tws[NUY]=1e10;
            return -1;
        }
    }
    else{
        // period
        if(result<0){
            tws[NUX]=1e10;
            tws[NUY]=1e10;
            return -1;
        }
        result2 = calculate_coupled_period_twiss( &line[length]->local[R11],tws2 );
        memcpy(tws,tws2,__SIZEOF_DOUBLE__*TWS_NUM );
        if(result2<0){
            if(result2 == -1){
                tws[NUX]=1e10;
                tws[NUY]=1e10;
                if(stat.printout)cout<<"Failed to decouple periodic transfer matrix. The linear matrix is unstable. "<<endl;
            }
            else if(result2 == -2){
                tws[NUX]=1e8;
                tws[NUY]=1e8;
                if(stat.printout)cout<<"Failed to decouple the periodic transfer matrix with mode 1"<<endl;
            }
            else if(result2 == -3){
                // tws[NUX]=
                if(stat.printout){
                    cout<<"Failed to get beta functions. The linear matrix is unstable."<<endl;
                    cout<<"Trx"<<line[length]->local[R11]+line[length]->local[R22]<<endl;
                    cout<<"Try"<<line[length]->local[R33]+line[length]->local[R44]<<endl;
                }
            }
            return -2;
        }
        if(local_twiss){
            compute_off_momentum_local_twiss(tws);
            tws[NUX]=line[length]->local[LOCAL_NUX];
            tws[NUY]=line[length]->local[LOCAL_NUY];
        }        
    }
    if(stat.printout){
        cout<<"TWS2[NUX]: "<<tws2[NUX]<<", TWS2[NUY]: "<<tws2[NUY]<<endl;
        cout<<"TWS[NUX]: "<<tws[NUX]<<", TWS[NUY]: "<<tws[NUY]<<endl;
        cout<<"TWS[BETAX]: "<<tws[BETAX]<<", TWS[BETAY]: "<<tws[BETAY]<<endl;
        cout<<"TWS[ALPHAX]: "<<tws[ALPHAX]<<", TWS[ALPHAY]: "<<tws[ALPHAY]<<endl;
        cout<<"TWS[GAMMAX]: "<<tws[GAMMAX]<<", TWS[GAMMAY]: "<<tws[GAMMAY]<<endl;
        cout<<"Delta: "<<line[0]->local[CODDP]<<", nux: "<<tws[NUX]<<", nuy: "<<tws[NUY]<<endl;
    }
    if(isnan(tws[NUX] ) ){ tws[NUX]=1e10; line[length]->local[LOCAL_NUX]=tws[NUX] ; } 
    if(isnan(tws[NUY] ) ){ tws[NUY]=1e10; line[length]->local[LOCAL_NUY]=tws[NUY] ; }

    // stat=stat0;
    // fill(line[0]->local+CODX,line[0]->local+CODDP+1,0);
    // line[0]->local[CODDP]=stat.dp0;
    // for(int j=1;j<=length;j++){
    //     fill(line[j]->local+CODX,line[j]->local+CODDP+1,0);
    // }
    // memcpy( &line[length]->local[R11],OneTurnTransferMatrixCache, 36*__SIZEOF_DOUBLE__ );
    return 0;

}


int CppBeamLine::compute_large_off_momentum_tunes(double dp){
    double tws1[TWS_NUM]={0};
    double int_nux0=floor(line[length]->tws(0,NUX) );
    double int_nuy0=floor(line[length]->tws(0,NUY) );

    if(stat.printout) cout<<"int_nux: "<<int_nux0<<",int_nuy: "<<int_nuy0<<endl;
    compute_off_momentum_twiss(tws1, stat.dp0-dp);
    globals[LOW_QX]=tws1[NUX]-int_nux0;
    globals[LOW_QY]=tws1[NUY]-int_nuy0;
    if(stat.printout) cout<<"low_nux: "<<globals[LOW_QX]<<",low_nuy: "<<globals[LOW_QY]<<endl;
    compute_off_momentum_twiss(tws1, stat.dp0+dp);
    globals[HIGH_QX]=tws1[NUX]-int_nux0;
    globals[HIGH_QY]=tws1[NUY]-int_nuy0;
    if(stat.printout) cout<<"high_nux: "<<globals[HIGH_QX]<<",high_nuy: "<<globals[HIGH_QY]<<endl;
    return 0;
}


int CppBeamLine::compute_off_momentum_sum_square(double dp){
    double tws1[TWS_NUM]={0};
    double nux0,nuy0, int_nux0, int_nuy0;
    double nux_bounds[2], nuy_bounds[2];
    double sum_sqr_qx=0,sum_sqr_qy=0, dpi=0, ddp;
    double dpmax1=0.0,dpmax2=1.0;
    double dQx1=0,dQx2=0, dQy1=0,dQy2=0;
    double prev_nux=0,prev_nuy=0;
    double predict_nux=0,predict_nuy=0;
    int dprange[]={4,3,3,2,2,2,1,1,1,1,0 };
    int cnt=0, iloop=0;
    int ngrid=10, ndiv=20;
    int ret1=0;
    if(stat.local_twiss ){
        nux0= line[length]->tws(0,NUX) ;
        nuy0= line[length]->tws(0,NUY) ;
    }
    else{
        nux0= fmod( line[length]->tws(0,NUX),1.0);
        nuy0= fmod( line[length]->tws(0,NUY),1.0);
    }
    int_nux0 = nux0 ;
    int_nuy0 = nuy0 ;
    prev_nux = int_nux0 ;
    prev_nuy = int_nuy0 ;
    predict_nux = int_nux0 ;
    predict_nuy = int_nuy0 ;

    nux_bounds[0] = 0.5*floor(nux0/0.5 ) ;
    nux_bounds[1] = nux_bounds[0]+0.5 ;
    nuy_bounds[0] = 0.5*floor(nuy0/0.5 ) ;
    nuy_bounds[1] = nuy_bounds[0]+0.5 ;
    
    ddp=abs(dp)/(ndiv);
    dpi=0;
    cnt=0;
    dpmax1=0.0;
    dpmax2=1.0;
    if(stat.printout){
        cout<<"compute_off_momentum_sum_square::nux0: "<<int_nux0<<", nuy0: "<<int_nuy0<<endl;
    }
    for(iloop=0;iloop<ngrid; iloop++){
        dpi+=dprange[iloop]*ddp;
        // if(stat.printout) cout<<"int_nux: "<<int_nux0<<",int_nuy: "<<int_nuy0<<endl;
        ret1=compute_off_momentum_twiss(tws1, dpi, stat.local_twiss );
        
        if( (tws1[NUX]<nux_bounds[0] || tws1[NUX]>nux_bounds[1]) || (tws1[NUY]<nuy_bounds[0] || tws1[NUY]>nuy_bounds[1] )  ){
            break;
        }
        if( (predict_nux<nux_bounds[0] || predict_nux>nux_bounds[1]) || (predict_nuy<nuy_bounds[0] || predict_nuy>nuy_bounds[1] )  ){
            dQx2=dQx1;
            dQy2=dQy1;
            break;
        }
        dQx1 = (tws1[NUX]-prev_nux)/(dprange[iloop]*ddp);
        dQy1 = (tws1[NUY]-prev_nuy)/(dprange[iloop]*ddp);

        sum_sqr_qx+=dprange[iloop]*sqr(tws1[NUX]-int_nux0) ;
        sum_sqr_qy+=dprange[iloop]*sqr(tws1[NUY]-int_nuy0) ;

        predict_nux = tws1[NUX]+dQx1*ddp*dprange[iloop+1] ;
        predict_nuy = tws1[NUY]+dQy1*ddp*dprange[iloop+1] ;
        prev_nux = tws1[NUX];
        prev_nuy = tws1[NUY];
        
        // if(tws1[NUX]>1e8 && tws1[NUY]>1e8 ) break;
        // if(tws1[NUX]<nux_bounds[0] || tws1[NUX]>nux_bounds[1]  ) {
        //     sum_sqr_qx+= tws1[NUX]>9.99e5 ? tws1[NUX] : 1e4*dprange[iloop]*sqr(tws1[NUX]-int_nux0) ;
        // }
        // else{
        //     sum_sqr_qx+=dprange[iloop]*sqr(tws1[NUX]-int_nux0) ;
        // }
        // if(tws1[NUY]<nuy_bounds[0] || tws1[NUY]>nuy_bounds[1] ) {
        //     sum_sqr_qy+= tws1[NUY]>9.99e5 ? tws1[NUY] : 1e4*dprange[iloop]*sqr(tws1[NUY]-int_nuy0) ;
        // }
        // else{
        //     sum_sqr_qy+=dprange[iloop]*sqr(tws1[NUY]-int_nuy0) ;
        // }
        cnt+=dprange[iloop];
        if(stat.printout){
            cout<<"compute_off_momentum_sum_square::nux: "<<tws1[NUX]<<", nuy: "<<tws1[NUY]<<endl;
        }
    }
    // sum_sqr_qx+=(ndiv-cnt)*1e10 ;
    // sum_sqr_qy+=(ndiv-cnt)*1e10 ;
    dpmax1+=cnt*ddp;
    dpmax2*=cnt*ddp;
    
    prev_nux = int_nux0 ;
    prev_nuy = int_nuy0 ;
    predict_nux = int_nux0 ;
    predict_nuy = int_nuy0 ;
    dQx1=dQx2=dQy1=dQy2=0.0;
    
    dpi=0;
    cnt=0;
    for(iloop=0;iloop<ngrid; iloop++){
        dpi-=dprange[iloop]*ddp;
        // if(stat.printout) cout<<"int_nux: "<<int_nux0<<",int_nuy: "<<int_nuy0<<endl;
        ret1=compute_off_momentum_twiss(tws1, dpi, stat.local_twiss );
        
        
        if( (tws1[NUX]<nux_bounds[0] || tws1[NUX]>nux_bounds[1]) || (tws1[NUY]<nuy_bounds[0] || tws1[NUY]>nuy_bounds[1] )  ){
            break;
        }
        if( (predict_nux<nux_bounds[0] || predict_nux>nux_bounds[1]) || (predict_nuy<nuy_bounds[0] || predict_nuy>nuy_bounds[1] )  ){
            dQx2=dQx1;
            dQy2=dQy1;
            break;
        }
        dQx1 = (tws1[NUX]-prev_nux)/(dprange[iloop]*ddp);
        dQy1 = (tws1[NUY]-prev_nuy)/(dprange[iloop]*ddp);

        sum_sqr_qx+=dprange[iloop]*sqr(tws1[NUX]-int_nux0) ;
        sum_sqr_qy+=dprange[iloop]*sqr(tws1[NUY]-int_nuy0) ;

        predict_nux = tws1[NUX]+dQx1*ddp*dprange[iloop+1] ;
        predict_nuy = tws1[NUY]+dQy1*ddp*dprange[iloop+1] ;
        prev_nux = tws1[NUX];
        prev_nuy = tws1[NUY];
        // if(tws1[NUX]>1e8 && tws1[NUY]>1e8 ) break;
        // if(tws1[NUX]<nux_bounds[0] || tws1[NUX]>nux_bounds[1]  ) {
        //     sum_sqr_qx+= tws1[NUX]>9.99e5 ? tws1[NUX] : 1e4*dprange[iloop]*sqr(tws1[NUX]-int_nux0) ;
        // }
        // else{
        //     sum_sqr_qx+=dprange[iloop]*sqr(tws1[NUX]-int_nux0) ;
        // }
        // if(tws1[NUY]<nuy_bounds[0] || tws1[NUY]>nuy_bounds[1] ) {
        //     sum_sqr_qy+= tws1[NUY]>9.99e5 ? tws1[NUY] : 1e4*dprange[iloop]*sqr(tws1[NUY]-int_nuy0) ;
        // }
        // else{
        //     sum_sqr_qy+=dprange[iloop]*sqr(tws1[NUY]-int_nuy0) ;
        // }
        cnt+=dprange[iloop];
        if(stat.printout){
            cout<<"compute_off_momentum_sum_square::nux: "<<tws1[NUX]<<", nuy: "<<tws1[NUY]<<endl;
        }
    }
    
    dpmax1+=cnt*ddp;
    dpmax2*=cnt*ddp;
    // sum_sqr_qx+=(ndiv-cnt)*1e10 ;
    // sum_sqr_qy+=(ndiv-cnt)*1e10 ;

    globals[SUM_SQR_QX]=sum_sqr_qx/(2*ndiv)+dpmax1*1e10 ;
    globals[SUM_SQR_QY]=sum_sqr_qy/(2*ndiv)+dpmax2*1e10 ;
    return 0;
}



int CppBeamLine::compute_theory_touscheklifetime_part(){
    int cnt=0, iloop, max_betax_index, max_etax_index;
    double h0=0, e_spread, emitx, couple,gamma, max_betax, max_etax, betaxi, betayi, etaxi ,ax1, ax2;
    double sigma_xi, sigma_yi, delta_p_accept, sigma_xpi, zetai, inv_tau=0;

    e_spread=globals[ESPREAD ];
    emitx=globals[EMITX ];
    gamma=globals[GAMMA ];
    couple=stat.mincouple;
    max_betax=0;
    max_etax=0;
    if(stat.max_betax<1e-8 && stat.max_etax<1e-8 ){
        for(iloop=0;iloop<length+1;++iloop){
            betaxi = line[iloop]->tws(-1,BETAX);
            etaxi = line[iloop]->tws(-1,ETAX);
            if(max_betax <  betaxi){
                max_betax = betaxi;
                max_betax_index=iloop;
            }
            if(max_etax < etaxi ){
                max_etax = etaxi;
                max_etax_index=iloop;
            }
        }
    }
    else{
        max_betax=stat.max_betax;
        max_etax=stat.max_etax;
        max_betax_index=-1;
        max_etax_index=-1;
    }
    
    for(iloop=0;iloop<length+1;++iloop){
        if(MARKER == line[iloop]->elem->kind  ) continue;
        h0 = 0.5*(line[iloop]->tws(0,H0)+line[iloop]->tws(-1,H0) ) ;

        betaxi = 0.5*(line[iloop]->tws(0,BETAX)+ line[iloop]->tws(-1,BETAX) ) ;
        betayi = 0.5*(line[iloop]->tws(0,BETAY)+ line[iloop]->tws(-1,BETAY) ) ;
        etaxi = 0.5*(line[iloop]->tws(0,ETAX)+ line[iloop]->tws(-1,ETAX) ) ;

        sigma_xi=sqrt(betaxi*emitx +sqr(e_spread*etaxi) );
        sigma_yi=sqrt(betayi*emitx*couple );
        sigma_xpi= emitx/sigma_xi*sqrt(1+ h0*sqr(e_spread)/emitx );
        if(max_betax_index<0 || max_etax_index<0 ){
            
            ax1=line[0 ]->local[AX ];
            delta_p_accept = ax1/(sqrt(h0*max_betax) + max_etax  );
        }
        else{    
            ax1=line[max_betax_index ]->local[AX ];
            ax2=line[max_etax_index ]->local[AX ];
            delta_p_accept=fmin( ax1/(sqrt(h0*max_betax)+ line[max_betax_index]->tws(-1,ETAX) )
                                ,ax2/(sqrt(h0*line[max_etax_index]->tws(-1,BETAX))+ max_etax )  ) ;
        }
        delta_p_accept=fmin(stat.rf_dp, delta_p_accept);
        
        zetai = sqr(delta_p_accept)/sqr(gamma*sigma_xpi );
        
        // if(stat.printout ){
        //     cout<<"delta_p: "<<delta_p_accept<<endl;
        //     cout<<"zeta: "<<zetai<<endl;
        // }
        inv_tau += touschekF(zetai )/(sigma_xi*sigma_yi*sigma_xpi*sqr(delta_p_accept) )*line[iloop]->elem->values.at(L) ;
    }

    globals[INV_TAU ]= inv_tau*sqr(2.81e-15)*2.99792458e8*stat.NP/(8*PI*gamma*gamma*gamma*globals[CIRCUMFERENCE ] ) ;
    return 0;
}


int CppBeamLine::computeSecondOrderChromaticities(const double dp){
    // Status stat0=stat;
    double dp0=stat.dp0;
    double Trx,Try;
    int result;
    double OneTurnTransferMatrixCache[6][6];
    double tws1[TWS_NUM],tws2[TWS_NUM];
    double cosmux=0, sinmux=0, cosmuy=0, sinmuy=0;
    Matrix<double,5,1> nux, nuy,etax, etapx, betax, betay ,alphax, alphay ,coeff_nux,coeff_nuy ;
    Matrix<double,5,5> fit_mat_x, fit_mat_y ;
    double i=0;
    int64_t ref_pt=stat.chrom_refpt ; 
    if(stat.chrom_refpt<0){
        ref_pt=length+1+stat.chrom_refpt ; 
    }
    else{
        ref_pt=stat.chrom_refpt ; 
    }
    double *tmpM= &line[ref_pt]->local[0];
    i=-3.0;
    nux.fill(0);
    nuy.fill(0);
    fit_mat_x.fill(0);
    fit_mat_y.fill(0);
    memcpy(OneTurnTransferMatrixCache, &line[ref_pt]->local[R11], 36*__SIZEOF_DOUBLE__ );
    memcpy(tws1, &line[0]->tws[0], TWS_NUM*__SIZEOF_DOUBLE__ );
    
    // if(stat.period){
    //     nux[2]=(line[ref_pt]->local[R12]>0? acos(0.5*stat.Trx):PIx2-acos(0.5*stat.Trx) )/PIx2;
    //     nuy[2]=(line[ref_pt]->local[R34]>0? acos(0.5*stat.Try):PIx2-acos(0.5*stat.Try) )/PIx2;
    // }
    // else{
    //     nux[2]= fmod( line[ref_pt]->tws(-1,NUX), 1.0);
    //     nuy[2]= fmod( line[ref_pt]->tws(-1,NUY), 1.0);
    // }
    
    nux[2]= line[ref_pt]->tws(-1,NUX);
    nuy[2]= line[ref_pt]->tws(-1,NUY);
    if(stat.printout) cout<<nux[2]<<","<<nuy[2]<<endl;
    
    betax[2]= line[ref_pt]->tws(-1,BETAX); 
    alphax[2]= line[ref_pt]->tws(-1,ALPHAX); 
    betay[2]= line[ref_pt]->tws(-1,BETAY); 
    alphay[2]= line[ref_pt]->tws(-1,ALPHAY); 
    etax[2] = line[ref_pt]->tws(-1,ETAX); 
    etapx[2]= line[ref_pt]->tws(-1,ETAPX); 

    // stat.computedrivingterms=false;
    if(stat.printout) cout<<"dp: "<<dp<<endl;
    for(int j=0;j<5;j++){
        i+=1.0;
        
        fit_mat_x(j,0)=1;
        fit_mat_y(j,0)=1;
        if(2==j){
            continue;
        }
        else if(!stat.third_order_chrom && (0==j || 4==j ) ){
            continue;
        }
        fill(line[0]->local+CODX,line[0]->local+CODDP+1,0);

        line[0]->local[CODDP]=dp0+i*dp;
        
        for(size_t k=1;k<5;++k){
            fit_mat_x(j,k)=pow(line[0]->local[CODDP],k);
            fit_mat_y(j,k)=pow(line[0]->local[CODDP],k);
        }
        result=compute_off_momentum_twiss(tws2, line[0]->local[CODDP], true);
        
        betax[j]= line[ref_pt]->tws(-1,BETAX); 
        alphax[j]= line[ref_pt]->tws(-1,ALPHAX); 
        betay[j]= line[ref_pt]->tws(-1,BETAY); 
        alphay[j]= line[ref_pt]->tws(-1,ALPHAY); 
        etax[j] = line[ref_pt]->tws(-1,ETAX); 
        etapx[j]= line[ref_pt]->tws(-1,ETAPX); 
        nux[j]= line[ref_pt]->tws(-1,NUX); 
        nuy[j]= line[ref_pt]->tws(-1,NUY); 

        // betax[j] = tws2[BETAX];
        // alphax[j]= tws2[ALPHAX];
        // betay[j] = tws2[BETAY];
        // alphay[j]= tws2[ALPHAY];
        // nux[j] = tws2[NUX];
        // nuy[j] = tws2[NUY];
        // etax[j] = tws2[ETAX];
        // etapx[j] = tws2[ETAPX];

        // if(nux[j]-nux[2]>0.5  ){ // nux0~0.0x, nux[j]~0.9x
        //     nux[j]-=1;
        // }
        // else if(nux[2]-nux[j]>0.5){ //nux0~0.9x, nux[j]~0.0x
        //     nux[j]+=1;
        // }
        // if(nuy[j]-nuy[2]>0.5  ){ // nuy~0.0x, nuyj]~0.9x
        //     nuy[j]-=1;
        // }
        // else if(nuy[2]-nuy[j]>0.5){ //nuy~0.9x, nuyj]~0.0x
        //     nuy[j]+=1;
        // }

        if(stat.printout){
            cout<<"Delta: "<<line[0]->local[CODDP]<<", nux: "<<nux[j]<<", nuy: "<<nuy[j]<<endl;
            cout<<" betax: "<<betax[j]<<", betay: "<<betay[j]<<endl;
        }
    }
    // stat=stat0;
    fill(line[0]->local+CODX,line[0]->local+CODDP+1,0);
    line[0]->local[CODDP]=stat.dp0;
    findClosedOrbit(line[0]->local+CODX );
    // for(int j=1;j<=length;j++){
    //     fill(line[j]->local+CODX,line[j]->local+CODDP+1,0);
    //     // line[j]->update_TransferMatrix(line[j-1]->local, &stat);
    // }
    // calculate();
    
    memcpy( &line[ref_pt]->local[R11],OneTurnTransferMatrixCache, 36*__SIZEOF_DOUBLE__ );
    globals[DQX]=0.5*(nux[3]-nux[1])/dp;
    globals[DQY]=0.5*(nuy[3]-nuy[1])/dp;
    // globals[DQX] = 0.33333333*(nux[3]+nux[2]-2.0*nux[1])/dp;
    // globals[DQY] = 0.33333333*(nuy[3]+nuy[2]-2.0*nuy[1])/dp;
    globals[D2QX]=(nux[1]+nux[3]-2*nux[2])/(2*dp*dp);
    globals[D2QY]=(nuy[1]+nuy[3]-2*nuy[2])/(2*dp*dp);
    globals[DETAX]=(etax[3] - etax[1])/dp ;
    globals[DETAPX]=(etapx[3] - etapx[1])/dp ;
    globals[DBETAX]=(betax[3] - betax[1])/dp ;
    globals[DBETAY]=(betay[3] - betay[1])/dp ;
    globals[DALPHAX]=(alphax[3] - alphax[1])/dp ;
    globals[DALPHAY]=(alphay[3] - alphay[1])/dp ;
    globals[DDBETAX]=(betax[1]+betax[3]-2*betax[2])/(2*dp*dp);
    globals[DDBETAY]=(betay[1]+betay[3]-2*betay[2])/(2*dp*dp);
    globals[DDETAX]=(etax[1]+etax[3]-2*etax[2])/(2*dp*dp);
    
    double A1X,B1X, A1Y,B1Y;
    A1X= globals[DALPHAX]- alphax[2]*globals[DBETAX]/betax[2] ;
    B1X= globals[DBETAX]/betax[2] ;
    A1Y= globals[DALPHAY]- alphay[2]*globals[DBETAY]/betay[2] ;
    B1Y= globals[DBETAY]/betay[2] ;

    globals[WX]=sqrt( sqr(A1X) + sqr(B1X)   );
    globals[WY]=sqrt( sqr(A1Y) + sqr(B1Y)   );

    if(stat.third_order_chrom){
        Matrix<double,5,5> tmp_fit_mat;
        tmp_fit_mat=fit_mat_x.transpose();
        coeff_nux= (tmp_fit_mat*fit_mat_x).householderQr().solve(tmp_fit_mat*nux);
        tmp_fit_mat=fit_mat_y.transpose();
        coeff_nuy= (tmp_fit_mat*fit_mat_y).householderQr().solve(tmp_fit_mat*nuy);
        if(stat.printout){
            cout<<fit_mat_x<<endl;
            cout<<fit_mat_y<<endl;

            cout<<coeff_nux<<endl;
            cout<<coeff_nuy<<endl;
        }
        globals[D3QX]=coeff_nux[3];
        globals[D3QY]=coeff_nuy[3];
        globals[D4QX]=coeff_nux[4];
        globals[D4QY]=coeff_nuy[4];
        // globals[D3QX]=(nux[4]-nux[0]+2*nux[1]-2*nux[3])/(12*dp*dp*dp);
        // globals[D3QY]=(nuy[4]-nuy[0]+2*nuy[1]-2*nuy[3])/(12*dp*dp*dp);
    }
    return 0;

}



int CppBeamLine::track(double* beams, size_t nbegin, size_t nend, size_t nturn_begin, size_t nturn_end ){
    size_t nturn0=nturn_begin, nturn1=nturn_end;
    size_t stat_no=0, first_end_pos=length;
    uint32_t i=0,j=0,k=0;
    if(nturn_begin > nturn_end){
        throw std::invalid_argument("nturn_begin should no more than nturn_end!");
    }
    if (0==nbegin && nend == length){/*绝大多数情况粒子从头到尾*/ }
    else{
        if(nbegin<nend+1 && nturn_end-nturn_begin >0){
            if(nend>length) throw std::invalid_argument("nend should no more than length of beamline!");
            stat_no=1;
            nturn0=nturn_begin+1;
            nturn1=nturn_end-1;
        }
        else if(nbegin<nend+1 ){ /*&& nturn_end-nturn_begin <1*/
            // 后续swith中不需要继续跟踪
            if(nend>length) throw std::invalid_argument("nend should no more than length of beamline!");
            stat_no=2;
            first_end_pos=nend;
            nturn0=nturn_begin+1;
        }
        else if(nbegin >nend){
            if(nbegin>length) throw std::invalid_argument("nbegin should no more than length of beamline!");
            stat_no=3;
            nturn0=nturn_begin+1;
        }
        for(j=nbegin;j<first_end_pos+1;j++){
            line[j]->track(beams,&stat);
            if(beams[LOSS]) break;
        }
        if(beams[LOSS]){
            beams[LOSSTURN]=nturn_begin;
            return 0;
        }
    }
    // 多圈，从头到尾的跟踪
    for( k=nturn0;k<nturn1+1;k++){
        for( j=1;j<length+1;j++){
            line[j]->track(beams, &stat);
            if(beams[LOSS]){
                if(j<nbegin ){
                    beams[LOSSTURN]=k-1;
                }
                else{
                    beams[LOSSTURN]=k;
                }
                return 0;
            }
        }
    }
    
    switch(stat_no){
        case 0:
        case 2:
            break;
        case 1:
        case 3:
            for(size_t j=1;j<nend+1;j++){
                line[j]->track(beams,&stat);
                if(beams[LOSS]){
                    beams[LOSSTURN]=nturn_end;
                    return 0;
                }
            }
            break;
        default:
            break;
    }
    beams[LOSSTURN]=nturn_end+1;
    return 0;
}

double CppBeamLine::get_DA_area(double* area_datas){
    
    size_t npara=stat.npara, max_machine_thread=1;
    double s=0;
    int alive_cnt[100]={0}; // maximum track lines is 100
    size_t nline=stat.track_lines, n_step=nline-1, track_turns=stat.track_turns; //
    matrix<double> beams(nline,TRK_NUM,0); //
    vector<double> max_radiuses(nline,0);  //
    size_t j=0,k=0;
    int i=0;

    if(stat.npara<1) npara=1;
    max_machine_thread=omp_get_num_procs() ;
    if(stat.npara>2*max_machine_thread ){
        std::cout<<"Warning: Parallel number is too big, npara="<<stat.npara<<std::endl;
        npara=2*max_machine_thread-1;
    }
    if(track_turns<2) throw std::invalid_argument("track_turns should be no less than 2!");

    double fact0=0;
    double betax1=line[0]->tws(-1,BETAX), betay1=line[0]->tws(-1,BETAY);
    double cos_theta=cos(i*PI/n_step), sin_theta=sin(i*PI/n_step) ;
    double mincoup=stat.mincouple ;
    double p_k=0,q_k=0,n_k=0, a2_yk=1, a2_xk=1, xo_k=0, xo_begin=0;
    double betax_k=0, betay_k=0;
    double min_sqrt_A=1e10;
    double max_radius=1;
    double inv_couple=1/(1+mincoup);
    double commax1=sqrt(betax1*1*inv_couple),commay1=sqrt(betay1*mincoup*inv_couple) ;
    double sigmax = commax1*sqrt(globals[EMITX]), sigmay = commay1*sqrt(globals[EMITX]);
    double ytol=1e-8;
    uint32_t max_da_range = fmin(500, fmax( 10, stat.max_da_range) ) ;
    for(j=1;j<length+1;j++){
        // ref OPA user manual
        a2_xk=pow(line[j]->local[AX],2);
        a2_yk=pow(line[j]->local[AY],2);
        betax_k=line[j]->tws(-1,BETAX);
        betay_k=line[j]->tws(-1,BETAY);
        xo_k=line[j]->local[CODX];
        n_k=a2_yk*betax_k*(1-mincoup)+a2_xk*betay_k*mincoup;
        p_k=a2_yk*fabs(xo_k)*sqrt((1-mincoup)*betax_k )/n_k;
        q_k=a2_yk*(pow(xo_k,2)- a2_xk)/n_k;
        min_sqrt_A=fmin(min_sqrt_A, -p_k+sqrt(pow(p_k,2)-q_k) );
    }
    xo_begin=line[0]->local[CODX];
    max_radius=min_sqrt_A;
    double dr = max_da_range/50.0;
    // r_step=0.1*max_radius;
    omp_set_num_threads(npara);
#pragma omp parallel for
    for(i=0;i<nline;i++){
        double r_step=0.1*max_radius, r_precision=r_step;
        double x,y,r,delta_x;
        double cos_theta=cos(i*PI/n_step), sin_theta=sin(i*PI/n_step);
        // double lose_cnt=0;
        // r=max_radius-r_step;
        r = max_da_range;
        max_radiuses[i]=r;
        while(r>=0.0 && alive_cnt[i] <5 ){ //
            fill(&beams[i], &beams[i]+TRK_NUM,0);
            beams(i,DP)=stat.dp0;
            x=r*cos_theta*sigmax+xo_begin;
            y=r*sin_theta*sigmay;
            beams(i,X)=x;
            beams(i,Y)=y+ ytol;
            track(&beams[i],0,length,1,track_turns );
            if( !beams(i,LOSS) ){
                r-=dr;
                alive_cnt[i]+=1;
                if(stat.printout){ cout<<i<<": "<<r<<","<<max_radiuses[i]<<endl; }
            }
            else{
                r-=dr;
                alive_cnt[i]=0;
                max_radiuses[i]=r;
            }
        }
    }
    double normal_area=0, sigma_r=0, ave_r=0, sum_r2_i=0;
    double r=0, r_prev=0, r_next=0, r_local_max=0;
    for(i=0;i<nline;i++){
        r=max_radiuses[i];
        ave_r+=r;
        sum_r2_i+=r*r;
        if(i>0){
            normal_area+=0.5*sin(PI/n_step)*max_radiuses[i]*max_radiuses[i-1];
        }
        area_datas[i]= r*cos(i*PI/n_step)*sigmax+ xo_begin ;
        area_datas[nline+i]= r*sin(i*PI/n_step)*sigmay+ytol ;
    }
    ave_r/=nline;
    sigma_r=sqrt(sum_r2_i/nline - ave_r*ave_r);
    // normal_area-=0.5*n_step*PI/nline*sigma_r*sigma_r;
    globals[DA]=normal_area;
    globals[DA_SIGMA]=sigma_r;
    return normal_area;
}

// int CppBeamLine::paralleltrack(int np, int nturn, double* beams){
//     int npara=1, max_machine_thread=1;
//     double s=0;
//     if(stat.npara<1) npara=1;
//     max_machine_thread=omp_get_num_procs() ;
//     if(stat.npara>2*max_machine_thread ){
//         std::cout<<"Warning: Parallel number is too big, npara="<<stat.npara<<std::endl;
//         npara=2*max_machine_thread-1;
//     }
//     omp_set_num_threads(npara);
//     if(stat.fma) trackdata.resize(np,nturn);
//     if(stat.fma){
// #pragma omp parallel for
//         for(int i=0;i<np;i++){
//             for(int k=0;k<nturn;k++){
//                 for(int j=0;j<=length;j++){
//                     line[j]->track(beams+i*TRK_NUM,&stat);
//                     if(beams[i*TRK_NUM+LOSS])break;
//                 }
//                 if(beams[i*TRK_NUM+LOSS]){
//                     beams[i+TRK_NUM+LOSSTURN]=k+1;
//                     break;
//                 }
//                 memcpy(trackdata.data+(4*nturn*i+4*k) ,beams+(i*TRK_NUM), 4*sizeof(double) );
//             }
//         }
//     }
//     else{
// #pragma omp parallel for
//         for(int i=0;i<np;i++){
//             for(int k=0;k<nturn;k++){
//                 for(int j=0;j<=length;j++){
//                     line[j]->track(beams+i*TRK_NUM,&stat);
//                     if(beams[i*TRK_NUM+LOSS])break;
//                 }
//                 if(beams[i*TRK_NUM+LOSS]){
//                     beams[i+TRK_NUM+LOSSTURN]=k+1;
//                     break;
//                 }
//             }
//         }
//     }
//     return 0;
// }


#endif