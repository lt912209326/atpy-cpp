#ifndef _CPPCOMPONENT_CPP_
#define _CPPCOMPONENT_CPP_

#include "cppcomponent.h"

CppComponent::CppComponent(){
    elem=nullptr;
    reverse=false;
    patch=false;
    e1=nullptr;
    e2=nullptr;
    fill(local,local+LOC_NUM,0);
    local[AX]=1;
    local[AY]=1;
    // tws=(double*)calloc(elem->nslice*TWS_NUM,__SIZEOF_DOUBLE__);
    // memset(local,0,elem->nslice*LOC_NUM*__SIZEOF_DOUBLE__);

}

CppComponent::CppComponent(CppElement* elem0, const bool reverse0, const size_t position0){
    // CppElementFactory factory;
    // elem=factory.CreateElement(elem0);
    elem=elem0;
    reverse=reverse0;
    patch=false;
    position=position0;
    // cout<<position<<endl;
    // cout<<this->position<<endl;
    e1=nullptr;
    e2=nullptr;
    // local=nullptr;
    // if(reverse && elem->kind==DIPOLE){
    //     if(!dynamic_cast<CppDipole*>(elem)->reverse_mat) elem->update_Matrix(reverse,0);
    // }
    fill(local,local+LOC_NUM,0);
    local[AX]=1;
    local[AY]=1;
    tws=matrix<double>(elem->nslice,TWS_NUM);
    // cache_tws=matrix<double>(elem->nslice,TWS_NUM);
    //  (double*)calloc(elem->nslice*LOC_NUM,__SIZEOF_DOUBLE__);
}

CppComponent::CppComponent(const CppComponent &comp0){
    if(elem){
        *elem=*(comp0.elem);
    }
    else{
        elem=comp0.elem;
    }
    reverse=comp0.reverse;
    patch=comp0.patch;
    position=comp0.position;
    e1=nullptr;
    e2=nullptr;
    // local=nullptr;
    memcpy(local,comp0.local,LOC_NUM*sizeof(double) );
    tws=comp0.tws;
    // cache_tws=comp0.cache_tws;
    // tws=matrix<double>(elem->nslice,TWS_NUM);
    // local=(double*)calloc(elem->nslice*LOC_NUM,__SIZEOF_DOUBLE__);
    // memset(local,0,elem->nslice*LOC_NUM*__SIZEOF_DOUBLE__);

}


CppComponent& CppComponent::operator= (const CppComponent& comp0){
    
    if(elem){
        *elem=*(comp0.elem);
    }
    else{
        elem=comp0.elem;
    }
    reverse=comp0.reverse;
    patch=comp0.patch;
    position=comp0.position;
    e1=nullptr;
    e2=nullptr;
    memcpy(local,comp0.local,LOC_NUM*sizeof(double) );
    tws=comp0.tws;
    // cache_tws=comp0.cache_tws;
    return *this;
}

CppComponent::~CppComponent(){
    reverse=false;
    patch=false;
    // memset(local,0,elem->nslice*LOC_NUM*__SIZEOF_DOUBLE__);
    elem=nullptr;
    fill(local,local+LOC_NUM,0);
    // free(local);
    // local=nullptr;
    if(e1)free(e1);
    if(e2)free(e2);
    e1=nullptr;
    e2=nullptr;

}


int CppComponent::display(const int disflag, const bool detail){
    cout<< std::setw(6)<< std::left  <<position<< std::setw(14)<< std::left  << elem->name;
    size_t ncount=0;
    
    vector<size_t> twiss={ COORD, BETAX, ALPHAX,GAMMAX, BETAY, ALPHAY,GAMMAY, ETAX, ETAPX, NUX, NUY,DCHROMX,DCHROMY};//,CHROMX,CHROMY };
    cout.precision(8);
    switch (disflag)
    {
    case KWD:
        /* code */
        // size_t keywords[]={ };
        for(int i=0;i<KWD_NUM;i++){
            ncount=count(elem->keywords.begin(),elem->keywords.end(), i);
            if(ncount>0){
                // cout << std::noshowpos<<std::scientific << std::setw(14) <<std::setprecision(6) << std::left << elem->values[i];
                cout << std::noshowpos<<std::defaultfloat <<std::setw(14) <<std::setprecision(6) << std::left << elem->values[i];
                // cout<<std::setiosflags(std::fixed)<<" "<<values[i];
            }
            else{
                // cout << std::noshowpos<<std::scientific << std::setw(14) <<std::setprecision(6) << std::left << "None";
                cout << std::noshowpos<<std::defaultfloat << std::setw(14) <<std::setprecision(6) << std::left << "None";
            }
        }
        break;
    case TWS:
        // size_t twiss[13]={ COORD, BETAX, ALPHAX, BETAY, ALPHAY, ETAX, ETAPX, NUX, NUY,DCHROMX,DCHROMY,CHROMX,CHROMY };
        if(detail && elem->nslice>1){
            cout<<":"<<endl;
            for(size_t j=0; j<elem->nslice;j++){
                cout<< std::setw(6)<< std::left  <<" "<< std::setw(14)<< std::left  << j+1;
                for(size_t i : twiss){
                    // cout << std::noshowpos<<std::scientific << std::setw(14) <<std::setprecision(6) << std::left << tws(j,i);
                    cout << std::noshowpos<<std::defaultfloat << std::setw(14) <<std::setprecision(6) << std::left << tws(j,i);
                }
                if(j<elem->nslice-1)cout<<endl;
            }
        }
        else{
            for(size_t i : twiss){
                // cout << std::noshowpos<<std::scientific << std::setw(14) <<std::setprecision(6) << std::left << tws(-1,i);
                cout << std::noshowpos<<std::defaultfloat <<std::setw(14) <<std::setprecision(6) << std::left << tws(-1,i);
            }
        }
        break;
    
    default:
        break;
    }
    cout<<endl;
    return 0;
}

int CppComponent::compute_TransferMatrix(const double* locin, const Status* stat){
    // if(stat->misaligment || fabs(stat->dp)>1e-8 ){
    //     memcpy(local+R11,locin+R11,36*sizeof(double));
    //     for(size_t i=0;i<6;i++){

    //     }
    // }
    elem->compute_TransferMatrix(locin, reverse, local  );
    return 0;
}


int CppComponent::update_TransferMatrix(double* rin, const Status* stat){
    elem->update_Matrix(reverse,rin+CODX ,stat);
    return 0;
}

int CppComponent::linearoptics(const double* twsin,const Status* stat, double* glbout){
    if(stat->slice && elem->nslice>1){
        for(int i=0;i<elem->nslice;i++){
            elem->linearoptics(twsin, ((float)(i+1))/elem->nslice,stat,reverse, &tws[i], glbout  );
            // cout<<"CppComponent::linearoptics: "<<elem->nslice<<endl;
        }
    }
    else{
        elem->linearoptics(twsin,1.0,stat,reverse, &tws[-1], glbout );
    }
    return 0;
}


int CppComponent::get_chromaticities(const double* twsin,const Status* stat){
    if( stat->slice && elem->nslice>1 ){
        for(int i=0;i<elem->nslice;i++){
            tws(i,CHROMX)=twsin[CHROMX] + tws(i,DCHROMX);
            tws(i,CHROMY)=twsin[CHROMY] + tws(i,DCHROMY);
        }
    }
    else{
        tws(-1,CHROMX)=twsin[CHROMX] + tws(-1,DCHROMX);
        tws(-1,CHROMY)=twsin[CHROMY] + tws(-1,DCHROMY);
    }
    return 0;
}



int CppComponent::get_geometry(const double* localin, const Status* stat){
    // 仅仅考虑在同一水平面的几何结构，需再完善
    double theta0=localin[THETAX];
    if(elem->kind !=DIPOLE){
        local[THETAX]=localin[THETAX];
    }
    else{
        theta0=localin[THETAX] + 0.5*elem->values[ANGLE];
        local[THETAX]=localin[THETAX]+elem->values[ANGLE];
    }
    local[S]=localin[S]+elem->values[L];
    local[GX]=localin[GX]+elem->values[L]*cos(theta0);
    local[GY]=localin[GY]+elem->values[L]*sin(theta0);
    return 0;
}

int CppComponent::track(double* rin, const Status* stat){
    if(rin[LOSS])return 0;
    if(patch) patch1(rin);
    // 入口判断是否丢失
    if(isnan(rin[X]) || isnan(rin[Y] )){
        rin[LOSS]=1;
        rin[LOSSPOS]=position;
    }
    if( (pow((rin[X]-local[DX])/local[AX],2)+pow((rin[Y]-local[DY])/local[AY],2))>1 ){
        rin[LOSS]=1;
        rin[LOSSPOS]=position;
    }
    elem->track(rin, stat,reverse);
    // 出口判断是否丢失
    if( (pow((rin[X]-local[DX])/local[AX],2)+pow((rin[Y]-local[DY])/local[AY],2))>1 ){
        rin[LOSS]=1;
        rin[LOSSPOS]=position;
    }
    if(patch) patch2(rin);
    return 0;

}


int CppComponent::patch1(double* rin){
    double theta1=-local[ROTATE1],  theta2=-local[ROTATE2],  theta3=-local[ROTATE3];
    double dx=local[DX],  dy=local[DY],  dz=local[DZ];
    double c1=cos(theta1), c2=cos(theta2), c3=cos(theta3);
    double s1=sin(theta1), s2=sin(theta2), s3=sin(theta3);
    double x,px,y,py,z,dp;
    x=rin[X]; px=rin[PX]; y=rin[Y]; py=rin[PY]; z=rin[Z]; dp=rin[DP];

    if(patch && !e1){
        double entrance[3][3]={{c1*c2,  c1*s2*s3-s1*c3,  c1*s2*c3+s1*s3},
                            {s1*c2,  s1*s2*s3+c1*c3,  s1*s2*c3-c1*s3},
                            {-s2  ,  c2*s3         ,  c2*c3         }
                            };

        double exits[3][3]={{c1*c2           ,  s1*c2           ,  -s2    },
                            {-c1*s2*s3-s1*c3 ,  -s1*s2*s3+c1*c3 ,  c2*s3 },
                            {c1*s2*c3+s1*s3  ,  s1*s2*c3-c1*s3  ,  -c2*c3 }
                            };

        e1=(double*)calloc(9,__SIZEOF_DOUBLE__ );
        e2=(double*)calloc(9,__SIZEOF_DOUBLE__ );
        memcpy(e1,entrance,9*__SIZEOF_DOUBLE__);
        memcpy(e2,entrance,9*__SIZEOF_DOUBLE__);
    }
    if(e1){
        rin[X] =e1[0*3+0]*x  + e1[0*3+1]*y  + e1[0*3+2]*z ;
        rin[PX]=e1[0*3+0]*px + e1[0*3+1]*py + e1[0*3+2]*dp;
        rin[Y] =e1[1*3+0]*x  + e1[1*3+1]*y  + e1[1*3+2]*z ;
        rin[PY]=e1[1*3+0]*px + e1[1*3+1]*py + e1[1*3+2]*dp;
        rin[Z] =e1[2*3+0]*x  + e1[2*3+1]*y  + e1[2*3+2]*z ;
        rin[DP]=e1[2*3+0]*px + e1[2*3+1]*py + e1[2*3+2]*dp;
    }
    return 0;
}


int CppComponent::patch2(double* rin){
    double x,px,y,py,z,dp;
    x=rin[X]; px=rin[PX]; y=rin[Y]; py=rin[PY]; z=rin[Z]; dp=rin[DP];
    if(e2){
        rin[X] =e2[0*3+0]*x  + e2[0*3+1]*y  + e2[0*3+2]*z ;
        rin[PX]=e2[0*3+0]*px + e2[0*3+1]*py + e2[0*3+2]*dp;
        rin[Y] =e2[1*3+0]*x  + e2[1*3+1]*y  + e2[1*3+2]*z ;
        rin[PY]=e2[1*3+0]*px + e2[1*3+1]*py + e2[1*3+2]*dp;
        rin[Z] =e2[2*3+0]*x  + e2[2*3+1]*y  + e2[2*3+2]*z ;
        rin[DP]=e2[2*3+0]*px + e2[2*3+1]*py + e2[2*3+2]*dp;
        rin[X] +=local[DX];
        rin[Y] +=local[DY];
    }
    return 0;
}



#endif