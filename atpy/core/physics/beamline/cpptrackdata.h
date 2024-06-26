#ifndef _CPPTRACKDATA_H_
#define _CPPTRACKDATA_H_

#include "cppconstants.h"

class TrackData{
    public:
    size_t nturn;
    size_t np;
    double *data;

    TrackData(){
        nturn=0;
        np=0;
        data=nullptr;
    }
    TrackData(size_t np0, size_t nturn0){
        nturn=nturn0;
        np=np0;
        data=(double*)calloc(4*np*nturn,sizeof(double) );
    }
    TrackData& operator= (TrackData& td){
        if(nturn*np==td.nturn*td.np){
            np=td.np;
            nturn=td.nturn;
        }
        else{
            if(data)free(data);
            np=td.np;
            nturn=td.nturn;
            data=(double*)calloc(4*np*nturn,sizeof(double) );
        }
        memcpy(data,td.data,4*nturn*np*sizeof(double));
        return *this;
    }

    int resize(const size_t np0, const size_t nturn0){
        if(nturn0*np0==nturn*np){
            np=np0;
            nturn=nturn0;
        }
        else{
            if(data)free(data);
            np=np0;
            nturn=nturn0;
            data=(double*)calloc(4*np*nturn,sizeof(double) );
        }
        return 0;
    }

    ~TrackData(){
        if(data){
            free(data);
            data=nullptr;
        }
    }
};

#endif