#ifndef _CPPELEMENT_H_
#define _CPPELEMENT_H_

#include "cppconstants.h"
#include "cppstructures.h"

class CppElement
{
    public:
        int kind;
        size_t nslice;
        string name;
        vector<size_t> position;
        vector<int> keywords;
        unordered_map<int,double> values;
        double M66[6][6];
        double T66[6][6];
        double *DynamicM66;
    CppElement();
    CppElement(const string);
    CppElement(const unordered_map<int,double> arg);
    CppElement(const CppElement&);
    CppElement& operator=(const CppElement& elem);
    inline int set_keywords(int kwd,double value0){
        values[kwd]=value0;
        return 0;
    }
    inline int display(){
        size_t ncount=0;
        cout<< std::setw(13)<< std::left  <<name;
        for(int i=0;i<KWD_NUM;i++){
            ncount=count(keywords.begin(),keywords.end(), i);
            cout << std::noshowpos << std::setw(13) <<std::setprecision(6) << std::left << (ncount>0?values[i]:0);
        }
        cout<<endl;
        return 0;
    }
        // virtual inline int get_keywords_ptr(int kwd,double value);
        virtual int linearoptics(const double* twsin, const double len_rate, const Status* stat,const bool reverse, double* twsout, double* glbout)=0;
        virtual int track(double* rin, const Status* stat, const bool reverse);
        virtual int update_Matrix(const bool reverse, const double* cod , const Status* stat)=0;
        virtual int compute_TransferMatrix(const double* R66in,const bool reverse,double* R66out);
        // virtual inline int update()=0;

        virtual ~CppElement(){
            position.clear();
            position.shrink_to_fit();
            keywords.clear();
            keywords.shrink_to_fit();
            if(DynamicM66){
                free(DynamicM66);
                DynamicM66=nullptr;
            }
            // cout<<"~CppElement:"<<name<<endl;
        }


};

int matmul66(double* C, double* A, double* B);


#endif