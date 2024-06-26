#ifndef _CPPOPTIMIZATION_H_
#define _CPPOPTIMIZATION_H_

#include "cppconstants.h"
#include "cppast.h"
#include <utility>

using std::make_pair;

class IdTable{
    public:
    unordered_map<string,AST*> id_dict;
    vector<string>   id_table;
    size_t num_id;
    IdTable(){
        num_id=0;
    }
    ~IdTable(){
        for(auto it: id_table){
            delete id_dict.at(it);
        }
    }
};


struct Variables{
    unordered_map<string,AST*>        independent_vars;
    unordered_map<string,AST*>        dependent_vars;
    vector<string>  ordered_independ_var_names;
    vector<string>  ordered_depend_var_names;
    set<string>     leftvalue_names;
    size_t  num_independent_vars;
    size_t  num_dependent_vars;
    Variables(){
        num_independent_vars=0;
        num_dependent_vars=0;
    }
    void append(string name0, AST* var,bool independent){
        if(independent){
            independent_vars.insert(make_pair(name0,var) );
            ordered_independ_var_names.push_back(name0);
            num_independent_vars+=1;
        }
        else{
            dependent_vars.insert(make_pair(name0,var) );
            ordered_depend_var_names.push_back(name0);
            num_dependent_vars+=1;
        }
    }
    ~Variables(){
        // cout<<"~Variables:before"<<endl;
        if (num_independent_vars>0){
            // cout<<"Variables::~Variables"<<endl;
            for(auto iter:independent_vars){
                delete iter.second;
            }
        }
        if (num_dependent_vars>0){
            // cout<<"Variables::~Variables"<<endl;
            for(auto iter:dependent_vars){
                delete iter.second;
            }
        }
        // cout<<"~Variables:after"<<endl;
    }
};





struct Constraints{
    vector<AST*>   values;
    vector<string> ordered_constraint_names;
    unordered_map<size_t,bool> time_consuming_terms;
    size_t num_constraint;
    Constraints(){
        num_constraint=0;
    }
    
    void append(string name0, AST* constraint){
        ordered_constraint_names.push_back(name0 );
        values.push_back(constraint );
        num_constraint+=1;
    }

    ~Constraints(){
        // cout<<"~Constraints:before"<<endl;
        if(num_constraint>0){
            // cout<<"Constraints::~Constraints"<<endl;
            for(auto iter:values){
                delete iter;
            }
        }
        // cout<<"~Constraints:after"<<endl;
    }
};



struct Optima{
    vector<AST*> values;
    vector<string> ordered_optima_names;
    unordered_map<size_t,bool> time_consuming_terms;
    vector<double> minormax;
    size_t num_optima;
    Optima(){
        num_optima=0;
    }
    
    void append(string name0, AST* optima,bool ismin){
        ordered_optima_names.push_back(name0 );
        values.push_back(optima );
        if(ismin){
            minormax.push_back(1.0);
        }
        else{
            minormax.push_back(-1.0);
        }
        num_optima+=1;
    }
    ~Optima(){
        // cout<<"~Optima:before"<<endl;
        if(num_optima>0){
            // cout<<"Optima::~Optima"<<endl;
            for(auto iter:values){
                delete iter;
            }
        }
        // cout<<"~Optima:after"<<endl;
    }
};


struct ChromCorrector{
    bool iscorr1;
    bool iscorr2;
    string corrector1;
    double aim_dQx;
    vector<string> dependent_corr1;
    vector<size_t> position1;
    vector<size_t> unique_position1;
    vector<double> coeff1;
    vector<double> unique_coeff1;
    // vector<size_t> elem_position1;
    // vector<double> elem_coeff1;
    string corrector2;
    double aim_dQy;
    vector<string> dependent_corr2;
    vector<size_t> position2;
    vector<size_t> unique_position2;
    vector<double> coeff2;
    vector<double> unique_coeff2;
    ChromCorrector(){
        iscorr1=false;
        iscorr2=false;
    }
};



#endif