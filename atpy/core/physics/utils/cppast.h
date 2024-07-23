#ifndef _CPPPARSER_H_
#define _CPPPARSER_H_

#include "cppcomponent.h"
#include <exception>

enum{
    ADD,    SUB,    MUL,    DIV,    POW,    MOD,    FLOOR,
    ABS,       SQRT,        SIN,    COS,    SINH,   COSH, 
    EXP,       
    DIM,       MAX,       MIN,       MAXABS,       MINABS, SUM,
    NUMBER,    PROPERTY,    ID,     REFER,
    DOT,       COMMA,
    DELAY,   ASSIGN,      
    VAR,      CONSTRAINT,   OPTIMIZATION,
    LPAREN,    RPAREN,    LBRA,    RBRA,
    POSITION,  INDEX,
    NEWLINE,     JOINLINE,      END
};
    




class AST{
    public:
    int token;
    double value;
    bool reducible;
    bool owner;
    AST(){
        token=0;
        value=0;
        reducible=false;
        owner=true;
    }
    virtual double calc()=0;
    virtual int simplify(){return 0;}
    // virtual void dealloc(){}
    virtual ~AST(){}
};




class Number:public AST
{
    public:
    Number(double number0){
        token=NUMBER;
        value=number0;
        reducible=false;
    }

    double calc()override{
        return value;
    }
};




class Node: public AST
{
    public:
    AST* left;
    AST* right;
    Node(AST* left0, int token0, AST* right0){
        left=left0;
        token=token0;
        right=right0;
        reducible=false;
    }
    
    double calc()override{
        switch(token){
            case ADD:
                value = left->calc() + right->calc();
                break;
            case SUB:
                value = left->calc() - right->calc();
                break;
            case MUL:
                value = left->calc()*right->calc();
                break;
            case DIV:
                value = left->calc()/right->calc();
                break;
            case POW:
                value = pow(left->calc(),right->calc() );
                break;
            case MOD:
                value = fmod(left->calc(), right->calc() );
                break;
            case FLOOR:
                value = floor( left->calc()/right->calc() );
                break;
            default:
                break;
        }  
        return value; 
    }

    int simplify()override{
        double tmp=0;
        left->simplify();
        if(left->reducible){
            tmp=left->calc();
            delete left;
            left=new Number(tmp);
        }
        right->simplify();
        if(right->reducible){
            tmp=right->calc();
            delete right;
            right=new Number(tmp);
        }
        if( left->token==NUMBER && right->token==NUMBER ){
            reducible=true;
        }
        return 0;
    }

    ~Node(){
        delete left;
        delete right;
        right=nullptr;
        left=nullptr;

    }
};



class Property:public AST{

    public:
    int  index, position;
    double* data;
    Property(int token0, int position0, int index0,  vector<CppComponent*>& line,double* GLB0){
        token = token0;
        position = position0;
        index = index0;
        if(token==KWD){
            if(count(line[position]->elem->keywords.begin(), line[position]->elem->keywords.end(), index) ){
                data=&(line[position]->elem->values.at(index) );
            }
            else{
                throw std::invalid_argument("Property::Property:"+line[position0]->elem->name+ "has no such keyword!");
            }
        }
        else if(token==TWS){
            data=&(line[position]->tws(-1,index) );
        }
        else if(token==LOC){
            data=&(line[position]->local[index] );
        }
        else if(token==GLB){
            data=&(GLB0[index] );
        }
        // value = database0;
        reducible=false;
    }
    
    double calc()override{
        value=data[0];
        return value;
    }

    int simplify(){
        return 0;
    }
    ~Property(){
        data=nullptr;
    }

};

class Var:public AST{
    public:
    Property* left;
    AST* right;
    double lb;
    double ub;
    double step;
    bool independent;
    Var(Property* left0,AST* right0){
        left=left0;
        right=right0;
        lb=-1e30;
        ub=1e30;
        independent=false;
        reducible=false;
    }
    Var(Property* left0, double lb0, double ub0, double step0){
        left=left0;
        lb=lb0;
        ub=ub0;
        step=step0;
        independent=true;
        reducible=false;
    }
    double calc()override{
        left->data[0]=right->calc();
        return 0;
    }
    int simplify()override{
        double tmp=0;
        if(!independent){
            right->simplify();
            if(right->reducible){
                tmp=right->calc();
                delete right;
                right=new Number(tmp);
            }
        }
        return 0;
    }
    ~Var(){
        if(independent){
            delete left;
            left=nullptr;
        }
        else{
            delete left;
            delete right;
            left=nullptr;
            right=nullptr;
        }
    }
};




class Identity:public AST{
    public:
    AST* expr;
    string name;
    size_t refcount;
    bool delayed;
    Identity(string& name0,const bool &delayed0, AST* expr0){
        token=ID;
        name=name0;
        delayed=delayed0;
        expr=expr0;
        reducible=false;
        if(!delayed){
            reducible=true;
        }
        calc();
    }

    double calc()override{
        value=expr->calc();
        return value;
    }

    int simplify()override{
        double tmp=0;
        expr->simplify();
        if(expr->token==NUMBER){
            tmp=expr->calc();
            expr=new Number(tmp);
            reducible=true;
        }
        else if(delayed==false){//delayed is false, then assign is true
            tmp=expr->calc();
            expr=new Number(tmp);
            reducible=true;
        }
        // value=expr->calc();
        return 0;
    }
    ~Identity(){
        delete expr;
        expr=nullptr;
    }
};
    

class Refer:public AST{
    public:
    AST* expr;
    string name;
    Refer(string& name0,Identity* ident){
        token=REFER;
        name=name0;
        expr=ident;
        reducible=false;
    }

    double calc()override{
        value=expr->value;
        return value;
    }

    int simplify()override{
        double tmp=0;
        expr->simplify();
        if(expr->token==NUMBER){
            tmp=expr->calc();
            expr=new Number(tmp);
        }
        reducible=false;
        return 0;
    }
    ~Refer(){
        expr=nullptr;
    }
};


class MonoFunction:public AST{
    public:
    AST* arg;
    MonoFunction(int token0, AST* arg0){
        token = token0;
        arg = arg0;
        reducible=false;
    }
            
    double calc()override{
        switch(token){
            case ABS:
                value=fabs(arg->calc() );
                break;
            case SQRT:
                value=sqrt(arg->calc() );
                break;
            case SIN:
                value=sin(arg->calc() );
                break;
            case COS:
                value=cos(arg->calc() );
                break;
            case SINH:
                value=sinh(arg->calc() );
                break;
            case COSH:
                value=cosh(arg->calc() );
                break;
            case EXP:
                value=exp(arg->calc() );
                break;
            default:
                break;
        }
        return value;
    }
    int simplify()override{
        double tmp;
        arg->simplify();
        if(arg->reducible){
            tmp=arg->calc();
            arg=new Number(tmp);
        }
        if(arg->token==NUMBER){
            reducible=true;
        }
        return 0;
    }
    ~MonoFunction(){
        delete arg;
        arg=nullptr;
    }

};

class BiFunction:public AST{
    public:
    AST* arg1;
    AST* arg2;
    BiFunction(AST* arg10, int token0, AST* arg20){
        arg1 = arg10;
        token = token0;
        arg2 = arg20;
        reducible=false;
    }
        
    double calc()override{
        switch(token){
            case DIM:
            value=fdim(arg1->calc(), arg2->calc() );
            break;
            default:
            break;
        }
        return value;
    }

    int simplify()override{
        double tmp=0;
        arg1->simplify();
        if(arg1->reducible){
            tmp=arg1->calc();
            delete arg1;
            arg1=new Number(tmp);
        }
        arg2->simplify();
        if(arg2->reducible){
            tmp=arg2->calc();
            delete arg2;
            arg2=new Number(tmp);
        }
        if( arg1->token==NUMBER && arg2->token==NUMBER ){
            reducible=true;
        }
        return 0;
    }
    ~BiFunction(){
        delete arg1;
        delete arg2;
        arg1=nullptr;
        arg2=nullptr;
    }

};


    
class RangeFunction:public AST{
    public:
    int     start, end, kind,index;
    vector<double*>    database;
    RangeFunction(int token0,int index0, int start0, int end0, int kind0, vector<CppComponent*>& line){
        token=token0;
        start = start0;
        index = index0;
        end = end0;
        kind=kind0;
        reducible=false;
        database={};
        switch(kind){
            case KWD:
                for(size_t i=start;i<end+1;i++){
                    if(count(line[i]->elem->keywords.begin(),line[i]->elem->keywords.end(),index)!=0){
                        database.emplace_back(&(line[i]->elem->values[index]));
                    }
                }
                break;
            case TWS:
                for(size_t i=start;i<end+1;i++){
                    database.emplace_back(&(line[i]->tws( int(line[i]->elem->nslice/2),index) ));
                }
                break;
            case LOC:
                for(size_t i=start;i<end+1;i++){
                    database.emplace_back(&(line[i]->local[index]));
                }
                break;
            default:
                
                throw std::invalid_argument("RangeFunction::RangeFunction: no match kind in RangeFunction!");
        }
    }
    
    
    double calc()override{
        double value0=0;
        if(token==MIN){
            value0=1e40;
            for(auto it : database){
                value0=fmin(it[0],value0);
            }
        }
        else if(token==MAX){
            value0=-1e40;
            for(auto it : database){
                value0=fmax(it[0],value0);
            }
        }
        else if(token==MINABS){
            value0=1e40;
            for(auto it : database){
                value0=fmin(abs(it[0]),value0);
            }
        }
        else if(token==MAXABS){
            value0=-1e40;
            for(auto it : database){
                value0=fmax(abs(it[0]),value0);
            }
        }
        value=value0;
        return value;
    }
    int simplify(){
        reducible=false;
        return 0;
    }

    ~RangeFunction(){
        for(auto &iter:database){
            iter=nullptr;
        }
    }
};





#endif
