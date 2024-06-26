from .cppcomponent cimport *

cdef extern from "../utils/cppast.h"nogil:


    cdef enum:
        ADD,    SUB,    MUL,    DIV,    POW,    MOD,    FLOOR,
        ABS,       SQRT,       SIN,    COS,    SINH,   COSH, 
        EXP,       DIM,       MAX,       MIN,
        NUMBER,    PROPERTY,    ID,     REFER,
        DOT,       SLICE,       COMMA,
        DELAY,   ASSIGN,      
        VAR,      CONSTRAINT,   OPTIMIZATION,
        LPAREN,    RPAREN,    LBRA,    RBRA,
        POSITION,  INDEX,
        NEXTLINE,     JOINLINE,      END


    cdef cppclass AST nogil:
        int token
        double value
        bint reducible
        bint owner
        AST()except +
        double calc()except+
        int simplify()except +


    cdef cppclass Number(AST)nogil:
        Number(double number0)except +
        double calc()except+


    cdef cppclass Node(AST)nogil:
        AST* left
        AST* right
        Node(AST* left0, int token0, AST* right0)except +
        double calc()except +
        int simplify()except +




    cdef cppclass Property(AST)nogil:
        int  index, position
        double* data
        Property(int token0, int position0, int index0,  vector[CppComponent*]& line,double* GLB)except +
        double calc()except+
        int simplify()except +

    cdef cppclass Var(AST)nogil:
        Property* left
        AST* right
        double lb
        double ub
        double step
        bint independent
        Var(Property* left0,AST* right0)except +
        Var(Property* left0, double lb0, double ub0, double step0)except +
        double calc()except+
        int simplify()except +


    cdef cppclass Refer(AST)nogil:
        Identity* expr
        string name
        Refer(string& name0, Identity*& ident)except +

        double calc()except+

        int simplify()except +

    cdef cppclass Identity(AST)nogil:
        AST* expr
        string name
        size_t refcount
        bint delayed
        Identity(string& name0,const bint &delayed0, AST* expr0)except +
        double calc()except+
        int simplify()except +
        
            


    cdef cppclass MonoFunction(AST)nogil:
        AST* arg
        MonoFunction(int token0, AST* arg0)except +
                
        double calc()except+
        int simplify()except +

    cdef cppclass BiFunction(AST)nogil:
        AST* arg1
        AST* arg2
        BiFunction(AST* arg10, int token0, AST* arg20)except +
        double calc()except+
        int simplify()except +



        
    cdef cppclass RangeFunction(AST)nogil:
        int     start, end, kind
        vector[double*]    database
        RangeFunction(int token0, int index, int start0, int end0, int kind0, vector[CppComponent*]& line)except +
        double calc()except+
        int simplify()except +
