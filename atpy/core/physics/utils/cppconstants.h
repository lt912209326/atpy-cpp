
#ifndef _CPPCONSTANT_H_
#define _CPPCONSTANT_H_
/* 
code style
macro                               :   ABC_NAME()
const and enum member               :   k_Keywords_Name
class and struct                    :   AbcName
general and class public  method    :   abc_name()
class private and protect method    :   _abs_name()
class public member                 :   abc_name
class private and protect member    :   _abc_name
*/

// #include <iostream>
#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include <cstring>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <map>

// using namespace std;
using std::cin;
using std::cout;
using std::endl;
// 数据类型
using std::vector;
using std::unordered_map;
using std::map;
using std::set;
using std::string;
using std::complex;
// 数学函数
using std::cos;
using std::acos;
using std::cosh;
using std::sin;
using std::asin;
using std::sinh;
using std::sqrt;
using std::fabs;
using std::conj;
using std::isnan;

using std::exp;
using std::memset;
using std::memcpy;
using std::fill;

using std::to_string;

#define sqr(x) ( (x)*(x) )


inline int SIGN(double x){
    if (x<0.){
        return -1;
    }
    else if(x>0.)
    {
        return 1;
    }
    else
    {
        return 0;
    }
};

#define M_PI		3.14159265358979323846
#define __SIZEOF_DOUBLE__ 8

// macro functions

#define NARG(...)       NARG_VCAT_(__VA_ARGS__,NARG_RSEQ_())
#define NARG_VCAT_(...) NARG_NSEQ_(__VA_ARGS__)
#define NARG_NSEQ_( \
     _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
    _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
    _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
    _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
    _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
    _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
    _61,_62,_63,  N, ...) N

#define NARG_RSEQ_() \
    63,62,61,60,                   \
    59,58,57,56,55,54,53,52,51,50, \
    49,48,47,46,45,44,43,42,41,40, \
    39,38,37,36,35,34,33,32,31,30, \
    29,28,27,26,25,24,23,22,21,20, \
    19,18,17,16,15,14,13,12,11,10, \
     9, 8, 7, 6, 5, 4, 3, 2, 1, 0


#define MKNAME(a,b) a##b

const double PI=M_PI;
const double INV_PI=1/M_PI;
const double PIx2=2.0*M_PI;
const double Cq=3.83e-13;
const double k_EPS=1e-8;

typedef unsigned short  Uindex16;
typedef unsigned short  Index32;
typedef unsigned int    Uint32;
typedef long long       Int64 ;
typedef double          Float64;



union DataType{
    void*   p08;
    char    c08;
    Uint32  ui32;
    Int64   i64;
    double  f64;
};

struct ValueBase{
    Uint32      type;
    DataType    value;
};

union AnyValue{
    double Float;
    int64_t Integer;
    uint64_t UInt;
    char Charactor;
    char String[8];
    void *ptr;

};


const double ZERO66[6][6]= {{0. ,0. ,0. ,0. ,0. ,0. },
                            {0. ,0. ,0. ,0. ,0. ,0. },
                            {0. ,0. ,0. ,0. ,0. ,0. },
                            {0. ,0. ,0. ,0. ,0. ,0. },
                            {0. ,0. ,0. ,0. ,0. ,0. },
                            {0. ,0. ,0. ,0. ,0. ,0. }};

const double EYE66[6][6]={{1. ,0. ,0. ,0. ,0. ,0. },
                          {0. ,1. ,0. ,0. ,0. ,0. },
                          {0. ,0. ,1. ,0. ,0. ,0. },
                          {0. ,0. ,0. ,1. ,0. ,0. },
                          {0. ,0. ,0. ,0. ,1. ,0. },
                          {0. ,0. ,0. ,0. ,0. ,1. }};

const double EYE67[6][7]={{1. ,0. ,0. ,0. ,0. ,0. ,0.},
                          {0. ,1. ,0. ,0. ,0. ,0. ,0.},
                          {0. ,0. ,1. ,0. ,0. ,0. ,0.},
                          {0. ,0. ,0. ,1. ,0. ,0. ,0.},
                          {0. ,0. ,0. ,0. ,1. ,0. ,0.},
                          {0. ,0. ,0. ,0. ,0. ,1. ,0.}};
enum Particle{
    ELECTRON,
    POSITRON,
    PROTON
};

        
//#   ELEMENTS CODE                   # 
enum {
    MARKER         = 0,
    DRIFT          = MARKER           +   1,
    DIPOLE         = DRIFT            +   1,
    QUADRUPOLE     = DIPOLE           +   2,
    SEXTUPOLE      = QUADRUPOLE       +   2,
    OCTUPOLE       = SEXTUPOLE        +   2,
    DECAPOLE       = OCTUPOLE         +   2,
    COMBINED       = DECAPOLE         +   2,
    TUNING         = COMBINED         +   2,
    EXACTDRIFT     = TUNING             + 2,
    GIRDER,
    RFCAVITY,
    WIGGLER,
    CRAB
};

enum {
    KWD, TWS, LOC, GLB, DVTs
};

//#   KEYWORD CODE                    #
enum KeywordToken: uint16_t
{
    L       ,ANGLE      ,K1         ,K2         ,K3         ,
    K4      ,E1         ,E2         ,NSLICE     ,TILT       ,
    DNUX    ,DNUY       ,
    BETAX1     ,ALPHAX1    ,BETAY1     ,ALPHAY1 ,ETAX1      ,ETAPX1     ,
    BETAX2     ,ALPHAX2    ,BETAY2  ,ALPHAY2    ,ETAX2      ,ETAPX2     ,
    KWD_NUM
};

const map<int,string> KEYWORDS_DICT={{L,"l"}, {ANGLE,"angle"}, {K1,"k1"}, {K2,"k2"}, {K3,"k3"}, 
                        {K4,"k4"}, {E1,"e1"}, {E2,"e2"}, {NSLICE,"nslice"}, {TILT,"tilt"},
                        {DNUX, "dnux"}, {DNUY, "dnuy"}, 
                        {BETAX1, "betax1"}, {ALPHAX1, "alphax1"}, {BETAY1, "betay1"}, {ALPHAY1, "alphay1"}, {ETAX1,"etax1"}, {ETAPX1,"etapx1"},
                        {BETAX2, "betax2"}, {ALPHAX2, "alphax2"}, {BETAY2, "betay2"}, {ALPHAY2, "alphay2"}, {ETAX2,"etax2"}, {ETAPX2,"etapx2"},
                        // {KWD_NUM, "kwd_num"}
                        };

//#   TWISS CODE                      #
enum TwissToken : uint16_t
{
    COORD      ,BETAX       ,ALPHAX     ,GAMMAX     ,BETAY      ,
    ALPHAY     ,GAMMAY      ,ETAX       ,ETAPX      ,ETAY       ,ETAPY      ,NUX        ,NUY         ,
    R1          ,R2         ,R3         ,R4        ,TWSMODE    ,H0        ,
    DCHROMX     ,DCHROMY    ,CHROMX      ,CHROMY     ,
    LH11001     ,LH00111     ,LH20001     ,LH00201     ,LH10110     ,LH21000     ,
    LH30000     ,LH10200     ,LH10002     ,LH10020,
    LH22000,     LH11110,        LH00220,
    LH31000,     LH40000,         LH20110,        LH11200,       LH20020,        LH20200,
    LH00310,     LH00400,
    TWS_NUM
};


const map<int,string> TWISS_DICT={{COORD,"coord"}, {BETAX,"betax"},  {ALPHAX,"alphax"},  {GAMMAX,"gammax"}, {BETAY,"betay"},  
                                {ALPHAY,"alphay"},  {GAMMAY,"gammay"}, {ETAX,"etax"},  {ETAPX,"etapx"}, {ETAY,"etay"}, {ETAPY,"etapy"}, {NUX,"nux"}, {NUY,"nuy"},  
                                {R1,"R1"}, {R2,"R2"}, {R3,"R3"},{R4,"R4"}, {TWSMODE,"mode"},  {H0, "H0"},
                                {DCHROMX,"dchromx"},  {DCHROMY,"dchromy"},  {CHROMX,"chromx"},  {CHROMY,"chromy"}, 
                                {LH11001, "lh11001"},  {LH00111, "lh00111"},  {LH20001, "lh20001"},  {LH00201, "lh00201"}, {LH10110, "lh10110"}, {LH21000, "lh21000"},
                                {LH30000, "lh30000"},  {LH10200, "lh10200"},  {LH10002, "lh10002"},  {LH10020, "lh10020"},
                                {LH22000, "lh22000"}, {LH11110, "lh11110"}, {LH00220, "lh00220"}, 
                                {LH31000, "lh31000"}, {LH40000, "lh40000"}, {LH20110, "lh20110"}, {LH11200, "lh11200"}, {LH20020, "lh20020"}, {LH20200, "lh20200"},
                                {LH00310, "lh00310"}, {LH00400, "lh00400"}, 
                                // {TWS_NUM, "tws_num"}  
                                };

// LOCAL CODE
enum LocalToken: uint16_t
{
    R11         ,R12        ,R13        ,R14        ,R15        ,R16        ,
    R21         ,R22        ,R23        ,R24        ,R25        ,R26        ,
    R31         ,R32        ,R33        ,R34        ,R35        ,R36        ,
    R41         ,R42        ,R43        ,R44        ,R45        ,R46        ,
    R51         ,R52        ,R53        ,R54        ,R55        ,R56        ,
    R61         ,R62        ,R63        ,R64        ,R65        ,R66        ,
    CODX        ,CODPX      ,CODY       ,CODPY      ,CODZ       ,CODDP      ,
    LOCAL_NUX   ,LOCAL_NUY,
    S           ,GX         ,GY         ,GZ         ,THETAX     ,THETAY     ,
    THETAZ      ,DX         ,DY         ,DZ         ,ROTATE1    ,ROTATE2    ,
    ROTATE3     ,AX         ,AY         ,
    LOC_NUM
};

enum ValueType:uint16_t
{
    kPOINT=0,
    kCHAR,
    kINT,
    kDOUBLE
};


const map<int,string> LOCAL_DICT=   {{R11, "R11"}, {R12, "R12"}, {R13, "R13"}, {R14, "R14"}, {R15, "R15"}, {R16, "R16"},
                                    {R21, "R21"}, {R22, "R22"}, {R23, "R23"}, {R24, "R24"}, {R25, "R25"}, {R26, "R26"},
                                    {R31, "R31"}, {R32, "R32"}, {R33, "R33"}, {R34, "R34"}, {R35, "R35"}, {R36, "R36"},
                                    {R41, "R41"}, {R42, "R42"}, {R43, "R43"}, {R44, "R44"}, {R45, "R45"}, {R46, "R46"},
                                    {R51, "R51"}, {R52, "R52"}, {R53, "R53"}, {R54, "R54"}, {R55, "R55"}, {R56, "R56"},
                                    {R61, "R61"}, {R62, "R62"}, {R63, "R63"}, {R64, "R64"}, {R65, "R65"}, {R66, "R66"},
                                    {CODX, "CODx"}, {CODPX, "CODpx"}, {CODY, "CODy"}, {CODPY, "CODpy"}, {CODZ, "CODz"}, {CODDP, "CODdp"},
                                    {LOCAL_NUX, "local_nux"}, {LOCAL_NUY, "local_nuy"},
                                    {S, "s"},  {GX, "Gx"},  {GY, "Gy"},  {GZ, "Gz"},  {THETAX, "thetax"},  {THETAY, "thetay"},
                                    {THETAZ, "thetaz"},  {DX, "Dx"},  {DY, "Dy"},  {DZ, "Dz"},  {ROTATE1, "rotate1"},  {ROTATE2, "rotate2"},
                                    {ROTATE3, "rotate3"},  {AX, "Ax"},  {AY, "Ay"},
                                    // {LOC_NUM,"loc_num"} 
                                    };
//#   GLOBAL PARAMETER CODE           #
enum GlobalToken: uint16_t
{
    MASS0           ,GAMMA          ,ENERGY         ,EMITX          ,
    NATURE_CHROMX   ,NATURE_CHROMY      ,TOTAL_K2L      ,
    CIRCUMFERENCE   ,RI1            ,RI2            ,RI3        ,RI3A       ,RI4    ,RI5             ,
    U0             ,ALPHAC         ,DAMP_FACTOR    ,
    JX             ,JY              ,JZ             ,ESPREAD   ,SPIN            ,
    TAUX       ,TAUY       ,TAUZ   ,
    QX             ,QY,                 DQX             ,DQY,     D2QX,          D2QY, D3QX,         D3QY,        D4QX,         D4QY,    
    H10010,        H10100,      H11001,        H00111,  
    //  Order==3     
    H21000,        H30000,     H10110,     H10020, H10200,   
    H20001,      H00201,     H10002,       
    DNUX_DJX,      DNUX_DJY,        DNUY_DJY,
    //  Order==4
        H22000,     H11110,        H00220,
    H31000,        H40000,         H20110,        H11200,       H20020,        H20200,
    H00310,         H00400,                DA    ,     DA_SIGMA ,
    DETAX ,         DETAPX,        DBETAX,        DBETAY,       DALPHAX,        DALPHAY,       
    DDETAX,     DDBETAX,    DDBETAY,    WX,      WY,  
    LOW_QX             ,LOW_QY,            HIGH_QX             ,HIGH_QY,        SUM_SQR_QX,     SUM_SQR_QY,
    INV_TAU,
    GLB_NUM
};


const map<int,string> GLOBALS_DICT= {{MASS0,"mass0"},{GAMMA,"gamma"},{ENERGY,"energy"},{EMITX,"emitx"}, 
                        {NATURE_CHROMX,"nature_chromx"}, {NATURE_CHROMY,"nature_chromy"}, {TOTAL_K2L,"total_K2L"},
                        {CIRCUMFERENCE,"circumference"}, {RI1,"RI1"},{RI2,"RI2"},{RI3,"RI3"},{RI3A,"RI3a"},{RI4,"RI4"},{RI5,"RI5"},
                        {U0,"U0"},{ALPHAC,"alphac"},{DAMP_FACTOR,"damp_factor"},
                        {JX,"Jx"}, {JY,"Jy"},{JZ,"Jz"},{ESPREAD,"e_spread"}, {SPIN,"spin"},
                        {TAUX,"taux"},{TAUY,"tauy"},{TAUZ,"tauz"},
                            {QX,"Qx"},{QY,"Qy"}, {DQX,"dQx"},{DQY,"dQy"}, {D2QX,"d2Qx"}, {D2QY,"d2Qy"}, {D3QX,"d3Qx"}, {D3QY,"d3Qy"}, {D4QX,"d4Qx"},  {D4QY,"d4Qy"}, 
                        {H10010,"H10010"},  {H10100,"H10100"}, {H11001,"H11001"},  {H00111,"H00111"},
                        {H21000,"H21000"},  {H30000,"H30000"},  {H10200,"H10200"},  {H10020,"H10020"}, {H10002,"H10002"}, 
                        {H20001,"H20001"},  {H00201,"H00201"},  {H10110,"H10110"},  
                        {DNUX_DJX,"dnux_dJx"},  {DNUX_DJY,"dnux_dJy"},  {DNUY_DJY,"dnuy_dJy"},
                        {H22000,"H22000"}, {H11110,"H11110"}, {H00220,"H00220"}, 
                        {H31000,"H31000"}, {H40000,"H40000"}, {H20110,"H20110"},   {H11200,"H11200"}, {H20020,"H20020"}, {H20200,"H20200"},
                        {H00310,"H00310"}, {H00400,"H00400"}, {DA, "DA"}, {DA_SIGMA, "DA_SIGMA"},
                        {DETAX, "detax"} , {DETAPX, "detapx"}, {DBETAX, "dbetax"}, {DBETAY, "dbetay"}, {DALPHAX, "dalphax"}, {DALPHAY, "dalphay"}, 
                        {DDETAX, "ddetax"}, {DDBETAX, "ddbetax"}, {DDBETAY, "ddbetay"}, {WX, "Wx"}, {WY, "Wy"}, 
                        {LOW_QX, "low_Qx"}  ,{LOW_QY, "low_Qy"},  {HIGH_QX, "high_Qx"}  ,{HIGH_QY, "high_Qy"}, 
                        {SUM_SQR_QX, "sum_sqr_Qx"},     {SUM_SQR_QY, "sum_sqr_Qy"}, {INV_TAU, "inv_tau"},
                        // { GLB_NUM, "glb_num"} 
                        };
enum {
    X           ,PX         ,Y          ,PY         ,Z          ,DP         ,LOSS           ,   
    LOSSTURN    ,LOSSPOS    ,
    TRK_NUM
};

enum {
    CHROMATIC_TERMS,    DRIVING_TERMS,  DA_TRACKING_TERMS,  MA_TRACKING_TERMS, MONITOR_OFF_MOMENTUM_TERMS, OFF_MOMENTUM_SUM_TERMS
};


//#   CONSTRAINT CONST                #
enum {
    ZERO    ,ONE    ,TWO    ,THREE  ,FOUR   ,
    FIVE    ,SIX    ,SEVEN  ,EIGHT  ,NINE   ,
    TEN
};
enum {
    LOWER,
    UPPER,
    BOTH,
    EQUIV
 };           

#endif