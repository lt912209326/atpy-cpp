from libcpp.vector cimport vector
from libcpp.algorithm cimport count as cpp_count, fill
from libcpp.map cimport map as c_map
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.complex cimport complex
from libcpp.set cimport set as cpp_set
# from cppstructures cimport *


#
cdef extern from "physics/utils/cppconstants.h":
    # double sqr(double x)
    cdef:
        const double PI
        const double INV_PI
        const double PIx2
        const double Cq
        const double ZERO66[6][6]
        const double EYE66[6][6]
        const double EYE67[6][7]
        
    cdef enum:
        ELECTRON,
        POSITRON,
        PROTON
    

            
    ##   ELEMENTS CODE                   # 
    cdef enum:
        MARKER         = 0
        DRIFT          = MARKER           +   1
        DIPOLE         = DRIFT            +   1
        QUADRUPOLE     = DIPOLE           +   2
        SEXTUPOLE      = QUADRUPOLE       +   2
        OCTUPOLE       = SEXTUPOLE        +   2
        DECAPOLE       = OCTUPOLE         +   2
        COMBINED       = DECAPOLE         +   2
        TUNING         = COMBINED         +   2
        EXACTDRIFT     = TUNING             + 2
        GIRDER
        RFCAVITY
        WIGGLER
        CRAB
    
    cdef enum:
        KWD, TWS, LOC, GLB, DVTs


    ##   KEYWORD CODE                    #
    cdef enum:
        L       ,ANGLE      ,K1         ,K2         ,K3         ,
        K4      ,E1         ,E2         ,NSLICE     ,TILT       ,
        DNUX    ,DNUY       ,BETAX1     ,ALPHAX1    ,BETAY1     ,
        ALPHAY1 ,ETAX1      ,ETAPX1     ,BETAX2     ,ALPHAX2    ,
        BETAY2  ,ALPHAY2    ,ETAX2      ,ETAPX2     ,
        KWD_NUM
    

    cdef const c_map[int,string] KEYWORDS_DICT
    ##   TWISS CODE                      #
    cdef enum:
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

    cdef const c_map[int,string] TWISS_DICT

    # LOCAL CODE
    cdef enum:
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
        

    cdef const c_map[int,string] LOCAL_DICT


    ##   GLOBAL PARAMETER CODE           #
    cdef enum:
        MASS0           ,GAMMA          ,ENERGY         ,EMITX          ,NATURE_CHROMX     ,
        NATURE_CHROMY      ,TOTAL_K2L      ,
        CIRCUMFERENCE   ,RI1            ,RI2            ,RI3        ,RI3A       ,RI4    ,
        RI5             ,U0             ,ALPHAC         ,DAMP_FACTOR    ,JX             ,
        JY              ,JZ             ,ESPREAD   ,
        SPIN            ,TAUX       ,TAUY       ,TAUZ   ,
        QX             ,QY,                 DQX             ,DQY,     D2QX,          D2QY, D3QX,         D3QY,        D4QX,         D4QY,    
        H10010,        H10100,      H11001,        H00111,  
        #  Order==3     
        H21000,        H30000,     H10110,     H10020, H10200,   H20001,      H00201,     H10002,       
        DNUX_DJX,      DNUX_DJY,        DNUY_DJY,
        #  Order==4
        H22000,     H11110,        H00220,
        H31000,        H40000,         H20110,        H11200,       H20020,        H20200,
        H00310,         H00400,                DA    ,     DA_SIGMA ,
        DETAX ,         DETAPX,        DBETAX,        DBETAY,       DALPHAX,        DALPHAY,       DDETAX,     DDBETAX,    DDBETAY,    WX,      WY,  
        LOW_QX             ,LOW_QY,            HIGH_QX             ,HIGH_QY,        SUM_SQR_QX,     SUM_SQR_QY,
        INV_TAU,
        GLB_NUM
    


    cdef const c_map[int,string] GLOBALS_DICT

    cdef enum:
        X           ,PX         ,Y          ,PY         ,Z          ,DP         ,LOSS           ,   
        LOSSTURN    ,LOSSPOS    ,
        TRK_NUM
    

    cdef enum :
        CHROMATIC_TERMS,    DRIVING_TERMS,  DA_TRACKING_TERMS,  MA_TRACKING_TERMS, MONITOR_OFF_MOMENTUM_TERMS, OFF_MOMENTUM_SUM_TERMS
    

    ##   CONSTRAINT CONST                #
    cdef enum:
        ZERO    ,ONE    ,TWO    ,THREE  ,FOUR   ,
        FIVE    ,SIX    ,SEVEN  ,EIGHT  ,NINE   ,
        TEN
    

    cdef enum:
        LOWER,
        UPPER,
        BOTH,
        EQUIV

    cdef cppclass matrix[T]:
        int nrow
        int ncol
        T* data

        matrix ()
        matrix (size_t nRow, size_t nCol)
        matrix (size_t nRow, size_t nCol, const T& e)
        matrix (const matrix[T]& m)

        matrix[T]& operator= (const T& e)
        matrix[T]& operator= (const matrix[T]& mat)

        size_t size ()
        T& operator() (int i, int j) const
        T& operator[] (int i) const





cdef dict KWD_INDEX
cdef dict INDEX_KWD

cdef dict TWS_INDEX
cdef dict INDEX_TWS

cdef dict LOC_INDEX
cdef dict INDEX_LOC

cdef dict GLB_INDEX
cdef dict INDEX_GLB

cdef dict ALL_INDEX

cdef dict ELEM_INDEX
cdef dict INDEX_ELEM


cdef list POSITION_DEPEND_NAMES
cdef dict DEFAULT_ELEM_KARGS

