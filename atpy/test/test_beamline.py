#-*- HOA 弧区 -*-
from atpy import *
MK2        =     Marker('MK2'     )
LALPHAD4   =      Drift('LALPHAD4',l    =1.611414e+00  )
NEWLALPHAQ3 = Quadrupole('NEWLALPHAQ3',l    =2.000000e-01  ,k1   =-2.320511e+00  )
LALPHAD3   =      Drift('LALPHAD3',l    =1.997493e-01  )
NEWLALPHAQ2 = Quadrupole('NEWLALPHAQ2',l    =3.000000e-01  ,k1   =4.154137e+00  )
LALPHAD2   =      Drift('LALPHAD2',l    =3.432117e-01  )
NEWLALPHAQ1 = Quadrupole('NEWLALPHAQ1',l    =2.000000e-01  ,k1   =-2.524451e+00  )
LALPHAD1   =      Drift('LALPHAD1',l    =2.962478e-01  )
ARCB3      =     Dipole('ARCB3'   ,l    =1.041136e+00  ,angle=9.168448e-02  ,e1   =5.000000e-01  ,e2   =5.000000e-01  )
ARAD2      =      Drift('ARAD2'   ,l    =1.128941e-01  )
ARAQ1      = Quadrupole('ARAQ1'   ,l    =2.000000e-01  ,k1   =-2.056375e+00  )
ARAD1      =      Drift('ARAD1'   ,l    =1.601213e-01  )
CHROMSY2   =  Sextupole('CHROMSY2',l    =1.000000e-01  ,k2   =5.233701e+01,nslice=4  )
PMD22L     =      Drift('PMD22L'  ,l    =1.029104e+00  )
PMQ3       = Quadrupole('PMQ3'    ,l    =2.000000e-01  ,k1   =2.088028e+00  )
PMD3       =      Drift('PMD3'    ,l    =1.000000e-01  )
CHROMSX1   =  Sextupole('CHROMSX1',l    =1.000000e-01  ,k2   =5.648354e+01,nslice=4  )
PMD2       =      Drift('PMD2'    ,l    =1.000000e-01  )
PMQ2       = Quadrupole('PMQ2'    ,l    =2.000000e-01  ,k1   =2.521066e+00  )
PMD21R     =      Drift('PMD21R'  ,l    =8.500382e-01  )
CHROMSY1   =  Sextupole('CHROMSY1',l    =1.000000e-01  ,k2   =5.421104e+00,nslice=4  )
PMD21L     =      Drift('PMD21L'  ,l    =1.570584e-01  )
PMQ1       = Quadrupole('PMQ1'    ,l    =2.000000e-01  ,k1   =-2.259560e+00  )
PMD1       =      Drift('PMD1'    ,l    =1.466716e-01  )
ARCB2      =     Dipole('ARCB2'   ,l    =8.989812e-01  ,angle=9.153650e-02  ,e1   =5.000000e-01  ,e2   =5.000000e-01  )
LAD5       =      Drift('LAD5'    ,l    =3.159382e-01  )
LADQ4      = Quadrupole('LADQ4'   ,l    =2.000000e-01  ,k1   =-1.724579e+00  )
LAD4       =      Drift('LAD4'    ,l    =4.183039e-01  )
LAFQ2      = Quadrupole('LAFQ2'   ,l    =3.000000e-01  ,k1   =4.380521e+00  )
LAD3       =      Drift('LAD3'    ,l    =5.386876e-01  )
LADQ3      = Quadrupole('LADQ3'   ,l    =2.000000e-01  ,k1   =-3.285513e+00  )
LAD01      =      Drift('LAD01'   ,l    =1.000000e-01  )
ARCB1      =     Dipole('ARCB1'   ,l    =8.317428e-01  ,angle=1.048069e-01  ,e1   =5.000000e-01  ,e2   =5.000000e-01  )
LADQ2      = Quadrupole('LADQ2'   ,l    =2.000000e-01  ,k1   =-2.017330e+00  )
LAD2       =      Drift('LAD2'    ,l    =5.428899e-01  )
LAFQ1      = Quadrupole('LAFQ1'   ,l    =3.000000e-01  ,k1   =4.808714e+00  )
LAD1       =      Drift('LAD1'    ,l    =4.262762e-01  )
LADQ1      = Quadrupole('LADQ1'   ,l    =2.000000e-01  ,k1   =-2.720660e+00  )
ARCB0      =     Dipole('ARCB0'   ,l    =8.688655e-01  ,angle=5.226279e-02  ,e1   =5.000000e-01  ,e2   =5.000000e-01  )

ARC        =    Line('ARC'     ,MK2       ,LALPHAD4  ,NEWLALPHAQ3,LALPHAD3  ,NEWLALPHAQ2,LALPHAD2  ,NEWLALPHAQ1,LALPHAD1  ,
                        ARCB3     ,ARAD2     ,ARAQ1     ,ARAD1     ,CHROMSY2  ,CHROMSY2  ,PMD22L    ,PMQ3      ,PMD3      ,CHROMSX1  ,
                        CHROMSX1  ,PMD2      ,PMQ2      ,PMD21R    ,CHROMSY1  ,CHROMSY1  ,PMD21L    ,PMQ1      ,PMD1      ,ARCB2     ,LAD5      ,
                     LADQ4     ,LAD4      ,LAFQ2     ,LAD3      ,LADQ3     ,LAD01     ,ARCB1     ,LAD01     ,LADQ2     ,LAD2      ,LAFQ1     ,
                     LAD1      ,LADQ1     ,LAD01     ,ARCB0     ,LAD01     ,LADQ1     ,
                     LAD1      ,LAFQ1     ,LAD2      ,LADQ2     ,LAD01     ,ARCB1     ,LAD01     ,LADQ3     ,LAD3      ,LAFQ2     ,LAD4      ,LADQ4     ,
                     LAD5      ,ARCB2     ,PMD1      ,PMQ1      ,PMD21L    ,CHROMSY2  ,CHROMSY2  ,PMD21R    ,PMQ2      ,PMD2      ,CHROMSX1  ,
                     CHROMSX1  ,PMD3      ,PMQ3      ,PMD22L    ,CHROMSY1  ,CHROMSY1  ,ARAD1     ,ARAQ1     ,ARAD2     ,ARCB3     ,
                     LALPHAD1  ,NEWLALPHAQ1,LALPHAD2  ,NEWLALPHAQ2,LALPHAD3  ,NEWLALPHAQ3,LALPHAD4  ,MK2       )
TOTALRING=Line("TOTALRING",10*ARC)

# from atpy import *
# translate(lat)
# from atpy.tools.translate import*



stat=Status(period=True,computedrivingterms=True, nonlineartermonly=True,third_order_chrom=True,dp=2e-4,nperiods=1,printout=True)

tws0 ={"betax": 2.1061, "betay":0.324,"etax":0} 


RING=BeamLine("RING", stat,ARC, **tws0)

token=r"""

    VAR, NAME=CHROMSX1[0].k2,LOWER= 0,  UPPER= 100;
    VAR, NAME=CHROMSY1[0].k2,LOWER=-100,  UPPER= 40;
    VAR, NAME=CHROMSY2[0].k2,LOWER=-100,  UPPER= 40;
    
    
    CONSTRAINT,EXPR:=DIM(ABS(dQx),1)+DIM(ABS(dQy),1);
    CONSTRAINT,EXPR:=DIM(ABS(d2Qx),100)+DIM(ABS(d2Qy),100);
    CONSTRAINT,EXPR:=DIM(ABS(d3Qx),500)+DIM(ABS(d3Qy),500);
    CONSTRAINT,EXPR:=DIM(ABS(dnux_dJx),5000)+DIM(ABS(dnux_dJy),5000) +DIM(ABS(dnuy_dJy),5000);
    
    
    OPTIMIZE,EXPR:=DIM(ABS(d2Qx),10)+DIM(ABS(d2Qy),10);
    OPTIMIZE,EXPR:=DIM(ABS(d3Qx),100)+DIM(ABS(d3Qy),100);
    OPTIMIZE,EXPR:=DIM(ABS(dnux_dJx),1000)+DIM(ABS(dnux_dJy),1000) +DIM(ABS(dnuy_dJy),1000);
    
    # CHROM,AIM_DQX=0,KNOB=CHROMSX1;
    # CHROM,SCCX01A:=-0.084*SCCX01;
    # CHROM,AIM_DQY=0,KNOB=CHROMSY1;
    # CHROM,SCCY01A:=-0.084*SCCY01;
    
    
    
"""

RING.parse(token)

RING.display("VAR")
RING.display("CONSTRAINT")
RING.display("OPTIMIZE")
da=RING.get_DA_area()
import matplotlib.pyplot as plt
fig,ax=plt.subplots(2,1)
ax[0].plot(da[0],da[1],".r")
fig.savefig("da.png",dpi=200)
print(da)
# problem=MyProblem(RING,1e-6)