from atpy.core.elements import*
from atpy.core.beamline import*
mk0=Marker("mk0")
d1=Drift("d1",l=0.5)
d2=Drift("d2",l=0.5)
d3=Drift("d3",l=0.5)

b2=Dipole("b2",l=0.5,angle=0.1,e2=0.8)
b3=Dipole("b3",l=0.5,angle=0.1,e2=0.8)
q3=Quadrupole("q3",l=0.3)
s4=Sextupole("s4",l=0.5,k2=1)
aa=Line("aa",mk0,d1,b3,d3,-2*d2,2*b2,3*q3,-4*s4)
stat=Status()

aa.export("../test_atpy.py","atpy")
aa.export("../test_atpy.sad","sad")
aa.export("../test_atpy.opa","opa")
aa.export("../test_atpy.madx","madx")
tws0={
    "betax":0.09,
    "betay":0.0006
}
# stat=Status(period=False)
line=BeamLine("line",stat,aa)


print(mk0.str("atpy"))
print(d1.str("atpy") )
print(b2.str("atpy") )
print(q3.str("atpy") )
print(s4.str("atpy") )
line.calc()
print(line)