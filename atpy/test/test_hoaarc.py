from atpy.core.elements import*
from atpy.core.beamline import*
from math import pi
sext_slices=30


D1 = Drift('D1', l=2.500000)
D2 = Drift('D2', l=0.240000)
D2A = Drift('D2A', l=0.120000)
D2B = Drift('D2B', l=0.120000)
D3 = Drift('D3', l=0.130000)
D4 = Drift('D4', l=0.345000)
D5 = Drift('D5', l=0.200000)
D6 = Drift('D6', l=0.120000)
D7 = Drift('D7', l=0.120000)
D8 = Drift('D8', l=0.200000)
D9 = Drift('D9', l=0.175000)
D10 = Drift('D10', l=0.120000)
D11 = Drift('D11', l=0.120000)
D12 = Drift('D12', l=0.120000)

Q1 = Quadrupole('Q1', l=0.210000, k1=6.381538)
Q2 = Quadrupole('Q2', l=0.210000, k1=-5.963180)
Q3 = Quadrupole('Q3', l=0.140000, k1=6.684036)

theta1 = 2.476187 * pi / 180
theta2 = 5.017160 * pi / 180
theta3 = -0.315127 * pi / 180

B1 = Dipole('B1', l=0.650000, angle=theta1, k1=-0.303251)
B2 = Dipole('B2', l=0.900000, angle=theta2, k1=-1.371681)
RB = Dipole('RB', l=0.160000, angle=theta3, k1=6.074757)

SF1 = Sextupole('SF1', l=0.160000, k2=312.223000 * 2, nslice=sext_slices)
SD1 = Sextupole('SD1', l=0.160000, k2=-291.244000 * 2, nslice=sext_slices)
SF2 = Sextupole('SF2', l=0.130000, k2=338.954000 * 2, nslice=sext_slices)
SD2 = Sextupole('SD2', l=0.130000, k2=-235.879000 * 2, nslice=sext_slices)
# SF1 = Sextupole('SF1', l=0.160000, k2=350.000000 * 2, nslice=sext_slices)
# SD1 = Sextupole('SD1', l=0.160000, k2=-250.0000 * 2, nslice=sext_slices)
# SF2 = Sextupole('SF2', l=0.130000, k2=400.0000 * 2, nslice=sext_slices)
# SD2 = Sextupole('SD2', l=0.130000, k2=-200.000 * 2, nslice=sext_slices)

half = Line("half",D1, Q1, D2A, D2B, Q2, D3, B1, D4, SD1, D5, Q3, D6, SF1, D7,
        RB, D8, SD2, D9, B2, D9, SD2, D8, RB, D10, SF2, D11, RB, D8, SD2, D9, B2,
        D9, SD2, D8, RB, D12)
cell = Line("cell",half,SF2,-half)



stat=Status(period=True,computedrivingterms=True,energy=2.2e9)

RING = BeamLine("RING",stat,cell)

RING.highorderchromaticity(0.001)
# RING.display("TWS",False)
RING.display("DVTs",False)

