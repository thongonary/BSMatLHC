# This file was automatically created by FeynRules 1.7.160
# Mathematica version: 9.0 for Microsoft Windows (64-bit) (January 25, 2013)
# Date: Tue 12 Nov 2013 12:31:00


from object_library import all_lorentz, Lorentz

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot
try:
   import form_factors as ForFac 
except ImportError:
   pass


UUV1 = Lorentz(name = 'UUV1',
               spins = [ -1, -1, 3 ],
               structure = 'P(3,2) + P(3,3)')

SSS1 = Lorentz(name = 'SSS1',
               spins = [ 1, 1, 1 ],
               structure = '1')

FFS1 = Lorentz(name = 'FFS1',
               spins = [ 2, 2, 1 ],
               structure = 'ProjM(2,1) + ProjP(2,1)')

FFV1 = Lorentz(name = 'FFV1',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,1)')

FFV2 = Lorentz(name = 'FFV2',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1)')

FFV3 = Lorentz(name = 'FFV3',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) - 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV4 = Lorentz(name = 'FFV4',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 2*Gamma(3,2,-1)*ProjP(-1,1)')

FFV5 = Lorentz(name = 'FFV5',
               spins = [ 2, 2, 3 ],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 4*Gamma(3,2,-1)*ProjP(-1,1)')

VVS1 = Lorentz(name = 'VVS1',
               spins = [ 3, 3, 1 ],
               structure = 'Metric(1,2)')

VVV1 = Lorentz(name = 'VVV1',
               spins = [ 3, 3, 3 ],
               structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')

SSSS1 = Lorentz(name = 'SSSS1',
                spins = [ 1, 1, 1, 1 ],
                structure = '1')

FFVS1 = Lorentz(name = 'FFVS1',
                spins = [ 2, 2, 3, 1 ],
                structure = 'P(-1,3)*P(3,4)*Gamma(-1,2,1) - P(-1,3)*P(-1,4)*Gamma(3,2,1)')

FFVV1 = Lorentz(name = 'FFVV1',
                spins = [ 2, 2, 3, 3 ],
                structure = 'P(4,3)*Gamma(3,2,1) - P(-1,3)*Gamma(-1,2,1)*Metric(3,4)')

FFVV2 = Lorentz(name = 'FFVV2',
                spins = [ 2, 2, 3, 3 ],
                structure = 'P(4,3)*Gamma(3,2,1) + P(3,4)*Gamma(4,2,1) - P(-1,3)*Gamma(-1,2,1)*Metric(3,4) - P(-1,4)*Gamma(-1,2,1)*Metric(3,4)')

VVSS1 = Lorentz(name = 'VVSS1',
                spins = [ 3, 3, 1, 1 ],
                structure = 'Metric(1,2)')

VVVV1 = Lorentz(name = 'VVVV1',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')

VVVV2 = Lorentz(name = 'VVVV2',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')

VVVV3 = Lorentz(name = 'VVVV3',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')

VVVV4 = Lorentz(name = 'VVVV4',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')

VVVV5 = Lorentz(name = 'VVVV5',
                spins = [ 3, 3, 3, 3 ],
                structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')

FFVSS1 = Lorentz(name = 'FFVSS1',
                 spins = [ 2, 2, 3, 1, 1 ],
                 structure = 'P(-1,3)*P(3,4)*Gamma(-1,2,1) + P(-1,3)*P(3,5)*Gamma(-1,2,1) - P(-1,3)*P(-1,4)*Gamma(3,2,1) - P(-1,3)*P(-1,5)*Gamma(3,2,1)')

FFVVS1 = Lorentz(name = 'FFVVS1',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'P(4,3)*Gamma(3,2,1) - P(-1,3)*Gamma(-1,2,1)*Metric(3,4)')

FFVVS2 = Lorentz(name = 'FFVVS2',
                 spins = [ 2, 2, 3, 3, 1 ],
                 structure = 'P(4,3)*Gamma(3,2,1) + P(3,4)*Gamma(4,2,1) - P(-1,3)*Gamma(-1,2,1)*Metric(3,4) - P(-1,4)*Gamma(-1,2,1)*Metric(3,4)')

FFVVSS1 = Lorentz(name = 'FFVVSS1',
                  spins = [ 2, 2, 3, 3, 1, 1 ],
                  structure = 'P(4,3)*Gamma(3,2,1) - P(-1,3)*Gamma(-1,2,1)*Metric(3,4)')

FFVVSS2 = Lorentz(name = 'FFVVSS2',
                  spins = [ 2, 2, 3, 3, 1, 1 ],
                  structure = 'P(4,3)*Gamma(3,2,1) + P(3,4)*Gamma(4,2,1) - P(-1,3)*Gamma(-1,2,1)*Metric(3,4) - P(-1,4)*Gamma(-1,2,1)*Metric(3,4)')

