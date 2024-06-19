import numpy as np
import DMConstants
from DMConstants import *
import ROOT

def GetVEarthAndSun(v0):
    vSun1 = vSolar_R*kmtom
    vSun2 = vSolar_Theta*kmtom + v0
    vSun3 = vSolar_Phi*kmtom

    vEavg_m = vEavg*kmtom
    
    epSpring = EpsilonSpringEquinox
    epSum = EpsilonSummerSolstice

    EarthSpring = vEavg_m * (epSpring[0]*vSun1 + epSpring[1]*vSun2 + epSpring[2]*vSun3)
    EarthSum = vEavg_m * (epSum[0]*vSun1 + epSum[1]*vSun2 + epSum[2]*vSun3)
    SunEarthSqr = vSun1*vSun1 + vSun2*vSun2 + vSun3*vSun3 + vEavg_m*vEavg_m

    delta = np.sqrt(EarthSpring*EarthSpring + EarthSum*EarthSum)

    Ymax = np.sqrt(SunEarthSqr + 2*delta)
    Ymin = np.sqrt(SunEarthSqr - 2*delta)
    Ymean = np.sqrt(SunEarthSqr)
    
    sinMax = EarthSpring/delta
    cosMax = -EarthSum/delta

    phase = np.arcsin(sinMax)

    if cosMax < 0.:
        phase = pi - phase
    phase = phase/twopi + earthT0
    return Ymax, Ymin, Ymean, phase
    
    
def MaxwellianVelocityDistribution(k0, k1, R0, r, E0, v0, vE, vEsc, vmin):
    rateval = 0
    
    kFactor = k1 / k0
    RateFactor = R0 / (E0 * r)
    
    x = vmin/v0
    y = vE/v0
    z = vEsc/v0

    exponential = np.exp( -(z*z) )

    normalization = 1. / (ROOT.TMath.Erf(z) - ((2.*z*exponential)/sqrtpi))
        
    if(x>=0. and x <= z-y):
        rateval = ROOT.TMath.Erf(x+y) - ROOT.TMath.Erf(x-y) - ((4. / sqrtpi) * y * exponential)
    if (x >z-y and x<= z+y):
        rateval =  ROOT.TMath.Erf(z) - ROOT.TMath.Erf(x-y) - ((2. / sqrtpi)*(z + y - x) * exponential)
    if (x > x+y):
        rateval = 0.
        
    return (RateFactor * kFactor * rateval * normalization / 2./ y), rateval, x

