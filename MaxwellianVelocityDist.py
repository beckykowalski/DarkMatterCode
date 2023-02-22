import DMConstants
import EventRateCalculations
import numpy as np

# Analytical calculation given maxwellian velocity distribuiton
# f(v, sigma) = 1/[sqrt(2pi)*sig] * exp[-|v|^2/(2*sigma^2)]
### sigma = velocity dispersion
### analytical form of vmin to vescape velocity integral f(v, sigma) / v

def MaxwellianVelocityDistribution(k0, k1, R0, r, E0, v0, vE, vEsc, vmin):
    rateval = 0
    
    kFactor = k1 / k0
    RateFactor = R0 / (E0 * r)
    
    x = vmin/v0
    y = vE/v0 # when time dependence is introduced, this becomes : (vsun+vearth(t))/v0
    z = vEsc/v0
    
    exponential = np.exp( -(z*z) )
    
    normalization = 1. / (sp.special.erf(z) - 2./np.sqrt(np.pi)*z*exponential)
    
    if(x>=0. and x <= z-y):
        rateval = sp.special.erf(x+y) - sp.special.erf(x-y) - (4. / np.sqrt(np.pi) * y * exponential)
    if (x >z-y and x<= z+y):
        rateval =  sp.special.erf(z) - sp.special.erf(x-y) - (2. / np.sqrt(np.pi)*(z + y - x) * exponential)
    if (x > x+y):
        rateval = 0.
                
    
    return (RateFactor * kFactor * rateval * normalization / 2. / y)

def VelocitySum(v0, vE):

    sunV = []
    sunV[0] = sunPeculiar_0
    sunV[1] = sunPeculiar_1 + v0
    sunV[2] = sunPeculiar_2

    E1dotSun = vEavg * ( vSun[0]*epsilon1_R + vSun[1]*epsilon1_Theta + vSun[2]*epsilon1_Phi ) 
    E2dotSun = vEavg * ( vSun[0]*epsilon2_R + vSun[1]*epsilon2_Theta + vSun[2]*epsilon2_Phi ) 
    
    vEvSunSq = sunV[0]*sunV[0] + sunV[1]*sunV[1] + sunV[2]*sunV[2] + vEavg*vEavg

    deltaV = np.sqrt(E1dotSun*E1dotSun + E2dotSun*E2dotSun)

    MaxModVel = np.sqrt(vEvSunSq + 2*deltaV) / v0
    MinModVel = np.sqrt(vEvSunSq - 2*deltaV) / v0
    AvgModVel = np.sqrt(vEvSunsq) / v0

    MaxSpringEquinox = E1dotSun / delta # should be purely directional
    MaxSummerSolstice = -E2dotSun / delta # should be purely directional

    phi = np.arcsin(MaxSpringEquinox)
    if cosMax < 0:
        phi = np.pi = phi
    phi = phi/(2.*np.pi) + earthT0 # units in years

    return MaxModVel, MinModVel, AvgModVel, phi
