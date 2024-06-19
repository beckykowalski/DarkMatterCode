import numpy as np
import DMConstants
from DMConstants import *


def ApplyResolution(Er, ResSigma):
    smear = np.random.normal(0, resSigma/2.355) # 2.355 to smear for full width half max, rather than standard deviation
    smearedE = Er + smear
    return smearedE

# all lists, either calculated or derived from data estimates (like exposure/bkg)
def SensitivityCalculation(bkg, delta2, exposure, k0, rho, MChi, MNElements, V0, VEsc, VEAvg, Estart, Estop, ENumEntries):

    # rate, energies, sm, s0, and ff will be arrays of size ENumEntries, ranging from Estart to Estop.
    # these arrays should be summed with bin width Ebin for sensitivity
    
    mSig = 1 # always set initial guess sigma to 1 for scalability
    rate, energies, Sm, s0, ff, xl = GetFullEventRate(k0, rho, M, mSig, MNElements, V0, VEsc, VEAvg, Estart, Estop, ENumEntries)
    Ebin = (Estop - Estart) / ENumEntries

    sig1 = (2*delta2-4)/exposure/Ebin

    Sm2 = 0
    
    for i in range(len(Sm)):
        Sm2 += Sm[i]*Sm[i]/(bkg[i])
        
    sig1 /= Sm2 
    sig1 = np.sqrt(sig1)

    return sig1
        
def NewtonsFunction(sig, Sm, S0, bkg, Ebin, exposure, delta2):
    # to get function, solve eqn 17 in https://arxiv.org/pdf/hep-ph/9912394.pdf, setting it to 0

    sum1 = 0.
    for i in range(len(Sm)):
        sum1 += Sm[i]*Sm[i] / (bkg[i] + sig)
    val = delta2 - 2 - sum1 * (((sig*sig)/2.)*exposure*Ebin) 
    return val

def NewtonsFunctionPrime(sig, Sm, S0, bkg, Ebin, exposure):
    # plug function from "Newtons Function" into derivative calculator solving for sig
    
    sum1 = 0.
    sum2 = 0.
    for i in range(len(Sm)):
        sum1 += Sm[i]*Sm[i] / (bkg[i] + sig*S0[i])
        sum2 += (S0[i]*Sm[i]*Sm[i])/((bkg[i]+S0[i]*sig)*(bkg[i]+S0[i]*sig))
    val = (exposure*Ebin/2.) * (2.*sum1*sig - sum2*sig*sig)
    return val
    
