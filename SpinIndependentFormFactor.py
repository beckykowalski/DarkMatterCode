import numpy as np
import scipy as sp
from scipy import special
import DMConstants
from DMConstants import *

def bessel(x):
    v = (np.sin(x) - x*np.cos(x))/x/x
    return v

def SIFormFactor(Er, A):
    # units of q in fm
    # Helms approximation 
    # recoil energy should be given in keV
    # nucleon mass 
    
    kgtoMeV = 5.6*10.**(29.)
    JtoMeV = 1./MeVtoJ
    keVtoMeV = 0.001

    mNucleon = 931.5 * A # matches number Maria used and gives better FF agreement
    
    q2 = (2. * (mNucleon) * (Er*keVtoMeV))  # momentum transfer in MeV/c
    q = np.sqrt(q2)
    s = 1. # units in fm 
    R = 1.2 * A**OneThird # units in fm
    R1 = np.sqrt((R*R)-(5*s*s))
    hbar_c = 197.2 # MeV*fm
    qr = q*R1/hbar_c # unitless 
    qs = q*s/hbar_c
    aux = R1*q/hbar_c
    j1 = (np.sin(qr) - (qr*np.cos(qr)))/(qr*qr)
#    aux2 = 3*bessel(aux)/aux * np.exp(-(s*s)*q2/hbar_c/hbar_c/2.)
    FF = 0.
    if Er == 0:
        FF = 1.
    else:
        FF = ((3.*j1)/qr) * np.exp(-1.*(s*s*q2/hbar_c/hbar_c/2.)) # factor of 2 missing
#        FF = aux2
    FF2 = FF*FF
    if FF2 > 1.:
        print("ERROR: SPIN INDEPENDENT FORM FACTOR IS GREATER THAN 1. VALUE IS: "+str(FF2)+" AT RECOIL ENERGY "+str(Er)+" J (momentum: "+str(q)+" MeV/c)")
    return FF2

