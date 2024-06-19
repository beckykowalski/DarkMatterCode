import numpy as np
import DMConstants
from DMConstants import *
from scipy import special
import VectorAddition
from VectorAddition import *

def MinEventRate(A, rho, Mchi, sigA, v0):

    # units should be: masses in GeV, distance in cm, velocity in cm/s

    txt = open("ValuesInEventRate.txt", "w")
    Na = 6.022 * 10.**26. # 10^26 kg*mol^-1 or 10^23 g*mol^-1

    NumberCount = Na / A
    MassNormalizedDensity = rho / Mchi
    numericalConst = 2. / sqrtpi

    Rate = numericalConst * NumberCount * MassNormalizedDensity * sigA * v0

    return Rate

def kConstants(k0, vesc, v0):

    erfFunc = special.erf(vesc / v0)
    numericalConst = 2. / sqrtpi
    velFunc = vesc / v0
    exponential = np.exp( - (vesc*vesc)/(v0*v0) )

    k1 = k0 * (erfFunc - (numericalConst * velFunc * exponential))

    return k1


def AtomicCross(MA, Mchi, A, sig0):

    # units should be in kg and cm^2
    MNucleon = .9315

    muNucleon = (MNucleon * Mchi) / (MNucleon + Mchi)
    muAtom = (MA * Mchi) / (MA + Mchi)

    sigA = sig0 * (muAtom*muAtom) / (muNucleon*muNucleon) * A*A

    return sigA
    
def MinVelocity(Er, E0, r, v0):

    vmin = v0 * np.sqrt(Er*keVtoJ / (E0*r))
    return vmin

def KineticFactor(Mchi, MA):
    #  Mchi and MA in kg

    numerator = 4 * MA * Mchi
    denominator = (Mchi + MA) * (Mchi + MA)

    kinFact = numerator / denominator
    return kinFact
    
def AtomicMassInkg(NumberOfNucleons):
#    mN = 1.67 * 10.**(-27)
    mN = 0.9315 * GeVtokg
    return mN * NumberOfNucleons

