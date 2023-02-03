import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import special
import DMConstants
from DMConstants import *
import pandas as pd

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
        


def MinEventRate(A, rho, Mchi, sigA, v0):

    # units should be: masses in GeV, distance in cm, velocity in cm/s

    txt = open("ValuesInEventRate.txt", "w")
    Na = 6.02 * 10.**26. # 10^26 kg*mol^-1 or 10^23 g*mol^-1

    NumberCount = Na / A
    MassNormalizedDensity = rho / Mchi
    numericalConst = 2. / np.sqrt(np.pi)

    Rate = numericalConst * NumberCount * MassNormalizedDensity * sigA * v0

    txt.write("Avagadros number is "+str(Na)+"\n")
    txt.write("atomic number of xe is "+str(A)+"\n")
    txt.write("dm density is "+str(rho)+"\n")
    txt.write("dm mass is "+str(Mchi)+"\n")
    txt.write("atomic cross section is "+str(sigA)+"\n")
    txt.write("min velocity is "+str(v0)+"\n")
    txt.close()
    return Rate

def kConstants(k0, vesc, v0):

    erfFunc = special.erf(vesc / v0)
    numericalConst = 2. / np.sqrt(np.pi)
    velFunc = vesc / v0
    exponential = np.exp( - (vesc*vesc)/(v0*v0) )

    k1 = k0 * (erfFunc - numericalConst * velFunc * exponential)

    return k1


def AtomicCross(MA, Mchi, A, sig0):

    # units should be in kg and cm^2

    MNucleon = 1.67 * 10.**(-27)

    muNucleon = (MNucleon * Mchi) / (MNucleon + Mchi)
    muAtom = (MA * Mchi) / (MA + Mchi)

    sigA = sig0 * (muAtom*muAtom) / (muNucleon*muNucleon) * A*A

    return sigA
    
#### confirmed all three calculations equal the same thing when all units are in joules     
def MinVelocity(Er, E0, r, v0):

    # v0 should be cm/s
    # E0 and Er should be in same units (say Joules)
    
    vmin = v0 * np.sqrt(Er / (E0*r))
    return vmin

def KineticFactor(Mchi, MA):

    #  Mchi and MA in kg

    numerator = 4 * MA * Mchi
    denominator = (Mchi + MA) * (Mchi + MA)

    kinFact = numerator / denominator
    return kinFact
    
def AtomicMassInkg(NumberOfNucleons):
    mN = 1.67 * 10.**(-27)
    return mN * NumberOfNucleons


def SIFormFactor(Er, MN, A):
    # units of q in fm
    # Helms approximation 

    kgtoMeV = 1./(1.79*10.**-30.)
    JtoMeV = 1./(1.609*10.**-13.)

    degtorad = np.pi/180.
    
    q = np.sqrt(2. * (MN*kgtoMeV) * (Er*JtoMeV))  # momentum transfer in MeV/c
#    q = np.sqrt(2. * (MN*kgtoMeV) * (Er/keVtoJ))  # momentum transfer in MeV/c
    s =  1. # units in fm 
    R = 1.2 * A **(1./3.) # units in fm
    R1 = np.sqrt(R*R-5*s*s)
    hbar_c = 197.3 # MeV*fm
    qr = q*R1/hbar_c # unitless 
    qs = q*s/hbar_c
    j1 = (np.sin(qr))/(qr*qr) - (np.cos(qr))/(qr)
    FF = 0.
    if Er == 0:
        FF = 1.
    else:
        FF = ((3.*j1)/qr) * np.exp(-1. * ((qs*qs)/2.))
    FF2 = FF*FF
    if FF2 > 1.:
        print("ERROR: SPIN INDEPENDENT FORM FACTOR IS GREATER THAN 1. VALUE IS: "+str(FF2)+" AT RECOIL ENERGY "+str(Er)+" J (momentum: "+str(q)+" MeV/c)")
    q2 = q * (3.*10.**23.)
    if Er == 1.602e-14:
        print("For recoil energy 100 keV: qr = "+str(qr)+", and FF is "+str(FF2))
    return FF2,qr
    
k0 = (np.pi * v0*v0)**(3./2.)

rho = 0.3 # GeV / (c^2 cm^3), DM local density


sig0 = 10.**(-45.) # cm^2, WIMP-Nucleon cross section 
MChi = 100. # GeV/c^2, mass of WIMP

# list of list: each element is a list of [mass_in_kg, atomic_number]
#MElements = [[AtomicMassInkg(132), 132]]
MElements = [[AtomicMassInkg(40), 40]]

E0 = 0.5 * (MChi*GeVtokg) * (v0 * kmtom)*(v0 * kmtom)

k1 = kConstants(k0, vEsc, v0) 

# units in keV
#energies = np.linspace(.0001, 110, 10000)
energies = np.linspace(.0001, 100, 100)

mariadata = pd.read_csv("/mnt/c/Users/bkow1/Downloads/ArgonEvR.csv")


mariaE = mariadata["Energy"].values.tolist()
mariaF = mariadata["DRDE"].values.tolist()

for Mel in MElements:

    sigA = AtomicCross(Mel[0], MChi*GeVtokg, Mel[1], sig0)
    R0 = MinEventRate(Mel[1], rho, MChi, sigA, v0*kmtocm) 
    
    Rate = []
    vels = []

    formfactorarray = []
    momentumarray = []
    
    r = KineticFactor(MChi * GeVtokg, Mel[0])
    ecounter = 0
    
    for Er in energies:
        
        ErJ = Er*keVtoJ
        if Er == 100.:
            print(ErJ)
        vm = MinVelocity(ErJ, E0, r, v0*kmtom)
        vels.append(vm)
        
        R = MaxwellianVelocityDistribution(k0, k1, R0/sectoyear, E0/keVtoJ, r, v0*kmtom, vE*kmtom, vEsc*kmtom, vm)
        
        FF,q = SIFormFactor(ErJ, Mel[0], Mel[1])
        Rate.append(R*FF)
        formfactorarray.append(FF)
        momentumarray.append(q)
    plt.plot(energies, Rate, label="Becky")
#    plt.plot(energies, formfactorarray, label="Becky")
    plt.plot(mariaE, mariaF, label="Maria")
#    plt.plot(energies, momentumarray)
#    plt.plot(momentumarray, formfactorarray)
    plt.yscale("log")
    plt.xlabel("Recoil Energy (keV)")
#    plt.xlabel("q")
    plt.ylabel("dR/dE (cpd/keV/kg)")
#    plt.ylabel("$F^{2}$")
#    plt.ylabel("Momentum (MeV/c)")
#    plt.xlabel("qr (unitless)")
 #   plt.ylim([.000000000000001, 1.])
    plt.grid(True)
    plt.ylim([.0000001, 1.])
    plt.legend(loc='best')
#    plt.savefig("XenonFFvsQR.png")
    plt.savefig("Argon_dRdEvE.png")
        
