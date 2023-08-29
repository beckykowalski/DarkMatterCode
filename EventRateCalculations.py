import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import special
import DMConstants
from DMConstants import *
import MaxwellianVelocityDist
from MaxwellianVelocityDist import *
import pandas as pd
from scipy.special import logsumexp

def MaxwellianVelocityDistribution(k0, k1, R0, r, E0, v0, vE, vEsc, vmin):
    rateval = 0
    
    kFactor = k1 / k0
    RateFactor = R0 / (E0 * r)
    
    x = vmin/v0
    y = vE/v0
    z = vEsc/v0
    
#    if mod == False:
#        y = vE/v0 # when time dependence is introduced, this becomes : (vsun+vearth(t))/v0
#    if mod == True:
#        y = vmod

    exponential = np.exp( -(z*z) )
        
    normalization = 1. / (sp.special.erf(z) - 2./np.sqrt(np.pi)*z*exponential)
        
    if(x>=0. and x <= z-y):
        rateval = sp.special.erf(x+y) - sp.special.erf(x-y) - (4. / np.sqrt(np.pi) * y * exponential)
    if (x >z-y and x<= z+y):
        rateval =  sp.special.erf(z) - sp.special.erf(x-y) - (2. / np.sqrt(np.pi)*(z + y - x) * exponential)
    if (x > x+y):
        rateval = 0.
        
        
    return (RateFactor * kFactor * rateval * normalization / 2. / y)
                                                                                                                    
def MaxwellianVelocityDistributionTimeDependence(k0, k1, R0, r, E0, v0, vEnom, vEsc, vmin, t, t0, omega, vEavg):
    rateval = 0
    
    kFactor = k1 / k0
    RateFactor = R0 / (E0 * r)

    vt = GetTimeDependentEarthVelocity(vEavg, t, t0, omega)

    x = vmin/v0
    # add time dependent velocity to y from rotation of earth (averaged at 29.8 km/s +/- sinusoidal component derived from sun's peculiar velocity at differnet 
    ## times of the year
    y = (vE + vt)/v0 
    z = vEsc/v0 
   
    exponential = np.exp( -(z*z) )
        
    normalization = 1. / (sp.special.erf(z) - 2./np.sqrt(np.pi)*z*exponential)
        
    if(x>=0. and x <= z-y):
        rateval = sp.special.erf(x+y) - sp.special.erf(x-y) - (4. / np.sqrt(np.pi) * y * exponential)
    if (x >z-y and x<= z+y):
        rateval =  sp.special.erf(z) - sp.special.erf(x-y) - (2. / np.sqrt(np.pi)*(z + y - x) * exponential)
    if (x > x+y):
        rateval = 0.
        
        
    return (RateFactor * kFactor * rateval * normalization / 2. / y), vt
                                                                                                                    

def GetTimeDependentEarthVelocity(vE, t, t0, omega):

    ### taken as equation 27 in https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.85.1561 (Katherine Freese Colloqium on Annual Modulation)

    
    w1 = (EpsilonSpringEquinox[0]*vSolar_R + EpsilonSpringEquinox[1]*vSolar_Theta + EpsilonSpringEquinox[2]*vSolar_Phi)
    w2 = (EpsilonSummerSolstice[0]*vSolar_R + EpsilonSummerSolstice[1]*vSolar_Theta + EpsilonSummerSolstice[2]*vSolar_Phi)
    b = np.sqrt(w1*w1 + w2*w2) 
    ####    vOb = vsun + b*vE*np.cos(omega * (t - t0))
    # this commponent will add with vy in Maxwellian Velocity distribution 
    vOb = b*vE*np.cos(omega * (t - t0))
    return(vOb)
'''
    # arbitrary DM Constant date set. Begin considering period begining March 22st 2018
    # this means t should incriment in units of DAYS
    omega = 0.0172 # units is 1/days

    # Set velocity due to revolution of Earth around Sun as function of time
    vRotationOfEarthAroundSun = []
    vRotationOfEarthAroundSun[0] = (0.9941 * np.cos(omega * t)) - (0.0504 * np.sin(omega * t))
    vRotationOfEarthAroundSun[1] = (0.1088 * np.cos(omega * t)) + (0.4946 * np.sin(omega * t))
    vRotationOfEarthAroundSun[2] = (0.0042 * np.cos(omega * t)) - (0.8677 * np.sin(omega * t))

    # setting peculiar velocity of sun relative to position in Milky Way
    SunPeculiar = [vSolar_R, vSolar_Theta, vSolar_Phi]

    # finding relative velocity (as function of time) for sun to Earth, with Sun's galactic position into account. 
    vFluxToSun = []
    vFluxToSun[0] = SunPeculiar[0] - vRotationOfEarthAroundSun[0]
    vFluxToSun[1] = SunPeculiar[1] + v0 - vRotationOfEarthAroundSun[1]
    vFluxToSun[2] = SunPeculiar[2] - vRotationOfEarthAroundSun[2] 

    # summar solstice: expected max of velocity flux. Spring equinox: expected midpoint of velocity flux. Getting average
    ## Earth velocity for modulating/static DM event rate normalzed to a max and mean value of velocity
    w1 = vEavg * (EpsilonSpringEquinox[0]*vFluxToSun[0] + EpsilonSpringEquinox[1]*vFluxToSun[1] + EpsilonSpringEquinox[2]*vFluxToSun[2])
    w2 = vEavg * (EpsilonSummerSolstice[0]*vFluxToSun[0] + EpsilonSummerSolstice[1]*vFluxToSun[1] + EpsilonSummerSolstice[2]*vFluxToSun[2])

    wT2 = vFluxToSun[0]*vFluxToSun[0] + vFluxToSun[1]*vFluxToSun[1] + vFluxToSun[2]*vFluxToSun[2]

    ExpectedModulationAmp = np.sqrt(w1*w1 + w2*w2)
    
    vEMax = np.sqrt(wT2 + 2*ExpectedModulationAmp) / vDis

    
    Max = MaxwellianVelocityDistribution(k0, k1, R0, r, E0, v0, vEMax, vEsc, vmin)
    Min = MaxwellianVelocityDistribution(k0, k1, R0, r, E0, v0, vEMin, vEsc, vmin)

    Sm = (Max-Min) / 2. / v0
'''

def MinEventRate(A, rho, Mchi, sigA, v0):

    # units should be: masses in GeV, distance in cm, velocity in cm/s

    txt = open("ValuesInEventRate.txt", "w")
    Na = 6.02 * 10.**26. # 10^26 kg*mol^-1 or 10^23 g*mol^-1

    NumberCount = Na / A
    MassNormalizedDensity = rho / Mchi
    numericalConst = 2. / np.sqrt(np.pi)

    Rate = numericalConst * NumberCount * MassNormalizedDensity * sigA * v0

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
    qfm = np.sqrt(2.*MN*Er) * (5.344*10.**-19.) / .1973
    if Er == 1.602e-14:
        print("For recoil energy 100 keV: qr = "+str(qr)+", and FF is "+str(FF2))
    return FF2,qfm


k0 = (np.pi * v0*v0)**(3./2.)

rho = 0.3 # GeV / (c^2 cm^3), DM local density


sig0 = 10.**(-45.) # cm^2, WIMP-Nucleon cross section 
MChi = 100. # GeV/c^2, mass of WIMP

# list of list: each element is a list of [mass_in_kg, atomic_number]
MElements = [[AtomicMassInkg(132), 132]]
#MElements = [[AtomicMassInkg(40), 40]]

E0 = 0.5 * (MChi*GeVtokg) * (v0 * kmtom)*(v0 * kmtom)

k1 = kConstants(k0, vEsc, v0) 

# units in keV
#energies = np.linspace(.0001, 110, 10000)
energies = np.linspace(.0001, 100, 80)

mariadata = pd.read_csv("MariaXenon.csv")
mariaE = mariadata["Xval"].values.tolist()
mariaR = mariadata["Yval"].values.tolist()


for Mel in MElements:

    # create array to plot the event rate (y axis)
    Rate = []

    # initialize not time/energy dependent variables for event rate 
    sigA = AtomicCross(Mel[0], MChi*GeVtokg, Mel[1], sig0)
    R0 = MinEventRate(Mel[1], rho, MChi, sigA, v0*kmtocm) 
    r = KineticFactor(MChi * GeVtokg, Mel[0])
    t_initial = 0.
#    omega = 1./365.
    omega = timeomega

    times = []
    modulatingvel = []
    # sum event rate for every day in a year  
    for Er in energies:
        Er_tbin = 0.
        ErJ = Er*keVtoJ
        FF,q = SIFormFactor(ErJ, Mel[0], Mel[1])
        vm = MinVelocity(ErJ, E0, r, v0*kmtom)
        for t in range(366): 
#        times.append(float(t))

            R,vt = MaxwellianVelocityDistributionTimeDependence(k0, k1, R0/sectoday, r, E0/keVtoJ, v0*kmtom, vE*kmtom, vEsc*kmtom, vm, t, t_initial, omega, vEavg)

#            if Er == energies[-1]:
#                modulatingvel.append(vt)
            #        R = MaxwellianVelocityDistribution(k0, k1, R0/sectoyear, r, E0/keVtoJ, v0*kmtom, vE*kmtom, vEsc*kmtom, vm)
            #        print(R)
#            Rate.append(R*FF)
            Er_tbin += R*FF
        Rate.append(Er_tbin)

            
    plt.plot(energies, Rate)
#    plt.plot(times, modulatingvel)
    plt.yscale("log")
    plt.xlabel("Recoil Energy (keV)")
#    plt.xlabel("Time (days)")
    plt.ylabel("dR/dE (cp/yr/keV/kg)")
#    plt.ylabel("Velocity (km/s)")
    plt.ylim([.0000001, .1])
    plt.grid(True)
    plt.savefig("EventRate_vs_energy_xenon_timedependent.png")
        

