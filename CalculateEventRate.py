import numpy as np
import DMConstants
from DMConstants import *
import CalculateMinimumEventRate
from CalculateMinimumEventRate import *
import IsothermalHaloVelocities
from IsothermalHaloVelocities import *
import SpinIndependentFormFactor
from SpinIndependentFormFactor import *
import ApplyDetectorResolution
from ApplyDetectorResolution import *
from scipy import ndimage as nd

def ReadFormFactorTextFile(FormFactorTxt, Er):
    txt = open(FormFactorTxt, "r")
    lines = txt.readlines()
    FF = 0
    for l in lines:
        energy_data = l.split()
        Energy = float(energy_data[0])
        if Energy == Er:
            FF = float(energy_data[1])
    return FF


# UNITS ARE IN COUNTS/KEV/KG/DAY
# if MNElements is a list > 1, will provide the rate summed over all elements (assuming phase 
# of DM is maximal value to have Sm at its largest)

#def HaloModelEventRate(MChi, sig0, MNElements, Estart, Estop, ENumEntries, resFWHM, FormFactorTxt):
def HaloModelEventRate(MChi, sig0, MNElements, Estart, Estop, ENumEntries, resFWHM):

    # returns in order:
    # 1: list for energy 
    # 2: list for mean event rate
    # 3: list for modulating event rate
    # 4: list for unmodulating event rate
    # 5: list for squared form factor
    # 6: list for minimum velocity corresponding to Er and Mchi
    
    E0 = 0.5 * (MChi*GeVtokg) * (v0*kmtom)*(v0*kmtom)
    k1 = kConstants(k0, vEsc, v0)
    
    Ymax, Ymin, Ymean, phase = GetVEarthAndSun(v0*kmtom)
    energies = []
    eIt = Estart
    for e in range(0, ENumEntries):
        energies.append(eIt)
        incriment = (Estop-Estart)/ENumEntries
        eIt += incriment

    SmTot = np.zeros(ENumEntries)
    S0Tot = np.zeros(ENumEntries)
        
    for Mel in MNElements:

        mu = 1/((1./Mel[0])+(1./MChi))
        Rate = []
        FFarr = []
        erfFunc = []
        xList = []
        sigA = AtomicCross(Mel[1]*0.9315, MChi, Mel[1], sig0)
        R0 = MinEventRate(Mel[1], rho, MChi, sigA, v0*kmtocm)
        r = KineticFactor(MChi, Mel[1]*0.9315)
        t_initial = 0.
        omega = timeomega
        t = 0.
        times = []
        modulatingvel = []

#        for Er in energies:
        for Er in range(len(energies)):
            Er_tbin = 0.
            ErJ = energies[Er]*keVtoJ
            FF = SIFormFactor(energies[Er], Mel[1])
#            FF = 

            vm = MinVelocity(energies[Er], E0, r, v0*kmtom)
            RMax,erf,xMax = MaxwellianVelocityDistribution(k0, k1, R0/sectoday, r, E0/keVtoJ, v0*kmtom, Ymax, vEsc*kmtom, vm)
            RMin,erf,xMin = MaxwellianVelocityDistribution(k0, k1, R0/sectoday, r, E0/keVtoJ, v0*kmtom, Ymin, vEsc*kmtom, vm)

            Sm = (RMax-RMin)/2. 
            S0 = (RMax+RMin)/2. 
            erfFunc.append(erf)

            print(str(Sm)+" for element "+str(Mel[1]))
#            FFarr.append(FF)
#            xList.append(xMax)

            SmTot[Er] += (Sm*FF*Mel[2])
            S0Tot[Er] += (S0*FF*Mel[2])


    if resFWHM == 0:
        return energies, SmTot, S0Tot, FFarr, xList

    else:
        
        # apply gaussian smearing to account for resolution to all rates
        # note: truncate sets range of array for building filter: [(-5*sigma + 0.5) to (5*sigma + 0.5+1)]
        
        newSmRate = nd.gaussian_filter1d(SmTot, resFWHM, truncate=5.0)
        newS0Rate = nd.gaussian_filter1d(S0Tot, resFWHM, truncate=5.0)

        return energies, newSmRate, newS0Rate, FFarr, xList

