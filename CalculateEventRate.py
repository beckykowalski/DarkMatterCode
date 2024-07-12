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


def HaloModelEventRate(k0, rho, MChi, sig0, MNElements, V0, Vesc, VEavg, Estart, Estop, ENumEntries, resFWHM):

    # returns in order:
    # 1: list for energy 
    # 2: list for mean event rate
    # 3: list for modulating event rate
    # 4: list for unmodulating event rate
    # 5: list for squared form factor
    # 6: list for minimum velocity corresponding to Er and Mchi
    
    E0 = 0.5 * (MChi*GeVtokg) * (V0*kmtom)*(V0*kmtom)
    k1 = kConstants(k0, Vesc, V0)
    
    Ymax, Ymin, Ymean, phase = GetVEarthAndSun(V0*kmtom)
    energies = []
    eIt = Estart
    for e in range(0, ENumEntries):
        energies.append(eIt)
        incriment = (Estop-Estart)/ENumEntries
        eIt += incriment
        
    for Mel in MNElements:

        print("Element mass number "+str(Mel[1]))
        
        mu = 1/((1./Mel[0])+(1./MChi))
        Rate = []
        FFarr = []
        erfFunc = []
        xList = []
        sigA = AtomicCross(Mel[1]*0.9315, MChi, Mel[1], sig0)
        R0 = MinEventRate(Mel[1], rho, MChi, sigA, V0*kmtocm)
        r = KineticFactor(MChi, Mel[1]*0.9315)
        t_initial = 0.
        omega = timeomega
        t = 0.
        times = []
        modulatingvel = []

        SmList = []
        S0List = []
        for Er in energies:
            Er_tbin = 0.
            ErJ = Er*keVtoJ
            FF = SIFormFactor(Er, Mel[1])

            vm = MinVelocity(Er, E0, r, V0*kmtom)
            RMax,erf,xMax = MaxwellianVelocityDistribution(k0, k1, R0/sectoday, r, E0/keVtoJ, V0*kmtom, Ymax, Vesc*kmtom, vm)
            RMin,erf,xMin = MaxwellianVelocityDistribution(k0, k1, R0/sectoday, r, E0/keVtoJ, V0*kmtom, Ymin, Vesc*kmtom, vm)
            RMean,erf,xMean = MaxwellianVelocityDistribution(k0, k1, R0/sectoday, r, E0/keVtoJ, V0*kmtom, Ymean, Vesc*kmtom, vm)

            Sm = (RMax-RMin)/2. * Mel[2] 
            S0 = (RMax+RMin)/2. * Mel[2]
            R = RMean * Mel[2]
            
            erfFunc.append(erf)
            FFarr.append(FF)
            xList.append(xMax)

            Rate.append(R*FF)
            SmList.append(Sm*FF)
            S0List.append(S0*FF)


    if resFWHM == 0:
        return energies, Rate, SmList, S0List, FFarr, xList

    else:
        
        # apply gaussian smearing to account for resolution to all rates
        # note: truncate sets range of array for building filter: [(-5*sigma + 0.5) to (5*sigma + 0.5+1)]
        
        newSmRate = nd.gaussian_filter1d(SmList, resFWHM, truncate=5.0)
        newS0Rate = nd.gaussian_filter1d(S0List, resFWHM, truncate=5.0)
        newRate = nd.gaussian_filter1d(Rate, resFWHM, truncate=5.0)
        
        return energies, newRate, newSmRate, newS0Rate, FFarr, xList

