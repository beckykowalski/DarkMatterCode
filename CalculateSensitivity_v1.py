import numpy as np
import DMConstants
from DMConstants import *
import CalculateEventRate
from CalculateEventRate import *


# all lists, either calculated or derived from data estimates (like exposure/bkg)
def SensitivityCalculation(sig0, MChi, exposure, delta2, MNElements, bkg, Estart, Estop, ENumEntries, resFWHM):

#    txt = open("Xe_rateForSensitivity_res"+str(resFWHM)+"_bkg"+str(bkg)+".txt", "w")

    # rate, energies, sm, s0, and ff will be arrays of size ENumEntries, ranging from Estart to Estop.
    # these arrays should be summed with bin width Ebin for sensitivity

    energies = []
    Sm = []
    S0 = []
    ff = []
    xl = []

    if ENumEntries == 1:

        numEntriesSens = 100
        energies, Sm, S0, ff, xl = HaloModelEventRate(MChi, sig0, MNElements, Estart, Estop, numEntriesSens, resFWHM)

    Ebin = (Estop - Estart) / ENumEntries
    Sm2 = 0
    Sm2 = Sm2*Sm2
    for i in range(len(Sm)):
        Sm2 += Sm[i]*Sm[i]/bkg[i]
    Sm2 /= bkg

    sig1 = np.sqrt((2*delta2-4)/Sm2/exposure/Ebin)

    # conversion from cm^2 to pb
#    sig1 *= 1e36

    return sig1
        
