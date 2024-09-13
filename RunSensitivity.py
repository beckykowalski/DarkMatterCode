import CalculateEventRate
from CalculateEventRate import *
import CalculateMinimumEventRate
from CalculateMinimumEventRate import *
import DMConstants
import IsothermalHaloVelocities
from IsothermalHaloVelocities import *
import SpinIndependentFormFactor
from SpinIndependentFormFactor import *
import VectorAddition
from VectorAddition import *
import RunEventRateCalculation
from RunEventRateCalculation import *
import CalculateSensitivity
from CalculateSensitivity import *
import matplotlib.pyplot as plt
import numpy as np

sig = 1.

MChis = np.logspace(0.1,3,1000)
#MElements = [[AtomicMassInkg(130), 130, .2725*.8], [AtomicMassInkg(128), 128, .2538*.8], [AtomicMassInkg(126), 126, .1506*.8], [AtomicMassInkg(125), 125, .0565*.8], [AtomicMassInkg(124), 124, .0379*.8], [AtomicMassInkg(123), 123, .0071*.8], [AtomicMassInkg(122), 122, .0204*.8], [AtomicMassInkg(121), 121, .0008*.8], [AtomicMassInkg(16), 16, .2]]
MElements = [[AtomicMassInkg(132), 132, 1.]]

# background is 2 counts/keV/kg/day
bkgRate = 500
#bkgRate = 1

sensitivities = []

delta2 = 12.8

txt = open("Sensitivity_Xe_0keVres_2Y_12.8_2bkg", "w")
txt.write("m,sig\n")

# calculating sensitivity for 10-15 keV energy bin, just one bin, 0 keV resolution SI 
for m in MChis:
#    sensitivity = SensitivityCalculation(sig, m, 2000*365, delta2, MElements, bkgRate, 10, 15, 1, 40.)
    sensitivity = SensitivityCalculation(sig, m, 2000*365, delta2, MElements, bkgRate, 10, 15, 1, 0.)
    sensitivity *= 1e36
    sensitivities.append(sensitivity)
    txt.write(str(m)+","+str(sensitivity)+"\n")

txt.close()
    
plt.plot(MChis, sensitivities)
plt.xlabel("$M_{\chi}$ (GeV)")
plt.ylabel("$\sigma_{SI}$ ($cm^{2}$)")
#plt.ylabel("$\sigma_{SI}$ ($\\frac{1}{cm^{2}}$)")
plt.title("Exp = 2TY, bkg = 2 c/kg/keV/d, Energy bin 10-15 keV, Reso = 2 keV")
plt.yscale("log")
plt.xscale("log")
#plt.ylim(.000001, .01)
plt.savefig("Xe_Sensitivity_10-15keV_2kgy_bkg500ckgkeVd_2keVRes.png")
