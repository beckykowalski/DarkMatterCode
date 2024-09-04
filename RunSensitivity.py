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


#MChis = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

MChis = np.logspace(0.1,3,1000)


MElements = [[AtomicMassInkg(130), 130, .2725*.8], [AtomicMassInkg(128), 128, .2538*.8], [AtomicMassInkg(126), 126, .1506*.8], [AtomicMassInkg(125), 125, .0565*.8], [AtomicMassInkg(124), 124, .0379*.8], [AtomicMassInkg(123), 123, .0071*.8], [AtomicMassInkg(122), 122, .0204*.8], [AtomicMassInkg(121), 121, .0008*.8], [AtomicMassInkg(16), 16, .2]]

# background is 2 counts/keV/kg/day
bkgRate = []
for i in range(20):
    bkgRate.append(2)

sensitivities = []

delta2 = 12.8

# calculating sensitivity for 10-15 keV energy bin, just one bin, 0 keV resolution SI 
for m in MChis:
    sensitivity = SensitivityCalculation(sig, m, 2000*365, delta2, MElements, bkgRate, 10, 20, 20, 400.)
    sensitivities.append(sensitivity)

plt.plot(MChis, sensitivities)
plt.xlabel("$M_{\chi}$ (GeV)")
plt.ylabel("$\sigma_{SI}$ ($\\frac{1}{cm^{2}}$)")
plt.title("Exp = 2TY, bkg = 2 c/kg/keV/d, Energy bin 10-20 keV")
plt.yscale("log")
plt.xscale("log")
plt.savefig("TeO2_Sensitivity_10to20keV_2kgy_bkg2ckgkeVd_2keVReso_20bin.png")
