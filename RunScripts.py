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


sig = 10.**(-45.)
mchi = 100.
MElements = [[AtomicMassInkg(132), 132, 1.]]
#MElements = [[AtomicMassInkg(130), 130, .2725*.8], [AtomicMassInkg(128), 128, .2538*.8], [AtomicMassInkg(126), 126, .1506*.8], [AtomicMassInkg(125), 125, .0565*.8], [AtomicMassInkg(124), 124, .0379*.8], [AtomicMassInkg(123), 123, .0071*.8], [AtomicMassInkg(122), 122, .0204*.8], [AtomicMassInkg(121), 121, .0008*.8], [AtomicMassInkg(16), 16, .2]]
#MElements = [[AtomicMassInkg(130), 130, 1.]]
#MElements = [[AtomicMassInkg(130), 130, .2725], [AtomicMassInkg(128), 128, .2538], [AtomicMassInkg(126), 126, .1506], [AtomicMassInkg(125), 125, .0565], [AtomicMassInkg(124), 124, .0379], [AtomicMassInkg(123), 123, .0071], [AtomicMassInkg(122), 122, .0204], [AtomicMassInkg(121), 121, .0008]]
#MElements = [[AtomicMassInkg(127.60), 127.6, 1.]]
#MElements = [[AtomicMassInkg(127.60), 127.6, .8], [AtomicMassInkg(16), 16, .2]]
#MElements = [[AtomicMassInkg(6), 6, .0758*.079855], [AtomicMassInkg(7), 7, 0.9241*0.079855], [AtomicMassInkg(16), 16, .368142], [AtomicMassInkg(100), 100, .552003]]
#MElements = [ [AtomicMassInkg(6.941), 6.941, 0.079855], [AtomicMassInkg(16), 16, .368142], [AtomicMassInkg(100), 100, .552003]]

#CalculateEventRatePlotter(sig, mchi, MElements, 1.4, 9.2, 1000, 200)
CalculateEventRatePlotter(sig, mchi, MElements, 0, 300, 3000, 0)
#CalculateEventRatePlotter(sig, mchi, MElements, 10, 15, 1, 20)
#CalculateEventRate(sig, mchi, MElements, 0, 100, 1000, 0)

# calculating sensitivity for 10-15 keV energy bin, just one bin, 0 keV resolution SI 
