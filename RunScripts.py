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

sig = 10.**(-45.)
mchi = 100.
#MElements = [[AtomicMassInkg(132), 132, 1.]]
#MElements = [[AtomicMassInkg(130), 130, .2725*.8], [AtomicMassInkg(128), 128, .2538*.8], [AtomicMassInkg(126), 126, .1506*.8], [AtomicMassInkg(125), 125, .0565*.8], [AtomicMassInkg(124), 124, .0379*.8], [AtomicMassInkg(123), 123, .0071*.8], [AtomicMassInkg(122), 122, .0204*.8], [AtomicMassInkg(121), 121, .0008*.8], [AtomicMassInkg(16), 16, .2]]
#MElements = [[AtomicMassInkg(127.60), 127.6, 1.]]
#MElements = [[AtomicMassInkg(127.60), 127.6, .8], [AtomicMassInkg(16), 16, .2]]
#MElements = [[AtomicMassInkg(16), 16, .2],[AtomicMassInkg(127.60), 127.6, .8]]
MElements = [[AtomicMassInkg(16), 16, .2], [AtomicMassInkg(127.60), 127.6, .8]]
#MElements = [[AtomicMassInkg(6.941), 6.941, .079855], [AtomicMassInkg(16), 16, .368142], [AtomicMassInkg(100), 100, 55.2003]]

#CalculateEventRate(sig, mchi, MElements, 0, 100, 1000, 20)
CalculateEventRate(sig, mchi, MElements, 0, 100, 1000, 0)
