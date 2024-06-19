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
MElements = [[AtomicMassInkg(132), 132, 1.]]

CalculateEventRate(sig, mchi, MElements, 0, 300, 3000)
