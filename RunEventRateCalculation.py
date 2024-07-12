import matplotlib.pyplot as plt
import DMConstants
from DMConstants import *
import parseDataFunc
from parseDataFunc import *
import CalculateEventRate
from CalculateEventRate import *

# units of sig0 = 1/cm^2
# ElementList syntax: [ [AtomicMassInkg(ElementMassNumber1), ElementMassNumber1, IsotopicFraction1],
#    [AtomicMassInkg(ElementMassNumber2), ElementMassNumber2, IsotopicFraction2], ...]

def CalculateEventRate(sig0, Mchi, ElementList, Estart, Estop, NumberEnergies, resoFWHM):

    ratio = []
    
    outputfile = open("EventRate_"+str(Mchi)+"GeV_sig"+str(sig0)+"_Range"+str(Estart)+"to"+str(Estop)+".csv", "w")
    outputfile.write("Energy,MeanRate,ModulatingRate,UnmodulatingRate,FF,MinV\n")

    print("WIMP MASS = "+str(Mchi)+" GeV")
    print("WIMP CROSS SECTIN = "+str(sig0)+" cm^-2")
    print("Generating event rate from "+str(Estart)+" - "+str(Estop) + " keV")

#    mariaE, mariaR, mariaFF, mariaSM, mariaS0 = PARSEMARIADATA("/home/becky/recopyDMCode/DarkMatterCode/mariaData/rate100GeV_SI1e-9_0to300keV.dat")
    mariaE, mariaR, mariaFF, mariaSM, mariaS0 = PARSEMARIADATA("rate100GeV_SI1e-9_0to300keV_2keVRes.dat")
        
    E, R, Sm, S0, FF, xList = HaloModelEventRate(k0, rho, Mchi, sig0, ElementList, v0, vEsc, vEavg, Estart, Estop, NumberEnergies, resoFWHM)
    
    for i in range(len(R)):
        outputfile.write(str(E[i])+","+str(R[i])+","+str(Sm[i])+","+str(S0[i])+","+str(FF[i])+","+str(xList[i])+"\n")
    
    for i in range(len(R)):
        ratio.append(Sm[i]/mariaSM[i])
    f, (a,b) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios':[4,1]})
    a.plot(E, S0, label="Becky")
    a.plot(mariaE, mariaS0, label="Maria")

    a.legend(loc="best")
    a.set_ylabel("Unmodulating Rate")
    a.set_yscale("log")
    b.plot(E, ratio)
    b.set_xlabel("Energy (keV)")
    b.set_ylabel("Becky/Maria")
    b.set_ylim([0.99,1.1])
    f.tight_layout()
    f.savefig("UnmodulatingRate_June19_100GeV_0to300keV_withReso.png")

    outputfile.close()
