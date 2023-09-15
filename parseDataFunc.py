import numpy as np
import matplotlib.pyplot as plt


def PARSEMARIADATA(datapath):
    datafile = open(datapath, "r")
    energy = []
    s0 = []
    sM = []
    phi = []

    for line in datafile:
        if line == "E:S0:Sm:Phi:FF2SI_1:FF2SD_1\n":
            continue
        else:
            data = line.split()
            energy.append(float(data[0]))
            s0.append(float(data[1]))
            sM.append(float(data[2]))
            phi.append(float(data[3]))
            
            
    rate = []

    for i in range(len(energy)):
        rateval = s0[i] + sM[i]*np.cos(phi[i])
        rate.append(rateval)

    return energy, rate
