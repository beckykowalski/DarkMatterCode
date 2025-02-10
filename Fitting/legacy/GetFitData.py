import numpy as np

def GetVar(datafile):

    data = open(datafile, "r")
    lines = data.readlines()

    S0 = 0.
    S0e = 0.
    SmS = 0.
    SmSe = 0.
    SmC = 0.
    SmCe = 0.
    T = 0.
    Te = 0.
    
    for l in lines:
        if l[0:2] == 'S0':
            split = l.split()

            S0 = split[2]
            S0e = split[4]
            print ("S0 in loop = "+str(S0))

        if l[0:3] == 'SmS':
            split = l.split()
           
            SmS = float(split[2])
            SmSe = float(split[4])
           
        if l[0:3] == 'SmC':
            split = l.split()

            SmC = float(split[2])
            SmCe = float(split[4])
           
        if l[0:1] == 'T' and l[0:2] != 'Tr' and l[0:2] != 'Ty':
            split = l.split()

            T = float(split[2])
            Te = float(split[4])



    daytorad = 2.*np.pi / 365
    radtoday = 1/daytorad
    
    print("T = "+str(T)+"+/-"+str(Te))
    print("S0 = "+str(S0)+"+/-"+str(S0e))

    print("SmS = "+str(SmS)+"+/-"+str(SmSe))
    print("SmC = "+str(SmC)+"+/-"+str(SmCe))

    
    Sm = np.sqrt(SmC**2+SmS**2)
    Sme = np.sqrt(SmCe**2+SmSe**2)
    print("Sm = "+str(Sm)+"+/-"+str(Sm))
        
    phi = np.arctan(abs(SmS/SmC))*radtoday
    derivSmC = - (SmS) / (SmC*SmC + SmS*SmS)
    derivSmS = (SmC) / (SmC*SmC + SmS*SmS)
        
    phie = np.sqrt((SmSe*derivSmS)**2 + (SmCe*derivSmC)**2)*radtoday
    print("phi = "+str(phi)+"+/-"+str(phie))

        
#GetVar("DC50_Sm10_phi0_T365_fitter.log")
#GetVar("DC50_Sm5_phi0_T365_fitter.log")
#GetVar("DC50_Sm10_phi91.25_T365_fitter.log")
#GetVar("DC30_Sm10_phi0_T365_fitter.log")
GetVar("DC50_Sm10_phi0_T182.5_fitter.log")
