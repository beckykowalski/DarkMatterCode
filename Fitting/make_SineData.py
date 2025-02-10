import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def sinFunc(DC, Sm, phi, T, x):
    val = 0.
    if T != 0:
        valt = DC + Sm*np.cos((2*np.pi/T) * (x + phi))
        val = np.random.poisson(valt,1)
        val = val[0]
 
    if T == 0:
        val = np.random.poisson(DC, 1)
        val = val[0]
    return val

def makeTextFile(DC, Sm, phi, T, it):

#    textfile = open("Sm13_phase_data/Sm13_phaseData_"+str(it)+".txt","w")
    textfile = open("flat_data/FlatData_"+str(it)+".txt","w")
#    textfile = open("DC"+str(DC)+"_Sm"+str(Sm)+"_phi"+str(phi)+"_T"+str(T)+"_flat.txt","w")
    textfile.write("time,rate\n")

    data = []
    time = []
    for i in range(1830):

        v = sinFunc(DC, Sm, phi, T, i)
        textfile.write(str(i)+","+str(v)+"\n")
        time.append(i)
        data.append(v)
        
    textfile.close()
    
#    plt.errorbar(time, data, fmt=".")
#    plt.xlabel("time (days)")
#    plt.ylabel("counts")
#    plt.savefig("DC_50_Sm_10_phase_0_T_365_testData.png")
#    plt.savefig("Sm13_phase_data/Sm13_phaseData_file"+str(it)+".png")
#    plt.cla()
#    plt.clf()
#    print("final time bin = "+str(time[-1]))
    return textfile


#makeTextFile(50,10,0,365)
#makeTextFile(50,10,91.25,365)
#makeTextFile(50,5,0,365)
#makeTextFile(50,10,0,182.5)
#makeTextFile(30,10,0,365)

for i in range(1000):
    print("it "+str(i))
    makeTextFile(30,0,0,0,i)

