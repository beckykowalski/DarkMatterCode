import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def sinFunc(DC, Sm, phi, T, x):
    val = 0.
    if T != 0:
        val = DC + Sm*np.cos((2*np.pi/T) * (x + phi))
    if T == 0:
        val = np.random.poisson(DC, 1)
        val = val[0]
    return val

def makeTextFile(DC, Sm, phi, T):

    textfile = open("DC"+str(DC)+"_Sm"+str(Sm)+"_phi"+str(phi)+"_T"+str(T)+".txt","w")
    textfile.write("time,rate\n")
    for i in range(365*5):

        v = sinFunc(DC, Sm, phi, T, i)
        textfile.write(str(i)+","+str(v)+"\n")

    textfile.close()
    return textfile


#makeTextFile(50,10,0,365)
#makeTextFile(50,10,91.25,365)
#makeTextFile(50,5,0,365)
#makeTextFile(50,10,0,182.5)
#makeTextFile(30,10,0,365)
makeTextFile(30,0,0,0)

