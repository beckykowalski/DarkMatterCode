import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

def sinFunc(DC, Sm, phi, T, x):
    return DC + Sm*np.cos(2*np.pi*x/T + phi)

def makeTextFile(DC, Sm, phi, T):

    textfile = open("DC"+str(DC)+"_Sm"+str(Sm)+"_phi"+str(phi)+"_T"+str(T)+".txt","w")
    textfile.write("time,rate\n")
    for i in range(365*5):

        v = sinFunc(DC, Sm, phi, T, i)
        textfile.write(str(i)+","+str(v)+"\n")

    textfile.close()
    return textfile


makeTextFile(50,10,0,365)
makeTextFile(50,10,91.25,365)
makeTextFile(50,5,0,365)
makeTextFile(50,10,0,182.5)
makeTextFile(30,10,0,365)

