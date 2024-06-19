import numpy as np

def AddVec(amp1, phi1, amp2, phi2):
# updates amp1 and phi1 
    if amp1 == 0.:
        if amp2 != 0:
            print("in addVec first loop")
            amp1 = amp2
            phi1 = phi2
        return amp1, phi1
    
    if phi1 == phi2:
        print("in second loop")
        amp1 += amp2
        return amp1, phi1

    if amp2 == 0:
        print("in third loop")
        return amp1, phi1

    
    CosVec = amp1*np.cos(phi1) + amp2*np.cos(phi2)
    SinVec = amp1*np.sin(phi1) + amp2*np.sin(phi2)
    
    amp1 = np.sqrt(CosVec*CosVec + SinVec*SinVec)
    phi1 = np.arcsin(CosVec/amp1)

    print("modifying amp1 with phi info")
    if SinVec < 0.:
        phi1 = pi - phi1
    if phi1 > 0:
        phi1 = phi1 + twopi
    return amp1, phi1

