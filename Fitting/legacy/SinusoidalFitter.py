import numpy as np
import ROOT
import matplotlib.pyplot as plt

# for a datafile that has x*365 entries:
#   for entries in file, average every 15 elements
#   returns two lists, average t (corresponding to 15 days), and average rate (average rate within 15 days)
def ReadData(datafile):
    textfile = open(datafile, "r")

    lines = textfile.readlines()

    t = []
    r = []

    for i in lines:
        if i == "time,rate\n":
            continue
        split = i.split(",")
        t.append(float(split[0]))
        r.append(float(split[1]))

        t_avg = []
        r_avg = []
        
        for i in range(len(t)):
        
            if i%15 == 0:
            
                sumt = sum(t[i:i+15])
                sumt /= 15
                sumr = sum(r[i:i+15])
                sumr /= 15

                t_avg.append(sumt)
                r_avg.append(sumr)
    return t, r, r_avg

# sinusoidal fit
# produces image of fit, logfile of fit output, and textfile containing relavent fit info
def Fitter(datafile, frameVal, isnll=False):

    ROOT.gSystem.RedirectOutput(datafile+"_fitter_testPhi.log")
    
    timeData, rateData, rateAvg = ReadData("SimulationForFitting/"+datafile+".txt")

    hist = ROOT.TH1D("hist", "", 120, 0, 1830)

    for i in range(120):
        hist.SetBinContent(i, rateAvg[i])
        err = np.sqrt(rateAvg[i])
        hist.SetBinError(i, err)
        
    nEvents = hist.Integral()

    # only gets histogram for events up to second to last bin, has a zero value for some reason 
    time = ROOT.RooRealVar("time", "time (days)", 0, 1795)
    data = ROOT.RooDataHist("dFhist", "dFhist", ROOT.RooArgList(time), ROOT.RooFit.Import(hist))

    # reasonable seed for S0
    firstbinCont = hist.GetBinContent(0)
    
    SmS = ROOT.RooRealVar("SmS", "SmS", 0, -15, 15)
    SmC = ROOT.RooRealVar("SmC", "SmC", 0, -15, 15)
    S0 = ROOT.RooRealVar("S0", "S0", firstbinCont, firstbinCont - 30, firstbinCont + 30)
#    phi = ROOT.RooRealVar("phi", "phi", 0, 2*np.pi)
#    periodfit = ROOT.RooRealVar("T", "T", 365, 15, 760)
    periodfit = ROOT.RooRealVar("T", "T", 365)

    # fit with maximum likelihood fitter
    pdfsig = ROOT.RooGenericPdf("sinpdf", "sinpdf", "S0 + (SmC*cos(2*pi*time/T) + SmS*sin(2*pi*time/T))", ROOT.RooArgSet(S0, SmC, time, periodfit, SmS))
#    pdfsig = ROOT.RooGenericPdf("sinpdf", "sinpdf", "S0 + (SmC*cos((2*pi/T) * (time + phi)))", ROOT.RooArgSet(S0, SmC, time, periodfit, phi))

    pdfsig.fitTo(data, ROOT.RooFit.SumW2Error(True), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.PrintLevel(1))

    # create negative log likelihood (can plot in different parameter spaces for debugging)
    nll = pdfsig.createNLL(data)

#    frame = time.frame(ROOT.RooFit.Title("toy data"))
    frame = time.frame(ROOT.RooFit.Title("toy data"))

    if isnll==True:
        if frameVal == "S0":
            frame = S0.frame(ROOT.RooFit.Title("toy data"))
        if frameVal == "T":
            frame = periodfit.frame(ROOT.RooFit.Title("toy data"))
        if frameVal == "SmS":
            frame = SmS.frame(ROOT.RooFit.Title("toy data"))
        if frameVal == "SmC":
            frame = SmC.frame(ROOT.RooFit.Title("toy data"))

    if isnll==False:
        data.plotOn(frame)
        pdfsig.plotOn(frame)
    if isnll==True:
        nll.plotOn(frame)
        
    can = ROOT.TCanvas("c", "", 800, 500)
    can.cd()
    frame.Draw()
    if isnll==False:
        can.SaveAs(datafile+".png")
    if isnll==True:
        can.SaveAs(datafile+"_nll_frame"+frameVal+".png")
    
if __name__ == "__main__":
#    Fitter("DC50_Sm10_phi0_T365", "time")
#    Fitter("DC50_Sm5_phi0_T365", "time")
#    Fitter("DC30_Sm10_phi0_T365", "time")
#    Fitter("DC50_Sm10_phi91.25_T365", "time")
    Fitter("DC30_Sm0_phi0_T0", "time")

#    Fitter("DC50_Sm10_phi0_T182.5", "time")


#    Fitter("DC50_Sm10_phi0_T365", "S0", True)
#    Fitter("DC50_Sm10_phi0_T365", "SmS", True)
#    Fitter("DC50_Sm10_phi0_T365", "SmC", True)
#    Fitter("DC50_Sm10_phi0_T365", "T", True)

####### splitting S0 as a flat pdf + Sm pdf didn't work becausde Sm would go to a negative number: can't have a negative pdf
