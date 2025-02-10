import numpy as np
import ROOT
import matplotlib.pyplot as plt


def FitterESlice(EStart, EStop):

    hist = ROOT.TH1D("timedata", "", 60, 0, 60)
    
    for t in range(60):
        f = ROOT.TFile("rootfile_data/SimulationHist_"+str(t)+".root")
        h1 = f.Get("toydata")

        counter = 0
        
        for e in h1:
            if e >= EStart and e < EStop:
                counter += 1

        if counter == 0:
            hist.SetBinContent(t, .0001)
        else:
            hist.SetBinContent(t, counter)

    nEvents = hist.Integral()
    
    time = ROOT.RooRealVar("time", "time", 0, 59)
    data = ROOT.RooDataHist("dFhist", "dFhist", ROOT.RooArgList(time), ROOT.RooFit.Import(hist))

    Sm = ROOT.RooRealVar("Sm", "Sm", 5, 0., 10)
    S0 = ROOT.RooRealVar("S0", "S0", 15, 5, 25)
    periodfit = ROOT.RooRealVar("T", "T", 12, .01, 25)
    phasefit = ROOT.RooRealVar("phi", "phi", 2, -13, 13)

    sigNum = ROOT.RooRealVar("sigNum", "sigNum", int(nEvents*0.2), 0, nEvents)
    bkgNum = ROOT.RooRealVar("bkgNum", "bkgNum", int(nEvents*0.8), 0, nEvents)

    pdfsig = ROOT.RooGenericPdf("sinpdf", "sinpdf", "S0 + (Sm*cos(2*pi*time/T + phi))", ROOT.RooArgSet(S0, Sm, time, periodfit, phasefit))
    # explicitly stating time as the observable - should be redundant because data histogram is set in dimension of time 
    pdfsig.getVal(ROOT.RooArgSet(time))
    
    pdfsig.fitTo(data, ROOT.RooFit.SumW2Error(True), ROOT.RooFit.Save(), ROOT.RooFit.Minimizer("Minuit2", "scan"), ROOT.RooFit.PrintLevel(1))
    
    nll = pdfsig.createNLL(data, ROOT.RooFit.Verbose())
    minimizer = ROOT.RooMinimizer(nll)
    minimizer.migrad()
    minimizer.hesse()
#    frame = time.frame(ROOT.RooFit.Title("toy time data"))
    frame = time.frame(ROOT.RooFit.Title("toy time data"))
#    nll.plotOn(frame)
    data.plotOn(frame)
    pdfsig.plotOn(frame)
    
    can = ROOT.TCanvas("c", "", 800, 500)
    can.cd()
    frame.Draw()
#    can.SaveAs("testdata_june21.png")
    can.SaveAs("testdata_june21.png")

if __name__ == "__main__":
    FitterESlice(10,15)
    
####### splitting S0 and Sm didn't work becausde Sm would go to a negative number: can't have a negative pdf
#pdfsig = ROOT.RooGenericPdf("sinpdf", "sinpdf", "S0 + Sm * (cos(2*TMath::Pi()*freq*time) + sin(2*TMath::Pi()*freq*time))", ROOT.RooArgSet(S0, Sm, freqfit, time))


#pdfsig = ROOT.RooGenericPdf("sinpdf", "sinpdf", "S0 + Sm * (cos(2*TMath::Pi()*freq*time + phase))", ROOT.RooArgSet(S0, Sm, freqfit, time, phasefit))
#pdfbkg = ROOT.RooGenericPdf("pol", "pol", "S0", ROOT.RooArgSet(S0))

#model = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(pdfsig, pdfbkg), ROOT.RooArgList(sigevents, bkgevents))

#model.fitTo(data, ROOT.RooFit.SumW2Error(False), ROOT.RooFit.PrintEvalErrors(2), ROOT.RooFit.Verbose(True))

