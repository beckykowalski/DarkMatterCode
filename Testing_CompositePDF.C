#include <iostream>
#include <fstream>
#include <math.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMath.h"

#include "RooCategory.h"
#include "RooMappedCategory.h"
#include "RooMultiCategory.h"
#include "RooSuperCategory.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
//#include "RooConst.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooUniform.h"
#include "RooAddModel.h"

using namespace std;
using namespace RooFit;

//std::vector< std::vector<double> > ReadData(const std::string& datafile) {
std::vector<std::vector<double> > ReadData(const char* datafile) {
  std::ifstream textfile(datafile);

  std::vector< std::vector<double> > output;

  
  if (!textfile.is_open()) {
    throw std::runtime_error("Could not open file");
  }

  std::vector<double> t;
  std::vector<double> r;
  std::string line;

  // Read each line from the file
  while (std::getline(textfile, line)) {
    // Skip the header line
    if (line == "time,rate") {
      continue;
    }

    std::istringstream ss(line);
    std::string time_str, rate_str;

    // Split the line by comma
    std::getline(ss, time_str, ',');
    std::getline(ss, rate_str, ',');

    t.push_back(std::stod(time_str));
    r.push_back(std::stod(rate_str));

  }

  textfile.close();
  output.push_back(t);
  output.push_back(r);
  return output;
}

//"../SimulationForFitting/flat_data/FlatData_188.txt"
void RooPlotter(int binWidth, const char* datafile, int it) {

  gROOT->SetBatch(true);
  double avgRate = 0;
  
  std::vector<std::vector<double> > Data = ReadData(datafile);  
  //  std::vector<std::vector<double> > Data = ReadData("DC30_Sm10_phi0_T365_poisson.txt");  

  TH1D* hist = new TH1D("hist", "", 1825, 0, 1825);
  hist->Rebin(binWidth);

  int nbins = hist->GetNbinsX();
  std::cout << "nbins = " << nbins << std::endl;
  for(int i=0; i<nbins;i++){
    hist->SetBinContent(i+1, Data[1][i]);
    avgRate += Data[1][i];    
  }

  std::cout << "number of entries in histogram is " << avgRate << std::endl;
  
  avgRate /= nbins;

  // initialize fit parameter seeds
  double maxHist = hist->GetMaximum();
  double minHist = hist->GetMinimum();
  double GuessBinCont = hist->GetBinContent(0);

  RooRealVar x("x", "time (days)", 0, 1825);
  RooDataHist* data = new RooDataHist("data", "data", RooArgSet(x), Import(*hist));  

  double TwoPi = 2*TMath::Pi();
  
  RooRealVar Sm("Sm", "Sm", .001, -1, 1);
  RooRealVar a1("a1", "a1", .001, -1, 1);
  RooRealVar a2("a2", "a2", .001, -1, 1);
  RooRealVar phi("phi", "phi", .001, 0, 365);
  RooRealVar frac("frac", "frac", .8, .5, 1.);
  RooRealVar N("N", "N", 52000, 10000, 100000);
  double omega = 0.0172;

  RooUniform bkg("bkg", "bkg", x);
  //  RooGenericPdf sig("sig", "sig", "Sm*(1 + cos(0.0172*(x + phi)))", RooArgList(Sm, x, phi));
  RooGenericPdf sig("sig", "sig", "1 + a1*cos(0.0172*x) + a2*sin(0.0172*x)", RooArgList(a1, x, a2));

  RooAddPdf sumPdf("sumPdf", "sumPdf", RooArgList(bkg, sig), RooArgList(frac));
  RooExtendPdf extendedsig("extendedsig", "extendedsig", sumPdf, N);
  
  
  RooFitResult *result = extendedsig.fitTo(*data, Strategy(2), Verbose(false), PrintLevel(0), Hesse(true), Save(true));
  //  RooAbsReal *nll = extendedPdf.createNLL(*data);
  RooPlot *frame = x.frame();
  data->plotOn(frame);
  extendedsig.plotOn(frame), RooFit::LineColor(kBlue-7);
  extendedsig.plotOn(frame, RooFit::Components(bkg), RooFit::LineStyle(kDashed), RooFit::LineColor(kViolet+1));
  extendedsig.plotOn(frame, RooFit::Components(sig), RooFit::LineStyle(kDashed), RooFit::LineColor(kTeal+1));

  
  TCanvas *can = new TCanvas("can", "", 800, 500);
  can->cd();
  frame->Draw();
  can->SaveAs(Form("flat_%d.png", it));
  can->SetBatch(true);
  std::vector<double> modulating_info;

  
  double tf = hist->GetBinCenter(1825);
  
  double fitStatus = (double)result->status();
 
  
}

int RunPlotter(){

  //  for(int i=0; i<1000; i++){
  for(int i=0; i<1; i++){
  
    //   std::vector<double> val = RooPlotter(1, Form("../SimulationForFitting/flat_data/FlatData_%d.txt",i), i);
    //    std::vector<double> val = RooPlotter(1, Form("../SimulationForFitting/Sm13_phase_data/Sm13_phaseData_%d.txt",i), i);
    //    RooPlotter(1, Form("../SimulationForFitting/Sm10_data/Sm10Data_%d.txt",i), i);
    RooPlotter(1, Form("../SimulationForFitting/Sm10_data/Sm10Data_%d.txt",i), i);
  }

  return 0;
}
