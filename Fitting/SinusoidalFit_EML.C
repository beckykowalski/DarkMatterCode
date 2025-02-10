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
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
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
  
  std::vector<std::vector<double> > Data = ReadData(datafile);  
  //  std::vector<std::vector<double> > Data = ReadData("DC30_Sm10_phi0_T365_poisson.txt");  

  TH1D* hist = new TH1D("hist", "", 1825, 0, 1825);
  hist->Rebin(binWidth);

  int nbins = hist->GetNbinsX();
  std::cout << "nbins = " << nbins << std::endl;
  for(int i=0; i<nbins;i++){
    hist->SetBinContent(i+1, Data[1][i]);
  }

  RooRealVar x("x", "time (days)", 0, 1825);
  RooDataHist* data = new RooDataHist("data", "data", RooArgSet(x), Import(*hist));  

  // a1 and a2: modulating amplitude components of PDF
  RooRealVar a1("a1", "a1", .001, -2, 2);
  RooRealVar a2("a2", "a2", .001, -2, 2);
  // N: total number of expected events in distribution (for extended likelihood)
  RooRealVar N("N", "N", 54750, 1000, 100000);

  // angular frequency of expected modulation 2*pi/365 - units: (days)^-1
  double omega = 0.0172;
  
  RooGenericPdf sig("sig", "sig", "1 + a1*cos(0.0172*x)+a2*sin(0.0172*x)", RooArgList(a1, x, a2));
  RooExtendPdf extendedsig("extendedsig", "extendedsig", sig, N);

  // fit an extended maximum likelihood with varaibles N, a1, and a2
  RooFitResult *result = extendedsig.fitTo(*data, Strategy(2), Verbose(false), PrintLevel(0), Hesse(true), Save(true));
  //  RooAbsReal *nll = extendedPdf.createNLL(*data);
  RooPlot *frame = x.frame();
  data->plotOn(frame);
  extendedsig.plotOn(frame);
  
  TCanvas *can = new TCanvas("can", "", 800, 500);
  can->cd();
  frame->Draw();
  can->SaveAs(Form("flat_%d.png", it));
  can->SetBatch(true);
  std::vector<double> modulating_info;

  // begin converting fit parameters to scaled expected fit parameters (b, a1_without_normalization, a2_without_normalization)
  double tf = hist->GetBinCenter(1825);

  double a1v = a1.getVal();
  double a2v = a2.getVal();
  double a1E = a1.getError();
  double a2E = a2.getError();
  double Nv = N.getVal();
  double NE = N.getError();

  // get b (and uncertainty on b) from integration of N, a1, and a2
  double b = (Nv - (a1v/omega)*(sin(omega*tf)) + (a2v/omega)*(cos(omega*tf) - 1) )/tf;
  double db_dN = 1/tf;
  double db_da1 = -sin(omega*tf)/(omega*tf);
  double db_da2 = (1/(omega*tf))*(cos(omega*tf) - 1);
  double db = sqrt((NE*db_dN)*(NE*db_dN) + (a1E*db_da1)*(a1E*db_da1) + (a2E*db_da2)*(a2E*db_da2));
  std::cout << "db_dN component = " << db_dN << std::endl;
     
  // get a1_without_normalization (and uncertainty on a1_without_normalization) from integration of N, b and a2
  double a1N = omega/sin(omega*tf) * (Nv - b*tf + (a2v/omega)*(cos(omega*tf) - 1) );
  double da1_dN = omega/sin(omega*tf);
  double da1_db = (omega*tf)/sin(omega*tf);
  double da1_da2 = (cos(omega*tf)-1)/sin(omega*tf);
  double da1 = sqrt((NE*da1_dN)*(NE*da1_dN) + (a2E*da1_da2)*(a2E*da1_da2) + (db*da1_db)*(db*da1_db));
  std::cout << "da1_dN component = " << da1_dN << std::endl;

  // get a2_without_normalization (and uncertainty on a2_without_normalization) from integration of N, b and a1
  double a2N = (omega/(cos(omega*tf)-1)) * (Nv - b*tf - (a1v/omega)*sin(omega*tf));
  double da2_dN = omega/(cos(omega*tf)-1);
  double da2_db = (omega*tf)/(cos(omega*tf)-1);
  double da2_da1 = sin(omega*tf)/(cos(omega*tf)-1);
  double da2 = sqrt((NE*da2_dN)*(NE*da2_dN) + (a1E*da2_da1)*(a1E*da2_da1) + (db*da2_db)*(db*da2_db));
  std::cout << "da2_dN component = " << da2_dN << std::endl;

  double Sm = sqrt(a1N*a1N + a2N*a2N);
  double dSm = sqrt( (da1*a1N/sqrt(a1N*a1N + a2N*a2N)) + (da2*a2N/sqrt(a1N*a1N + a2N*a2N)) );

  double phi = atan(a2N/a1N);
  double dphi = sqrt( (da1*-a2N/(a1N*a1N + a2N*a2N)) + (da2*a1N/(a1N*a1N + a2N*a2N)) );
  
  
  std::cout << "b = " << b << " +/- " << db << std::endl;
  std::cout << "a1 = " << a1N << " +/- " << da1 << std::endl;
  std::cout << "a2 = " << a2N << " +/- " << da2 << std::endl;
  std::cout << "Sm = " << Sm << " +/- " << dSm << std::endl;
  std::cout << "phi = " << phi << " +/- " << dphi << std::endl;
  
  double fitStatus = (double)result->status();
 
  
}

int RunPlotter(){

  //  for(int i=0; i<1000; i++){
  for(int i=0; i<1; i++){
  
    //   std::vector<double> val = RooPlotter(1, Form("../SimulationForFitting/flat_data/FlatData_%d.txt",i), i);
    //    std::vector<double> val = RooPlotter(1, Form("../SimulationForFitting/Sm13_phase_data/Sm13_phaseData_%d.txt",i), i);
    RooPlotter(1, Form("../SimulationForFitting/Sm10_data/Sm10Data_%d.txt",i), i);
  }

  return 0;
}
