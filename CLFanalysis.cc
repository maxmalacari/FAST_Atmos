// Max Malacari - 05/11/2018
// FAST CLF analysis
// ./CLFanalysis ./TestData/trace__2018_05_10_07h58m19s.root 0

// C/C++ classics
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sys/stat.h>
#include <numeric>
#include <time.h>
#include <utility>

// ROOT
#include "TMinuit.h"
#include "TROOT.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TVector.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

// Utility functions
bool LoadFitTrace(string fileName);
bool LoadFitTraceByPMT(string fileName);
bool LoadTelescopeEfficiency(string fileName);
bool LoadWLEfficiency(string fileName);
bool LoadDetectionEfficiency(string fileName);

// Fitting functions
double AtmosChi2(double Haer, double Laer, double norm, bool saveTrace=false);
void MinFunc(int &, double *, double &f, double *par, int );
double FitAtmosphere(double& Haer, double& dHaer, double& Laer, double& dLaer, double& norm, double& dNorm, bool fitRayleigh=false);

// Atmospheric functions
double GetMieTrans(double h1, double h2, double theta, double Haer, double Laer, double Hmix);
double GetRayTrans(double h1, double h2, double theta, double Hmol, double Lmol);
double GetScatteredFraction(double height, double theta, double Haer, double Laer, double Hmix, double Hmol, double Lmol);
double EvaluateHGFunc(double theta);
double EvaluateMolFunc(double theta);

// Global variables for fit and other functions
int gN_CLF = 0;
string gTraceTitle;
TH1F gTrace;
TH1F gPMTTrace[4];
TH1F gBestFit;
TH2F gTelEff;
TH1F gWLEff, gDetEff;

// Relative calibration constants for FAST telescope 2
// See John's FAST meeting presentation 5/11/2018
//const double gCalib[4] = {1.407, 1.042, 1.352, 1.};
const double gCalib[4] = {1., 1., 1., 1.};

// Fixed simulation parameters -------
#warning "Are all these parameters correct? Incl. pointing direction, energy, etc."
const double gLaserEnergy = 6.; // energy in mJ
const double gLaserWl = 355.; // wavelength in nm (note: molecular atmos. tuned to 355 nm)
const double gHmix = 0.; // aerosol mixing layer height in km
const double gHmol = 8.; // molecular scale height in km
const double gLmol = 14.2; // molecular horizontal attenuation length at sea level in km (Bucholz)
const double gZen = 0.0 * TMath::Pi()/180.; // rad - CLF pointing zenith
const double gAzi = 234.0 * TMath::Pi()/180.; // rad - CLF pointing azimuth
const double gCorex = 0.0; // km - CLF location
const double gCorey = 0.0; // km - CLF location
//const double gTelX = 17.0; // km - FAST position TA
//const double gTelY = -12.1; // km - FAST position TA
const double gTelX = -10.12; // km - FAST position Auger
const double gTelY = -23.95; // km - FAST position Auger
const double gHeight = 1.4; // km - ground height at TA/Auger
//const double gFovZen = 74.9 * TMath::Pi() / 180.0; // FAST pointing direction
//const double gFovAzi = (-59.8) * TMath::Pi() / 180.0;
const double gFovZen = 75.03 * TMath::Pi() / 180.0; // FAST pointing direction Auger
const double gFovAzi = (14.8) * TMath::Pi() / 180.0;
const double gJitterWidth = 22; // mean jitter width in 100 ns bins taken from old data
const int gSlideBins = 50; // number of bins either side to slide the template trace for fitting - nominally 100
#warning "Add laser energy jitter"
// -----------------------------------

// physical constants ------------------
TRandom3 gRandGen(time(NULL));
const double c = 299792458.; // Speed of light
const double h = 6.62607e-34; // Planck's constant
// ----------------------------

int main(int argc, char** argv){

  // Flag for simulated trace
  bool simTrace = false;

  // Set the initial fit values and step sizes
  double Haer = 1.; // first guess in km
  double Laer = 30.; // first guess in km
  double dHaer = 0.1; // initial step size in km
  double dLaer = 10.; // initial step size in km
  double norm = 1.; // allow absolute calib. to float
  double dNorm = 0.1; // initial step size

  if (argc < 3){
    cout << endl;
    cout << "Usage: ./CLFanalysis FAST_trace_file.root simTrace (0,1)" << endl;
    return 0;
  }

  simTrace = argv[2];

  // Load trace to fit using individual PMTs to adjust relative calibration
  string fitTrace = argv[1];
  if (simTrace){
    if(LoadFitTrace(fitTrace)){
      cout << "Loaded " << gN_CLF << " measured CLF shots from file " << fitTrace << endl;
    }
    else{
      cout << "Failed to load trace from file " << fitTrace << endl;
      return 0;
    }
  }
  else {
    if(LoadFitTraceByPMT(fitTrace)){
      cout << "Loaded " << gN_CLF << " measured CLF shots from file " << fitTrace << endl;
    }
    else{
      cout << "Failed to load trace from file " << fitTrace << endl;
      return 0;
    }
  }

  // Load spectrally-independent angular efficiency
  string telEffFile = "./Eff/telEff.dat";
  if (LoadTelescopeEfficiency(telEffFile)){
    cout << "Loaded telescope angular efficiency" << endl;
  }
  else {
    cout << "Failed to load telescope angular efficiency" << endl;
    return 0;
  }

  // Load wavelength efficiency (filter transmission and mirror reflectivity)
  string wlEffFile = "./Eff/wlEff.dat";
  if (LoadWLEfficiency(wlEffFile)){
    cout << "Loaded spectral efficiency" << endl;
  }
  else {
    cout << "Failed to load spectral efficiency" << endl;
    return 0;
  }

  // Load PMT detection efficiency
  string detEffFile = "./Eff/detEff.dat";
  if (LoadDetectionEfficiency(detEffFile)){
    cout << "Loaded PMT detection efficiency" << endl;
  }
  else {
    cout << "Failed to load PMT detection efficiency" << endl;
    return 0;
  }


  // Fit aerosol + molecular atmosphere!
  double chi2 = FitAtmosphere(Haer, dHaer, Laer, dLaer, norm, dNorm, false);
  
  // Output the result
  cout << "------ BEST FIT ------" << endl;
  cout << "Chi2: " << chi2 << endl;
  cout << "Haer: " << Haer << " +/- " << dHaer << " km" << endl;
  cout << "Laer: " << Laer << " +/- " << dLaer << " km" << endl;
  cout << "VAOD: " << 1./Laer*(Haer+gHmix) << endl;
  cout << "Calib: " << norm << endl;

  
  // Plot the results ------------------------------------------

  TCanvas summaryCanvas("summaryCanvas","test",1000,1300);
  summaryCanvas.Divide(1,2);

  summaryCanvas.cd(1);

  if (simTrace){
    gTrace.SetStats(0);
    gTrace.SetTitle(gTraceTitle.c_str());
    gTrace.SetXTitle("Time bins [100 ns]");
    gTrace.SetYTitle("N_{p.e.} / 100 ns");
    gTrace.GetXaxis()->CenterTitle();
    gTrace.GetYaxis()->CenterTitle();
    gTrace.GetXaxis()->SetTitleSize(0.05);
    gTrace.GetXaxis()->SetLabelSize(0.05);
    gTrace.GetYaxis()->SetTitleSize(0.05);
    gTrace.GetYaxis()->SetLabelSize(0.05);
    gTrace.GetXaxis()->SetTitleOffset(0.9);
    gTrace.GetYaxis()->SetTitleOffset(0.9);
    gTrace.SetLineColorAlpha(kBlack,0.9);
    gTrace.GetYaxis()->SetRangeUser(-6,30);
  }
  gTrace.Draw();
  for (int pmt=0; pmt<4; pmt++) gPMTTrace[pmt].Draw("same");

  TLegend leg1(0.58,0.58,0.88,0.88);
  leg1.AddEntry("gTrace","Summed trace","l");
  if (!simTrace){
    leg1.AddEntry("CLFTracePMT0","PMT 4","l");
    leg1.AddEntry("CLFTracePMT1","PMT 5","l");
    leg1.AddEntry("CLFTracePMT2","PMT 6","l");
    leg1.AddEntry("CLFTracePMT3","PMT 7","l");
  }
  leg1.SetBorderSize(0);
  leg1.Draw("same");

  summaryCanvas.cd(2);

  // Clone trace to plot untitled version on same canvas
  TH1F traceClone = *(TH1F*)gTrace.Clone();
  traceClone.SetTitle("");
  traceClone.Draw();

  // Generate best fit trace for plotting
  double tmp = AtmosChi2(Haer, Laer, norm, true);
  gBestFit.Draw("same");
  gBestFit.SetLineWidth(2);
  gBestFit.SetLineColorAlpha(kRed,0.7);

  TLegend leg2(0.58,0.78,0.88,0.88);
  leg2.AddEntry("gTrace","Measured trace","l");
  leg2.AddEntry("gBestFit","Best fit","l");
  leg2.SetBorderSize(0);
  leg2.Draw("same");

  TPaveText textBox(0.58,0.50,0.88,0.74);
  textBox.AddText(Form("H_{aer} = %.2f km",Haer));
  textBox.AddText(Form("L_{aer} = %.2f km",Laer));
  textBox.AddText(Form("Norm. = %.2f",norm));
  textBox.AddText(Form("VAOD = %.2f",1./Laer*(Haer+gHmix)));
  textBox.AddText(Form("#chi^{2}/ndf = %.2f",chi2/(1000.-2.-2.*gSlideBins)));
  textBox.Paint("NDC");
  textBox.SetBorderSize(0);
  textBox.SetFillColor(kWhite);
  textBox.SetTextFont(42);
  textBox.Draw();
  
  summaryCanvas.SaveAs("./summary.pdf");
  
  TCanvas fitResultCanvas;
  traceClone.Draw();
  gBestFit.Draw("same");
  leg2.Draw("same");
  textBox.Draw();

  fitResultCanvas.SaveAs("./result.pdf");

  return 0;
}

// Load a saved CLF trace - full telescope
bool LoadFitTrace(string fileName){

  // Open a measured laser trace
  TFile traceFile(fileName.c_str());
  if (!traceFile.IsOpen()) return false;
  
  gTrace = *(TH1F*)traceFile.Get("CLFTraceSum");
  // Number of stacked shots in the trace
  TVectorD *NshotVec = (TVectorD*)traceFile.Get("N_CLF");
  gN_CLF = (int)(*NshotVec)[0];

  traceFile.Close();

  return true;
}

// Load a saved CLF trace - individual PMTs
// Allows for rescaling of PMT relative calibration
bool LoadFitTraceByPMT(string fileName){

  // Open a measured laser trace
  TFile traceFile(fileName.c_str());
  if (!traceFile.IsOpen()) return false;
  
  // Number of stacked shots in the trace
  TVectorD *NshotVec = (TVectorD*)traceFile.Get("N_CLF");
  gN_CLF = (int)(*NshotVec)[0];

  // Open individual PMTs and add to total trace w/ rescaling applied
  gTrace = TH1F("gTrace","",1000,0,1000);
  for (int pmt=0; pmt<4; pmt++){
    stringstream plotName;
    plotName << "CLFTracePMT" << pmt;
    gPMTTrace[pmt] = *(TH1F*)traceFile.Get(plotName.str().c_str());
    gPMTTrace[pmt].Scale(1./gCalib[pmt]/gN_CLF);
    for (int bin=1; bin<=1000; bin++){
      gTrace.Fill(bin, gPMTTrace[pmt].GetBinContent(bin));
    }
  }

  // Retrieve the date and time of this stack
  // Not an ideal way of doing it
  TH1F gTrace2 = *(TH1F*)traceFile.Get("CLFTraceSum");
  gTraceTitle = gTrace2.GetTitle();
  

  TCanvas initCanvas;
  gTrace.SetStats(0);
  gTrace.SetTitle(gTraceTitle.c_str());
  gTrace.SetXTitle("Time bins [100 ns]");
  gTrace.SetYTitle("N_{p.e.} / 100 ns");
  gTrace.GetXaxis()->CenterTitle();
  gTrace.GetYaxis()->CenterTitle();
  gTrace.GetXaxis()->SetTitleSize(0.05);
  gTrace.GetXaxis()->SetLabelSize(0.05);
  gTrace.GetYaxis()->SetTitleSize(0.05);
  gTrace.GetYaxis()->SetLabelSize(0.05);
  gTrace.GetXaxis()->SetTitleOffset(0.9);
  gTrace.GetYaxis()->SetTitleOffset(0.9);
  gTrace.SetLineColorAlpha(kBlack,0.9);
  gTrace.Draw();
  //gTrace.SetLineWidth(2);
  gPMTTrace[0].SetLineColorAlpha(kRed,0.6);
  gPMTTrace[1].SetLineColorAlpha(kGreen+2,0.6);
  gPMTTrace[2].SetLineColorAlpha(kAzure+2,0.6);
  gPMTTrace[3].SetLineColorAlpha(kYellow+2,0.6);
  for (int pmt=0; pmt<4; pmt++){
    gPMTTrace[pmt].Draw("same");
  }
  TLegend leg(0.58,0.58,0.88,0.88);
  leg.AddEntry("gTrace","Summed trace","l");
  leg.AddEntry("CLFTracePMT0","PMT 4","l");
  leg.AddEntry("CLFTracePMT1","PMT 5","l");
  leg.AddEntry("CLFTracePMT2","PMT 6","l");
  leg.AddEntry("CLFTracePMT3","PMT 7","l");
  leg.SetBorderSize(0);
  leg.Draw("same");
  initCanvas.SaveAs("./input_trace.pdf");

  traceFile.Close();

  return true;
}


bool LoadTelescopeEfficiency(string fileName){

  gTelEff = TH2F("gTelEff","",399,-20,20,399,-20,20);

  // Load the telescope efficiency
  ifstream effFile(fileName.c_str());
  if (!effFile.is_open()) return false;

  double el, az, pmt1, pmt2, pmt3, pmt4;
  
  while(!effFile.eof()){
    effFile >> el >> az >> pmt1 >> pmt2 >> pmt3 >> pmt4;
    gTelEff.Fill(az, el, pmt1+pmt2+pmt3+pmt4);
  }

  TCanvas telEffCanvas;
  gTelEff.SetStats(0);
  gTelEff.Draw("COL4Z");
  telEffCanvas.SaveAs("./telEff.pdf");

  effFile.close();

  return true;
}

// Load the filter transmission and mirror reflectivity
bool LoadWLEfficiency(string fileName){

  gWLEff = TH1F("gWLEff","",341,260,600);

  // Load the telescope efficiency
  ifstream effFile(fileName.c_str());
  if (!effFile.is_open()) return false;

  double wl, eff, useless;
  
  while(!effFile.eof()){
    effFile >> wl >> eff >> useless >> useless;
    gWLEff.Fill(wl, eff/100.);
  }

  effFile.close();

  return true;
}

// Load the detection efficiency
bool LoadDetectionEfficiency(string fileName){

  gDetEff = TH1F("gDetEff","",17,200,600);

  // Load the telescope efficiency
  ifstream effFile(fileName.c_str());
  if (!effFile.is_open()) return false;

  double wl, eff, err;
  
  while(!effFile.eof()){
    effFile >> wl >> eff >> err;
    gDetEff.Fill(wl, eff/100.);
  }

  effFile.close();

  return true;
}

// Chi-squared function to compare simulated and measured trace
// Function runs the simulated trace along the measured trace in time
// to account for timing offset
double AtmosChi2(double Haer, double Laer, double norm, bool saveTrace){

  // Constants
  static const double laser[] = { sin(gZen) * cos(gAzi), sin(gZen) * sin(gAzi), cos(gZen) }; // laser axis unit direction (points up)
  static const double core[] = { gCorex-gTelX, gCorey-gTelY, 0.0 }; // laser position vector
  static const double center[] = { sin(gFovZen) * cos(gFovAzi), sin(gFovZen) * sin(gFovAzi), cos(gFovZen) }; // FAST FOV vector
  static const double detDistance = sqrt(core[0]*core[0] + core[1]*core[1] + core[2]*core[2]);
  static const int nBins = 1000; // 1000 * 100 ns bin trace

  
  // Simulate the CLF -----------------------------------------------------
  double nLaserPhotons = gLaserEnergy*1.e-3*gLaserWl*1.e-9/(c*h); // laser photons at ground level
       
  double binWidth = 100.e-9*c/1000.; // km - 30m ~ 100 ns binning at shower axis

  double prevLaserPhotons;
  
  double currentHeight = 0.;
  double signal;

  TH1F sim("sim","",1000,0,1000); // use a vector instead?
  for (int n = 0; n < nBins; n++) { // n steps up the laser axis

    prevLaserPhotons = nLaserPhotons;
    double nextHeight = currentHeight + binWidth;
    
    double Taer = GetMieTrans(currentHeight, nextHeight, TMath::Pi()/2., Haer, Laer, gHmix);
    double Tmol = GetRayTrans(currentHeight, nextHeight, TMath::Pi()/2., gHmol, gLmol);
    double T = Taer * Tmol;
    nLaserPhotons *= T;
    //cout << Taer << " " << Tmol << " " << nLaserPhotons << endl;
    double nRemovedLaserPhotons = prevLaserPhotons - nLaserPhotons;

    double binHeight = nextHeight - binWidth/2.; // Height in middle of bin
    double elevAngle =  atan(binHeight/detDistance); // Elevation in middle of bin
    double scattFrac = GetScatteredFraction(binHeight, elevAngle, Haer, Laer, gHmix, gHmol, gLmol);
    double Taer2Det = GetMieTrans(0., binHeight, elevAngle, Haer, Laer, gHmix);
    double Tmol2Det = GetRayTrans(0., binHeight, elevAngle, gHmol, gLmol);

    double step = currentHeight+binWidth/2.;
    double r[] = { core[0] + laser[0] * step, core[1] + laser[1] * step, core[2]
		   + laser[2] * step }; // point in the atmosphere
    double distance = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    double solidAngle = 1. / (distance * distance * 1e6); // 10^-6 so denominator is also in m^2
    
    double r_zen = (gFovZen - acos(r[2] / distance)) * 180.0 / TMath::Pi();
    double r_azi = atan(r[0]/r[1]) * 180.0 / TMath::Pi() + 59.8;
    #warning "Check this hard-coded 59.8 degrees..."

    double laserPhotons = nRemovedLaserPhotons*Taer2Det*Tmol2Det*scattFrac*solidAngle; // nLaserPhotons at this height
    
    if (r_zen < 20.0 && r_azi < 20.0 && r_zen > -20.0 && r_azi > -20.0) { // Falls within FOV
      double telEff = gTelEff.Interpolate(r_azi, r_zen);
      double wlEff = gWLEff.Interpolate(gLaserWl);
      double detEff = gDetEff.Interpolate(gLaserWl);
      signal = laserPhotons * telEff * wlEff * detEff * norm; // photons -> p.e.
      #warning "Check the raytracing map, and add mirror efficiency, etc. are we double counting something?"
      if (n<794) sim.SetBinContent(n+206, signal); // offset to put trace in middle of view
    }
    currentHeight = nextHeight;
    binWidth = 100.e-9*c/(1.+currentHeight/distance)/1000.; // in km - see notebook page 63
  }
  // ------------------------------- End CLF simulation

  // Stack N_CLF shots (include trigger and energy jitter)
  TH1F stackedSim("stackedSim","",1000,0,1000); // 1000 * 100 ns bins - use a vector instead
  gN_CLF = 10000; // stack many times to reduce fluctuations
  for (int k=0; k<gN_CLF; k++){
    double jitter = gRandGen.Gaus(0,gJitterWidth); // random jitter
    for (int bin=1; bin<=sim.GetNbinsX(); bin++){ // loop through bins in simulated trace
      if (bin+jitter <= stackedSim.GetNbinsX()) stackedSim.Fill(bin+jitter, sim.GetBinContent(bin)/gN_CLF);
    }
  }
  // -------------- Finished stacking -------------------

  
  // Calculate fit statistic for the best time shift by sliding the measured and simulated traces relative to each other
  double bestchi2 = 999999.;
  double bestShift = 0.;
  for (int shift = -gSlideBins; shift < gSlideBins; shift++){ // loop shifts
    double chi2 = 0;
    for (int bin=gSlideBins; bin<=(stackedSim.GetNbinsX()-gSlideBins); bin++){ // loop bins in template, minus the slideBins on either end (So we always compare the same number of bins)
      chi2 += pow( gTrace.GetBinContent(bin-shift) - stackedSim.GetBinContent(bin), 2.) / fabs(gTrace.GetBinContent(bin-shift));
    }
    if (chi2 < bestchi2){
      bestchi2 = chi2;
      bestShift = shift; // save this in case we want to save the trace
    }
  }

  if (saveTrace){
    gBestFit = TH1F("gBestFit","",1000,0,1000);
    for (int bin=gSlideBins; bin<=(stackedSim.GetNbinsX()-gSlideBins); bin++) gBestFit.SetBinContent(bin, stackedSim.GetBinContent(bin+bestShift));
  }

  //cout << bestShift << endl;
  
  return bestchi2;
  
}

// Minuit minimizing function
void MinFunc(int &, double *, double &f, double *par, int ){
  f = AtmosChi2(par[0], par[1], par[2], false);
  return; 
}

// Find best-fit atmosphere
double FitAtmosphere(double& Haer, double& dHaer, double& Laer, double& dLaer, double& norm, double& dNorm, bool fitRayleigh){

  int nmax_iterations = 100;

  //Minuit Initialization
  Int_t npar = 3;

  TMinuit* ptMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of npar

  ptMinuit->SetPrintLevel(1);//-1 no output - 1 std output

  double arglist[10];  int ierflg = 0;
  ptMinuit->mnexcm("SET NOW", arglist, 1, ierflg);//No warnings
  ptMinuit->SetErrorDef(1.);//1 for chi2, 0.5 for -logL
  ptMinuit->mnparm(0,"Haer", Haer, dHaer, 0.5, 1.5, ierflg);
  ptMinuit->mnparm(1,"Laer", Laer, dLaer, 1., 500., ierflg);
  ptMinuit->mnparm(2,"Norm", norm, dNorm, 0.1, 2.0, ierflg);

  // Fix normalization?
  ptMinuit->FixParameter(0);
  ptMinuit->FixParameter(1);
  ptMinuit->FixParameter(2);

  // Rayleigh atmosphere
  if (fitRayleigh){
    ptMinuit->mnparm(0,"Haer", 0.0000001, 0, 0.5, 1.5, ierflg);
    ptMinuit->mnparm(1,"Laer", 9999999., 0, 1., 500., ierflg);
    ptMinuit->FixParameter(0);
    ptMinuit->FixParameter(1);
  }
  
  ptMinuit->SetFCN(MinFunc);
  ptMinuit->SetMaxIterations(nmax_iterations);
  arglist[0]=1;//1 std, 2 try to improve (slower)
  ptMinuit->mnexcm("SET STR",arglist,1,ierflg);
  //#warning "Remove these arguments, just testing the tolerance."
  //arglist[0] = nmax_iterations;
  //arglist[1] = 10;
  ptMinuit->mnexcm("SIMPLEX",arglist,2,ierflg);
  //ptMinuit->mnexcm("MINUIT2",arglist,2,ierflg);	
  
  double funMin = 0.;
  double fedm, errdef;
  int npari, nparx, istat;
  ptMinuit->mnstat(funMin, fedm, errdef, npari, nparx, istat);
  ptMinuit->GetParameter(0, Haer, dHaer);
  ptMinuit->GetParameter(1, Laer, dLaer);
  ptMinuit->GetParameter(2, norm, dNorm);

  delete ptMinuit;

  return funMin;
}

// Get aerosol transmission between two heights at an angle theta to the vertical
double GetMieTrans(double h1, double h2, double theta, double Haer, double Laer, double Hmix)
{

  double tau;

  if (h1 < Hmix && h2 < Hmix)
    {
      tau = h2/Laer - h1/Laer;
    }
  if (h1 < Hmix && h2 >= Hmix)
    {
      tau = (Hmix-h1)/Laer + (Haer/Laer)*(1-exp(-(h2-Hmix)/Haer));
    }
  if (h1 >= Hmix && h2 >= Hmix)
    {
      tau = (Haer/Laer)*(exp(-(h1-Hmix)/Haer)-exp(-(h2-Hmix)/Haer));
    }
  
  return exp(-tau/sin(theta));
}

// Get molecular transmission between two heights at an angle theta to the vertical
double GetRayTrans(double h1, double h2, double theta, double Hmol, double Lmol)
{
  
  double deltaVMOD = (Hmol/Lmol)*( exp(-(h1+gHeight)/Hmol) - exp(-(h2+gHeight)/Hmol));

  return exp(-fabs(deltaVMOD)/sin(theta));
}

// Fraction of the beam scattered into dOmega due to aerosol and molecular scattering
double GetScatteredFraction(double height, double theta, double Haer, double Laer, double Hmix, double Hmol, double Lmol)
{

  double alphaAer;
  if (height < Hmix) alphaAer = (1./Laer);
  if (height >= Hmix) alphaAer = (1./Laer)*exp(-(height-Hmix)/Haer);
  double alphaMol = (1./Lmol)*exp(-(height+gHeight)/Hmol);
  double alphaTotal = alphaAer + alphaMol;
  
  double fractionAer = alphaAer/alphaTotal * EvaluateHGFunc(TMath::Pi()/2.+theta);
  double fractionMol = alphaMol/alphaTotal * EvaluateMolFunc(TMath::Pi()/2.+theta);
  double scatteredFraction = fractionAer + fractionMol;

  return scatteredFraction;
}

// Henyey-Greenstein scattering phase function
double EvaluateHGFunc(double theta)
{
  // Forward and backward scattering coefficients
  static double g = 0.6;
  static double f = 0.4;

  double dSigmadTheta = pow(1.-f,2.)/(4.*TMath::Pi()) * ( 1. / pow(1.+f*f-2.*f*cos(theta), 3./2.) + g*(3.*pow(cos(theta),2.) -1.)/(2.*pow(1.+f*f,3./2.)));
  
  return dSigmadTheta;
}

// Molecular scattering phase function
double EvaluateMolFunc(double theta)
{
  double dSigmadTheta = (3./(16.*TMath::Pi())) * (1.+cos(theta)*cos(theta));

  return dSigmadTheta;
}
