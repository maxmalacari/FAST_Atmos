// Turn single Offline-simulated laser trace into a stacked set of N_clf traces with included noise
// for use with the CLF analysis software
// M.M. - 8/2/2019

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

int main(int argc, char** argv){

  int nCLF = 250; // number of CLF shots to stack
  double noiseWidth = 98./5.; // RMS noise / 100 ns to add to trace
  double jitterWidth = 22.; // RMS trigger jitter in units of 100 ns

  TRandom3 randGen(time(NULL));

  if (argc < 2){
    cout << endl;
    cout << "Usage: ./makeDataFromSim FASTsim_file.root" << endl;
    return 0;
  }

  // Load trace to fit using individual PMTs to adjust relative calibration
  string fileName = argv[1];

  // Load the simulated FAST trace
  TFile simFile(fileName.c_str());
  if (!simFile.IsOpen()){
    cout << "Failed to load FAST simulation file " << fileName << endl;
    return 0;
  }

  TH1F simTrace = *(TH1F*)simFile.Get("traceOut");
  simFile.Close();

  // New stacked trace
  TH1F CLFTraceSum("CLFTraceSum","",1000,0,1000);

  for (int n=0; n<nCLF; n++){
    double jitter = randGen.Gaus(0,jitterWidth);
    for (int i=1; i<= simTrace.GetSize()-2; i++){
      CLFTraceSum.Fill(simTrace.GetBinCenter(i)+jitter, simTrace.GetBinContent(i)+randGen.Gaus(0,noiseWidth));
    }
  }
  CLFTraceSum.Scale(1./double(nCLF));

  TCanvas testCanvas("testCanvas","",1200,1400);
  testCanvas.Divide(1,2);
  testCanvas.cd(1);
  simTrace.Draw();
  testCanvas.cd(2);
  CLFTraceSum.Draw();
  testCanvas.SaveAs("test.pdf");

  // Save in format suitable for CLF analysis code
  TVectorD N_clf(1);
  N_clf[0] = nCLF;
  TFile fOut("./simtrace.root","RECREATE");
  CLFTraceSum.Write();
  N_clf.Write("N_CLF");
  fOut.Close();
  
  // Create a stack of N copies of the simulation, adding a mean timing jitter of 22 bins (maybe)
  // add noise to each trace
  // stack traces
  // save in same format as real data

  


  return 0;

}
