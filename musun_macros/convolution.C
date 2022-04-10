/*
Step 1: Just created a function that convolutes Klein Nishina cross section with a gaussian
Step 2: Create a fit function with three parameters - max kinetic energy calibrated with channel,width of the gaussian and normalization
Currently just to test if it gives same output as a T-max of 480, I tried the first parameter with 480
 */
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TRandom.h"
#include "TROOT.h"
// 137Cs 667 keV energy of photon emmited which has max KE as ~ 480 keV of scattered e



double crossSection(double energyK, double eMax)
{
  //The energies must be in keV
  double cos = 1-(energyK*511.0)/(eMax*(eMax-energyK));
  double P = 1/(1+(eMax/511.0)*(1-cos)); //ratio of scattered to incident photon energy
  double KlenCross = P*P*(P+(1/P)-1+cos*cos);
  //printf("Cross section is %f \n",KlenCross);
  // h4->Fill(180*TMath::ACos(cos)/3.14);
  return KlenCross;
   
}


double KlenGauss(double *KE, double *par)
{
  //par[0] is the calibration factor
  //par[1] is the width of the Gaussian (sigma)
  //par[2] is some normalization constant
  double eGamma = par[0];//Constant calibration factor
  int eKE = (int)KE[0]*par[0];
  
  /****** THIS FACTOR (par[0]) MUST BE MULTIPLED TO THE MAX KE ******/
  double MaxKE = 480*par[0]; // 480 keV is chosen for 137Cs - will change to 1030 for 60Co
  double sum = 0.0;
  for (int j=0; j<(int)MaxKE;j++)
    sum += TMath::Gaus(eKE,j,par[1])*crossSection((double)eKE,668.0);
  return sum/par[2];
}

void convolution()
{
  TF1 *g1 = new TF1("convolution",KlenGauss,0,668,3);
  g1->SetParameter(0,1.0);
  g1->SetParameter(1,20.0);
  g1->SetParameter(2,1.0);
  g1->Draw(); 
}
/*
void myfit(int run)
   {
     char name[20];
     sprintf(name,"try%d.root", run);
     TFile *_file0 = TFile::Open(name);     
     h1=hFadcNeutronEnergyVsChannel->ProjectionY("h1",1,1);
     TF1 *g1=gROOT->GetFunction("convolution");
     //TF1 *g1=new TF1("conv",convolution,0,668,3);
     double norm = h1->Integral(132,470);
     //g1->SetParameters(7.0,900.0,norm);
     h1->SetLineColor(kRed);
     h1->Fit("convolution","R","",1350.,4700.);
   }
*/
