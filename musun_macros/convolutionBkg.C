/*THIS IS THE LATEST UPDATED VERSION
Step 1: Just created a function that convolutes Klen Nishina cross section with a gaussian
The Gaussian is centered around zero spreads 3*sigma both sides but convolution starts with the content in the centre bin of gaussian with 1st bin of Klen Nishina distribution + bins of gaussian move left (till +-3*sigma) and Klen N move right
2nd bin of convolution = Move the Gaussian by a bin and repeat like before till end of the Gaussian (+- 3*sigma)
Step 2: Create a fit function with three parameters - max kinetic energy calibrated with channel,width of the gaussian and normalization
Currently just to test if it gives same output as a T-max of 480, I tried the first parameter with 480

 */

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
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

double sigma, sigma1;

TF1 *f1, *f2,*f3;
/* This was used before PeterK's correction
double crossSection(double energyK, double eMax)
{
  //The energies must be in keV
  if (energyK>2*eMax*eMax/(511.0+2*eMax)) 
    return 0; 
  double cos = 1-(energyK*511.0)/(eMax*(eMax-energyK));
  double P = 1/(1+(eMax/(511.0))*(1-cos)); //ratio of scattered to incident photon energy
  double KlenCross = P*P*(P+(1/P)-1+cos*cos);
  return KlenCross;
   
}
*/
//This is modified after PeterK's correction
double crossSection(double energyK, double eMax)
{
  //The energies must be in keV
  if (energyK>2*eMax*eMax/(511.0+2*eMax)) 
    return 0; 
  // double cos = 1-(energyK*511.0)/(eMax*(eMax-energyK));
  double P = energyK*511.0/(eMax*(eMax-energyK));
  double KlenCross = 511.0*(2-2*P+P*P+energyK/(eMax*(eMax-energyK)))/(eMax*eMax);
  return KlenCross;
   
}
//Unconvoluted function of cross section only so, this has negligible sigma from f2
double Klen(double *KE, double *par)
{
  int eKE = (int)(KE[0]*480)/par[0];//3298 is the channel # read from the graph for max. KE = 480 for Cs
  //int eKE = (int)(KE[0]*1030)/par[0];//7554 is the channel # read from the graph for max. KE = 1030 for Co
  double eGamma = 667.0;// For Cs
  //double eGamma = 1250.0;// For Co
  double sum = 0.0, sumgauss = 0.0;
  double gauss;
  /****** THIS FACTOR (par[0]) MUST BE MULTIPLED TO THE MAX KE ******/
  /* original value of KEmax = 480 and photon enery =668*/
  // for (int j=0; j<(int)MaxKE;j++)
  // {
     sum = 0.0;
     sumgauss = 0.0;
     for(int i=0;i<(int)eGamma;i++)
       {
	 int k = -i+(int)eKE;
	 if (k<-3*par[1]) continue;
	 if (k>+3*par[1]) continue;
	   {
	     
	     gauss = f2->Eval(k);
             sumgauss += gauss;//This is done to normalize the Gaussian for parts which were not normalized
	     sum += gauss*crossSection((double)i,eGamma);
	   }
	   
	     
       }
     //printf("Sumgauss:%f \n",sumgauss);
     if(sumgauss>1.0)
       //return sum/sumgauss*402*par[2]*2.4*par[0]/1030+par[3];//For Co
       return sum/sumgauss*200*par[2]*par[0]/480+par[3];//For Cs
     else 
       return sum*200*par[2]*par[0]/480+par[3]; // For Cs
     //return sum*402*par[2]*par[0]/1030+par[3];// For Co
}

double back(double *KE, double *par){
 
  return par[0]*TMath::Exp(par[1]*KE[0])+par[2]*TMath::Exp(par[1]*KE[3])+par[4]*TMath::Exp(par[1]*KE[5]);
  //return par[0]/KE[0]+par[1]/(KE[0]*KE[0]);
 
  //return 2000+par[0]*KE[0]+par[1]*(KE[0]*KE[0]);
}

double KlenGauss(double *KE, double *par)
{
  //par[0] is the max KE channel
  //par[1] is the width of the Gaussian (sigma)
  //par[2] is some normalization constant gives the count in the channel for max KE
  //par[3] background

  int eKE = (int)(KE[0]*480)/par[0];//3298 is the channel # read from the graph for max. KE = 480 for Cs
  //int eKE = (int)(KE[0]*1030)/par[0];//7554 is the channel # read from the graph for max. KE = 1030 for Co
  double eGamma = 667.0;// For Cs
  //double eGamma = 1250.0;// For Co
  double sum = 0.0, sumgauss = 0.0;
  double gauss;
  /****** THIS FACTOR (par[0]) MUST BE MULTIPLED TO THE MAX KE ******/
  /* original value of KEmax = 480 and photon enery =668*/
  // for (int j=0; j<(int)MaxKE;j++)
  // {
     sum = 0.0;
     sumgauss = 0.0;
     for(int i=0;i<(int)eGamma;i++)
       {
	 int k = -i+(int)eKE;
	 if (k<-3*par[1]) continue;
	 if (k>+3*par[1]) continue;
	   {
	     
	     gauss = f1->Eval(k);
             sumgauss += gauss;//This is done to normalize the Gaussian for parts which were not normalized
	     sum += gauss*crossSection((double)i,eGamma);
	   }
	   //if (i<10)
	   //printf("i, k, gaus, cs, sum eKE = %d %d %f %f %f %d\n",i,k,gauss,crossSection((double)i,667.0),sum,(int)eKE);
	  
       }
     if (sumgauss > 1.0)
       return sum/sumgauss*200*par[2]*par[0]/480+par[3];//For Cs
     //return sum/sumgauss*402*par[2]*par[0]/1030+par[3];
     else
       return sum*200*par[2]*par[0]/480+par[3];// For Cs
     //return sum*402*par[2]*par[0]/1030+par[3];//For Co
   
}

double fitBack(double *KE, double *par){
  return (back(KE, par)+KlenGauss(KE, &par[6]));
}
  

void convolutionBkg(int run, int channel, double p0=2500., double p1=-1.15e-3, double p2=1500., double p3=-3.8e-4, double p4=1000., double p5=3.9e-4, double sigma=45, double channelKE=3100)
//void convolutionBkg(int run, int channel, double p0=-3.2, double p1=8.2e-4 )
{
  char name[20];
  sprintf(name,"hist%d.root", run);
  TFile *_file0 = TFile::Open(name);    
  TH3D *h3 = (TH3D*)_file0->Get("hFadcNeutronEnergyVsChannel");
  TH1D *h1 = (TH1D*)h3->ProjectionY("h",channel,channel);
  h1->SetTitle("Pulse Energy NU6 with 137Cs ");
  sigma1 = 0.01;

  //sigma = 90.0;// 
  //sigma = 65; // (3320,40.0,5500,200) gives an ok shape for Cs
  //double channelKE = 4970.0; for Co
  //double channelKE = 3200.0;
  f1 = new TF1("f1","TMath::Gaus(x,[0],[1])",-3*sigma,+3*sigma);
  f1->SetParameters(0,sigma);  
  //f1->Draw();
  f2 = new TF1("f2","TMath::Gaus(x,[0],[1])",-3*sigma1,+3*sigma1);
  //f3 = new TF1("f3","expo",0.,10000.);
  f2->SetParameters(0,sigma1);
  //f3->SetParameters(8.4,-1.2e-3);
  //f3->Draw();
  //These are for Cs------------- 
  //TF1 *g1 = new TF1("convolution1",KlenGauss,0,10000,4);// For Cs spectra
  TF1 *g1 = new TF1("convolution1",fitBack,0,10000,10);
  
  TF1 *g3 = new TF1("back",back,0,10000,6);
  TF1 *g4 = new TF1("conv",KlenGauss,0,10000,4);
  TF1 *g2 = new TF1("klen",Klen,0,10000,4);
  //g3->FixParameter(0,p0);
  //g3->FixParameter(1,-1.6e-3);
  
  double norm = h1->GetBinContent(320);
  double par[9];
  g2->SetParameters(channelKE,sigma1,norm,10.0);
  g4->SetParameters(channelKE,sigma,norm,10.0);
  g1->SetParameters(p0,p1,p2,p3,p4,p5,channelKE,sigma,norm,28.0);// These parameters give a good fit for Cs
  //Fix sigma and try to fit the bump for NU6 and ND6
  g3->SetParameters(p0,p1,p2,p3,p4,p5);
  g1->SetParNames("Amp","Slope1","Channel for KE","Sigma","Norm","Background");
  h1->Draw();
  /*  
      
  h1->Fit("back","R");
  h1->Fit("conv","R+");// For Cs 
  g3->GetParameters(par);
  g4->GetParameters
  */
  //g1->ReleaseParameter(3); 
  // h1->Fit("convolution1","R","",500,4200);// For Cs 

  // These are for Co -------------
  /*TF1 *g1 = new TF1("convolution1",KlenGauss,0,7000,4);// For Co spectra 
  double norm = h1->GetBinContent(520);// For Co  
  g1->SetParameters(channelKE,sigma,norm,0.0);
  g1->SetParNames("Channel for KE","Sigma","Norm","Background");
  
  h1->Fit("convolution1","R","",4500,6300);// For Co 
  TF1 *g2 = new TF1("klen1",Klen,0,7000,4);
  g2->SetParameters(channelKE,sigma1,norm,0.0);
  */

  h1->SetLineColor(kRed);
  g4->SetLineColor(kBlue);
  g1->SetLineColor(kGreen);
  // f3->SetLineColor(kGreen);
  //printf("%f\n",h1->GetBinContent(340));
  
  g4->Draw("same");
  g1->Draw("same");
  g2->Draw("same");
  //f1->Draw("same");
  g3->Draw("same");
}
