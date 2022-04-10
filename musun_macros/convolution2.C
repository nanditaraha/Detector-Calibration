/*
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

TF1 *f1, *f2;

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
//Unconvoluted function of cross section only
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
       //return sum/sumgauss*par[2]*2.4*par[0]/1030+par[3];//For Co
       return sum/sumgauss*2.4*par[2]*par[0]/480+par[3];//For Cs
     else 
       return sum*2.4*par[2]*par[0]/480+par[3];
       //return sum*2.4*par[2]*par[0]/1030+par[3];
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
       return sum/sumgauss*par[2]*par[0]/480+par[3];//For Cs
     //return sum/sumgauss*par[2]*par[0]/1030+par[3];
     else
       return sum*par[2]*par[0]/480+par[3];// For Cs
     //return sum*par[2]*par[0]/1030+par[3];//For Co
   
}

void convolution2(int run, int run1, int channel)
{
  char name[20], name1[20];
  sprintf(name,"hist%d.root", run);
  sprintf(name1,"hist%d.root", run1);
  TFile *_file0 = TFile::Open(name);    
  TFile *_file1 = TFile::Open(name1);    
  TH3D *h3 = (TH3D*)_file0->Get("hFadcNeutronEnergyVsChannel");
  TH3D *h2 = (TH3D*)_file1->Get("hFadcNeutronEnergyVsChannel");
  TH1D *h = (TH1D*)h3->ProjectionY("h",channel,channel);
  TH1D *h1 = (TH1D*)h2->ProjectionY("h1",channel,channel);
  h->SetTitle("Pulse Energy ND14 with 137Cs ");
  sigma1 = 60.0;

  sigma = 70.0;// sigma = 40; // (3320,40.0,5500,200) gives an ok shape for Cs
  double channelKE = 2379.0;
  double channelKE1 = 2318.0;
  f1 = new TF1("f1","TMath::Gaus(x,[0],[1])",-3*sigma,+3*sigma);
  f1->SetParameters(0,sigma);  
  f2 = new TF1("f2","TMath::Gaus(x,[0],[1])",-3*sigma1,+3*sigma1);
  f2->SetParameters(0,sigma1);
  //These are for Cs------------- 
  TF1 *g1 = new TF1("convolution2",KlenGauss,0,9000,4);// For Cs spectra
  TF1 *g2 = new TF1("convolution2",KlenGauss,0,9000,4);
  double norm = h->GetBinContent(237);//For Cs
  double norm1 = h1->GetBinContent(231);//For Cs
  g2->SetParameters(channelKE1,sigma1,norm1,300.0);
  g1->SetParameters(channelKE,sigma,norm,64.0);// These parameters give a good fit for Cs
  //g1->FixParameter(3,91.0);
  g1->SetParNames("Channel for KE","Sigma","Norm","Background");
  g2->SetParNames("Channel for KE","Sigma","Norm","Background");
  h1->Draw();
  h->Draw("same");

 
  h->Fit("convolution2","R","",1900,3600);// For Cs 
  h1->Fit("convolution2","R","",1900,3600);// For Cs 
  //g1->ReleaseParameter(3); 
  //h1->Fit("convolution1","R","",2250,4000);// For Cs 

  // These are for Co -------------
  /*TF1 *g1 = new TF1("convolution1",KlenGauss,0,7000,4);// For Co spectra 
  double norm = h1->GetBinContent(520);// For Co  
  g1->SetParameters(channelKE,sigma,norm,0.0);
  g1->SetParNames("Channel for KE","Sigma","Norm","Background");
  
  h1->Fit("convolution1","R","",4500,6300);// For Co 
  TF1 *g2 = new TF1("klen1",Klen,0,7000,4);
  g2->SetParameters(channelKE,sigma1,norm,0.0);
  */
  h->SetLineColor(kRed);
  //h->SetLineColor(kMagenta);
  g2->SetLineColor(kBlue);
  g1->SetLineColor(kGreen);
  //printf("%f\n",h1->GetBinContent(340));
  
  g2->Draw("same");
  g1->Draw("same");
 
  h->Draw("same");
  h1->Draw("same");
}
