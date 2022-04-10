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

double sigma, sigma1, normalize =1, normalize1 = 1;

TF1 *f1, *f2;
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
double crossSection(double energyK, double eMax) //simple C++ function
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
double Klen(double *KE, double *par)  //Root function
{
  int eKE = (int)(KE[0]*480)/par[0];//3298 is the channel # read from the graph for max. KE = 480 for Cs
  //int eKE = (int)(KE[0]*1030)/par[0];//7554 is the channel # read from the graph for max. KE = 1030 for Co
  double eGamma = 662.0;// For Cs
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
       return sum/sumgauss*150*par[2]*par[0]/480+par[3];//For Cs a factor of 150 for channel 7- default constructor
       //return sum/sumgauss*2800*par[2]*par[0]/1030+par[3];//For Co
       //return normalize*sum/sumgauss*par[2]*par[0]/480+par[3];// For Cs with peak bin as 340 for channel 7
     else 
       return sum*150*par[2]*par[0]/480+par[3]; // For Cs
       //return normalize*sum*par[2]*par[0]+par[3]; // For Cs
     //return sum*2800*par[2]*par[0]/1030+par[3];// For Co
}

double KlenGauss(double *KE, double *par)
{
  //par[0] is the max KE channel
  //par[1] is the width of the Gaussian (sigma)
  //par[2] is some normalization constant gives the count in the channel for max KE
  //par[3] background

  int eKE = (int)(KE[0]*480)/par[0];//3298 is the channel # read from the graph for max. KE = 480 for Cs
  //int eKE = (int)(KE[0]*1030)/par[0];//7554 is the channel # read from the graph for max. KE = 1030 for Co
  //double eGamma = KE[0]*667.0/par[0];// For Cs
  double eGamma = 662.0;// For Cs
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
       //return normalize1*sum/sumgauss*par[2]*par[0]+par[3];//For Cs
       return sum/sumgauss*150*par[2]*par[0]/480+par[3];//For Cs - default constructor ch 7
       //return sum/sumgauss*2800*par[2]*par[0]/1030+par[3];// For Co
     else
       //return normalize1*sum*par[2]*par[0]+par[3];// For Cs
       return sum*par[2]*150*par[0]/480+par[3];// For Cs - default constructor ch 7
     //return sum*2800*par[2]*par[0]/1030+par[3];//For Co
   
}

void convolution1(int run=32541, int channel=7, double sigma=60., double channelKE=4100., double binKE=218., double back=5.)
//convolution1(32541, 3,60.,2330.,236.,5.)
//convolution1(32541, 7,60.,4100.,218.,5.)
//convolution1(36868, 1,60.,3630.,218.,5.)

{
  char name[20];
  sprintf(name,"try%d.root", run);
  TFile *_file0 = TFile::Open(name);    
  TH3D *h3 = (TH3D*)_file0->Get("hFadcNeutronEnergyVsChannel");
  TH1D *h1 = (TH1D*)h3->ProjectionY("h",channel,channel);
  h1->SetTitle("Pulse Energy ND14 with 60Co ");
  sigma1 = 0.01;//Used for plotting only Klein Nishina

  //sigma = 90.0;// 
  //sigma = 50; // (2320,40.0,5500,200) gives an ok shape for Cs
  //double channelKE = 4970.0; for Co
  //double channelKE = 2300.0; for cs;
  f1 = new TF1("f1","TMath::Gaus(x,[0],[1])",-3*sigma,+3*sigma);
  f1->SetParameters(0,sigma);  
  f1->SetLineColor(kMagenta);
  f1->Draw();
  f2 = new TF1("f2","TMath::Gaus(x,[0],[1])",-3*sigma1,+3*sigma1);
  f2->SetParameters(0,sigma1);
  //These are for Cs------------- 
  TF1 *g1 = new TF1("convolution1",KlenGauss,0,5500,4);// For Cs spectra
  TF1 *g4 = new TF1("convolution2",KlenGauss,0,5500,4);// For Cs spectra
  TF1 *g2 = new TF1("klen1",Klen,0,5500,4);
  TF1 *g3 = new TF1("klen2",Klen,0,5500,4);   
  double norm = h1->GetBinContent(binKE);//For Cs
  printf("Before:%f and klein %f normalize %f\n",h1->GetBinContent(340),g2->Eval(h1->GetBinCenter(340)),normalize);
  g1->SetParameters(channelKE,sigma,norm,back);
  g2->SetParameters(channelKE,sigma1,norm,back);  
  double x = g2->Eval(h1->GetBinCenter(340));
  normalize = h1->GetBinContent(340)/x ;//for Cs bin # 340 is a good const to normalize areas for channel 7.
  printf("x %f\n",x);
  g3->SetParameters(channelKE,sigma1,norm,back);
  normalize1 = h1->GetBinContent(340)/x ;//for Cs bin # 340 is a good const to normalize areas.
  g4->SetParameters(channelKE,sigma,norm,back);// These parameters give a good fit for Cs
  printf("Later:x %f\n",x);
  //Fix sigma and try to fit the bump for NU6 and ND6
  //g1->FixParameter(1,sigma);
  g1->SetParNames("Channel for KE","Sigma","Norm","Background");
  h1->Draw();
  

  h1->Fit("convolution2","R","",3400,5500);// For Cs 

  //g1->ReleaseParameter(3); 
  //h1->Fit("convolution1","R","",1420,2500);// For Cs 

  // These are for Co -------------
  /*TF1 *g1 = new TF1("convolution1",KlenGauss,0,7000,4);// For Co spectra 
  double norm = h1->GetBinContent(520);// For Co  
  g1->SetParameters(channelKE,sigma,norm,0.0);
  g1->SetParNames("Channel for KE","Sigma","Norm","Background");
  
  h1->Fit("convolution1","R","",4500,6300);// For Co 
  TF1 *g3 = new TF1("klen2",Klen,0,7000,4);
  g3->SetParameters(channelKE,sigma1,norm,0.0);
  */
  //h1->SetLineColor(kRed);
  
  g3->SetLineColor(kBlue);
  g4->SetLineColor(kGreen);
  printf("%f and klein %f ,normalize %f\n",h1->GetBinContent(340),g2->Eval(h1->GetBinCenter(340)),normalize);
  
  g3->Draw("same");
  g4->Draw("same");
  h1->Draw("same");
  f1->Draw("same");
}
