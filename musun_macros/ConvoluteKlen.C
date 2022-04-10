/*
Step 1: Just created a function that convolutes Klein Nishina cross section with a gaussian
Step 2: Create a fit function with three parameters - max photon energy (required for calibrating with channel),width of the gaussian and normalization
THIS IS THE RIGHT ONE TO BE USED...
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
// 137Cs 667 keV energy of CE (max photon energy) which has max KE as 482 keV

//TF1 *g1 = new TF1("g1",);
//*g2,*g3,*g4,*g5;

//TF1 *g2;
//double par[3];
//double err[3];
//double err1[3];
/*
double crossSection(double energyK, double eMax)
{
  //The energies must be in keV
  
  if (energyK>2*eMax*eMax/(511.0+2*eMax)) return 0; 
  double cos = 1-(energyK*511.0)/(eMax*(eMax-energyK));
  double P = 1/(1+(eMax/511.0)*(1-cos)); //ratio of scattered to incident photon energy
  double KlenCross = P*P*(P+(1/P)-1+cos*cos);
  //printf("Cross section is %f \n",KlenCross);
  // h4->Fill(180*TMath::ACos(cos)/3.14);
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

void ConvoluteKlen()
{
  double sum = 0.0, sumgauss = 0.0;
  double cs;
  TH1F *h1 = new TH1F("h1","h1",600,0,600);

  
  TH1F *h2 = new TH1F("h2","h2",600,0,600);
  TH1F *h3 = new TH1F("h3","h3",700,0,700);


  
  for (int i=0; i<481;i++)
    {
      //if (i<20 || i>500)h1->SetBinContent(i,0);
      h1->SetBinContent(i,crossSection((double)i,667.0));
      }

  h1->Draw();  

  for (int i=0; i<481;i++)
   {
     h2->SetBinContent(i,TMath::Gaus(i,20,20));
     //Gaus function takes Parameters - 1st plotting pt. 2nd mean (peak)& 3rd sigma
      //h2->SetBinContent(i,2*TMath::Exp(-0.5*i)); 
   }
  
  //h2->SetLineColor(kRed);
  //h2->Draw("same");
  
  int sigma= 20;

  TF1 *f1 = new TF1("f1","TMath::Gaus(x,[0],[1])",-3*sigma,+3*sigma);
  f1->SetParameters(0,sigma);

  for (int j=0; j<520;j++)
   {
     sum=0.0;
     sumgauss = 0.0;
     //int x = j-3*sigma;
     //printf("x = %d\t",x);
     for(int i=0;i<520;i++)
       {
	 int k = -i+j;
	 if (k<-3*sigma) continue;
	 if (k>+3*sigma) continue;
	   {
	     
	     double gauss = f1->Eval(k);
	     sumgauss += gauss;//This is done to normalize the Gaussian for parts which were not normalized
	     sum += gauss*crossSection((double)i,667.0);
	     //sum += TMath::Gaus(j,0,sigma)*crossSection((double)i,667.0);
	   }
	   // printf("i, j, k, gaus, cs, sum = %d %d %d %f %f %f\n",i,j,k,gauss,crossSection((double)i,667.0),sum);
       } 

     //if (j<10)  printf("j,sum %d %f\n",j,sum);
     h3->SetBinContent(j,sum/sumgauss);
     
   }

  //return sum;
  h3->Scale(h1->Integral(0,520)/h3->Integral(0,520),"");
  h3->SetLineColor(kBlue);
  h3->Draw("same");
 }
