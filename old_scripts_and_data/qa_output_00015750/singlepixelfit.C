#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <sstream>
#include <TFile.h>
#include <assert.h>

// I like to have include directives @Tim

Double_t fitf(Double_t *x, Double_t *par)
{
  Double_t arg = 0;

  Double_t err_func = 1+TMath::Erf(x[0]-par[2]);

  Double_t landau = par[0]*TMath::Landau(x[0],par[1],par[2]);

  Double_t gaussians = 0;
  gaussians += par[5] * TMath::Exp(-TMath::Power(x[0]-par[4],2)/(2*par[6]*par[6]));

  for (int i = 1;i<4;i++)
    {
      // par[2*i+6] = par[6];
      gaussians += par[2*i+5] * TMath::Exp(-TMath::Power(x[0]-(i*par[3]+par[4]),2)/(2*par[2*i+6]*par[2*i+6]));
    }
  //par[0] = landau amplitude
  //par[1] = landau mpv
  //par[2] = landau sigma
  //par[3] = gap spacing
  //par[4] = first peak mean
  //par[5] = first peak amplitude
  //par[6] = first peak sigma
  //par[7,9,11,13] = other peak amplitudes
  //par[8,10,12,14] = other peak sigmas

  Double_t fitval = landau+gaussians;

  return fitval;
}



void singlepixelfit(int chnlnum)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);   //make the plot list all the fit information


  TFile * infile  = TFile::Open("histograms.root");
  TH1D* h_alladc[384];

  for (int i = 0; i < 384;i++)
    {
      std::stringstream aa;
      aa <<"h_alladc_"<<i;
      h_alladc[i] = (TH1D*)infile->Get(aa.str().c_str());
    }
  

  TF1* f_singlepixels[384];

  {
    int i = chnlnum;
  
    std::stringstream aa;
    aa <<"bob_"<<i;
    f_singlepixels[i] = new TF1(aa.str().c_str(),fitf,0,160,13);
    




      TF1* f_tmp = new TF1("f_tmp","gaus",20,40);
      h_alladc[i]->Fit(f_tmp,"Q0","",20,40);

      TF1* f_landautmp = new TF1("f_landautmp","landau",0,15);
      h_alladc[i]->Fit(f_landautmp,"Q0","",1.5,18);


      TF1* f_landaugaus = new TF1("f_landaugaus","landau(0) + gaus(3)",0,30);
      f_landaugaus->SetParameters(f_landautmp->GetParameter(0),f_landautmp->GetParameter(1),f_landautmp->GetParameter(2),f_tmp->GetParameter(0),f_tmp->GetParameter(1),f_tmp->GetParameter(2));
      h_alladc[i]->Fit(f_landaugaus,"Q0","",0.5,45);


      // f_singlepixels[i]->SetParameters(f_tmp->GetParameter(0),30,5,f_tmp->GetParameter(0)/3,28,f_tmp->GetParameter(0)/10, f_tmp->GetParameter(0)/20,f_tmp->GetParameter(0)/100);
      
      f_singlepixels[i]->FixParameter(0,f_landaugaus->GetParameter(0));
      f_singlepixels[i]->FixParameter(1,f_landaugaus->GetParameter(1));
      f_singlepixels[i]->FixParameter(2,f_landaugaus->GetParameter(2));
      f_singlepixels[i]->FixParameter(4,f_landaugaus->GetParameter(4));
      f_singlepixels[i]->FixParameter(5,f_landaugaus->GetParameter(3));
      f_singlepixels[i]->FixParameter(6,f_landaugaus->GetParameter(5));


      f_singlepixels[i]->SetParameter(3,28);
      // f_singlepixels[i]->SetParameter(4,f_tmp->GetParameter(1));
      // f_singlepixels[i]->SetParameter(5,f_tmp->GetParameter(0));
      // f_singlepixels[i]->SetParameter(6,f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParameter(7,f_tmp->GetParameter(0)/2);
      f_singlepixels[i]->SetParameter(8,f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParameter(9,f_tmp->GetParameter(0)/5);
      f_singlepixels[i]->SetParameter(10,f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParameter(11,f_tmp->GetParameter(0)/10);
      f_singlepixels[i]->SetParameter(12,f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParameter(13,f_tmp->GetParameter(0)/15);
      f_singlepixels[i]->SetParameter(14,f_tmp->GetParameter(2));

      
      // f_singlepixels[i]->SetParLimits(1,0,15);
      // f_singlepixels[i]->SetParLimits(2,0,100);
      f_singlepixels[i]->SetParLimits(3, 20,40);
      // f_singlepixels[i]->SetParLimits(4,10,40);
      // f_singlepixels[i]->SetParLimits(5, f_landautmp->Eval(f_tmp->GetParameter(1))/2,1000000);
      // f_singlepixels[i]->SetParLimits(6, 0,100);
      f_singlepixels[i]->SetParLimits(7, f_landautmp->Eval(35+f_tmp->GetParameter(1)),1000000);
      f_singlepixels[i]->SetParLimits(9, f_landautmp->Eval(70+f_tmp->GetParameter(1)),1000000);
      f_singlepixels[i]->SetParLimits(11, f_landautmp->Eval(105+f_tmp->GetParameter(1)),1000000);
      f_singlepixels[i]->SetParLimits(13, f_landautmp->Eval(140+f_tmp->GetParameter(1)),1000000);


      f_singlepixels[i]->SetParLimits(8,0.5*f_tmp->GetParameter(2),1.15*f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParLimits(10,0.5*f_tmp->GetParameter(2),1.15*f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParLimits(12,0.5*f_tmp->GetParameter(2),1.15*f_tmp->GetParameter(2));
      f_singlepixels[i]->SetParLimits(14,0.5*f_tmp->GetParameter(2),1.15*f_tmp->GetParameter(2));




  //par[0] = landau amplitude
  //par[1] = landau mpv
  //par[2] = landau sigma
  //par[3] = gap spacing
  //par[4] = first peak mean
  //par[5] = first peak amplitude
  //par[6] = first peak sigma
  //par[7,9,11,13] = other peak amplitudes
  //par[8,10,12,14] = other peak sigmas






      h_alladc[i]->Fit(f_singlepixels[i],"Q","",1.5,140);

  }


TCanvas*c1 = new TCanvas("c1","",700,500);
 gPad->SetLogy();
h_alladc[chnlnum]->Draw();
}
