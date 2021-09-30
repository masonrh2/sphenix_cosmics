#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <sstream>
#include <TFile.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TGraph.h>
#include <assert.h>

/**
 * The function used for fitting
 */
Double_t fitf (Double_t *x, Double_t *par) {
  Double_t arg = 0;
  Double_t err_func = 1+TMath::Erf(x[0]-par[2]);
  Double_t landau = par[0]*TMath::Landau(x[0],par[1],par[2]);
  Double_t gaussians = 0;
  gaussians += par[5] * TMath::Exp(-TMath::Power(x[0]-par[4],2)/(2*par[6]*par[6]));
  for (int i = 1; i < 4; i++) {
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

  Double_t fitval = landau + gaussians;
  return fitval;
}

void all_fits () {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);   //make the plot list all the fit information

  TFile *infile = TFile::Open("histograms.root");
  
  // retrieve all histograms from the root file
  TH1D *h_alladc[384];
  for (int i = 0; i < 384;i++) {
    h_alladc[i] = (TH1D*) infile->Get(Form("h_alladc_%i", i));
  }
  
  TF1 *f_singlepixels[384];
  Double_t sp_gaps[384]; 
  Double_t chisqr_ndfs[384]; 

  for (int i = 0; i < 384; i++) {
    // this scheme could likely be proved upon by finding more suitable parameters for each histogram
    // rather than using the same constants for each histogram and sector...
    
    f_singlepixels[i] = new TF1(Form("f_%i", i), fitf,0,160,13);
    
    TF1* f_tmp = new TF1("f_tmp","gaus",20,40);
    h_alladc[i]->Fit(f_tmp,"Q0","",20,40);

    TF1* f_landautmp = new TF1("f_landautmp","landau",0,15);
    h_alladc[i]->Fit(f_landautmp,"Q0","",1.5,18);

    TF1* f_landaugaus = new TF1("f_landaugaus","landau(0) + gaus(3)",0,30);
    f_landaugaus->SetParameters(f_landautmp->GetParameter(0), f_landautmp->GetParameter(1) ,f_landautmp->GetParameter(2) ,f_tmp->GetParameter(0) ,f_tmp->GetParameter(1) ,f_tmp->GetParameter(2));
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

    TFitResultPtr sp_fit = h_alladc[i]->Fit(f_singlepixels[i],"Q S","",1.5,140);
    Double_t sp_gap = sp_fit->Parameter(3);
    sp_gaps[i] = sp_gap;
    chisqr_ndfs[i] = sp_fit->Chi2() / sp_fit->Ndf();
    //std::cout << Form("Channel %i, Parameter 3: %f", i, sp_gap) << std::endl;
    
    TCanvas* c1 = new TCanvas(Form("c%i", i), "", 700, 500);
    gPad->SetLogy();
    h_alladc[i]->Draw();

    TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
    legend->AddEntry(f_landaugaus, Form("SP Gap: %.3f", sp_gap), "l");
    legend->AddEntry("", Form("ChiSqr/NDF: %.3f", chisqr_ndfs[i]));
    legend->Draw();

    c1->SaveAs(Form("./fit_channel_hists/channel_%i.png", i));
  }
  

  Double_t min_chisqr_ndf = chisqr_ndfs[0];
  int min_chisqr_idx = 0;
  Double_t max_chisqr_ndf = chisqr_ndfs[0];
  int max_chisqr_idx = 0;
  // now save all single pixel gaps to csv file
  FILE *fp = fopen("all_gaps.csv", "w+");
  fprintf(fp, "channel_num, sp_gap, chisqr/ndf"); // header
  for (int i = 0; i < 384; i++) {
    Double_t chisqr_ndf = chisqr_ndfs[i];
    fprintf(fp, "\n%i, %f, %f", i, sp_gaps[i], chisqr_ndf);
    if (chisqr_ndf > max_chisqr_ndf) {
      max_chisqr_idx = i;
      max_chisqr_ndf = chisqr_ndf;
    }
    if (chisqr_ndf < min_chisqr_ndf) {
      min_chisqr_idx = i;
      min_chisqr_ndf = chisqr_ndf;
    }
  }
  fclose(fp);

  // plot chisqr / ndf by channel
  Double_t bins[384];
  for (int i = 0; i < 384; i++) {
    bins[i] = i;
  }
  TCanvas *c_chi = new TCanvas("c_chi", "", 700, 500);
  TGraph *chi_graph = new TGraph(384, bins, chisqr_ndfs);
  chi_graph->SetTitle("Single Pixel Fit Chi Sqr by Channel");
  chi_graph->GetXaxis()->SetTitle("Channel");
  chi_graph->GetYaxis()->SetTitle("ChiSqr/NDF");
  chi_graph->SetMarkerStyle(33);
  chi_graph->SetMarkerSize(1.0);
  chi_graph->GetXaxis()->SetLimits(0.0, 383.0);
  chi_graph->Draw("AP");
  // canvas->Update();
  c_chi->SaveAs("./chisqr_by_channel.png");

  printf("max chisqr_ndf: %f (channel %i); min chisqr_ndf: %f (channel %i)\n", max_chisqr_ndf, max_chisqr_idx, min_chisqr_ndf, min_chisqr_idx);

  // compute averages for IBs 0-2 and 3-5
  Double_t avg_gap = 0;
  for (int i = 0; i < 192; i++) {
    avg_gap += sp_gaps[i];
  }
  Double_t ib0_2_avg_gap = avg_gap / 192;
  avg_gap = 0;
  for (int i = 192; i < 384; i++) {
    avg_gap += sp_gaps[i];
  }
  Double_t ib3_5_avg_gap = avg_gap / 192;
  printf("IBs 0-2: avg sp gap = %f; IBs 3-5 avg sp gap = %f\n", ib0_2_avg_gap, ib3_5_avg_gap);
}