#include <TROOT.h>
#include <TSystem.h>
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
#include <TLine.h>

#include <cmath>
#include <stdexcept>
#include <sys/stat.h>

// divides to errorful arrays of numbers x + dx / y + dy and puts result in z, dz
void divide_err(size_t n, double *x, double *dx, double *y, double *dy, double *z, double *dz) {
  for (int i = 0; i < n; i++) {
    double x0 = x[i];
    double dx0 = dx[i];
    double y0 = y[i];
    double dy0 = dy[i];
    double z0 = x0/y0;
    double dz0 = z0*pow(pow(dx0/x0, 2.0) + pow(dy0/y0, 2.0), 0.5);
    z[i] = z0;
    dz[i] = dz0;
  }
}

void sp_fit_get_peaks(TF1* sp_fit, double p1_0, double p2_0, double *p1, double *p2) {
  *p1 = sp_fit->GetMaximumX(p1_0 - 5, p1_0 + 5);
  *p2 = sp_fit->GetMaximumX(p2_0 - 5, p2_0 + 5);
}

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

int all_fits (int run_num, char *path_to_runs, TFile *histograms) {
  int n_channels = 0;
  
  char *file_prefix = strdup(Form("%s/%i", path_to_runs, run_num));
  //printf("FILE_PREFIX IS [%s]\n", file_prefix);
  struct stat file_stat;
  if (stat(file_prefix, &file_stat) == -1) {
    //printf("TRIED TO MKDIR\n");
    mkdir(file_prefix, 0700);
  } else {
    //printf("DID NOT TRY TO MKDIR\n");
  }
  if (stat(Form("%s/fit_channel_hists", file_prefix), &file_stat) == -1) {
    //printf("TRIED TO MKDIR\n");
    mkdir(Form("%s/fit_channel_hists", file_prefix), 0700);
  } else {
    //printf("DID NOT TRY TO MKDIR\n");
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);   //make the plot list all the fit information
  
  // retrieve all histograms from the root file
  TH1D *h_alladc[384];
  for (int i = 0; i < 384;i++) {
    h_alladc[i] = (TH1D*) histograms->Get(Form("h_alladc_%i", i));
  }
  
  TF1 *f_singlepixels[384];
  Double_t sp_gaps[384];
  Double_t sp_gaps_err[384];
  Double_t chisqr_ndfs[384];
  Double_t first_gauss_mean[384];
  Double_t first_gauss_mean_err[384];
  Double_t first_gauss_amplitude[384];
  Double_t first_gauss_amplitude_err[384];
  Double_t first_gauss_sigma[384];
  Double_t first_gauss_sigma_err[384];
  Double_t other_gauss_amplitudes[384][4];
  Double_t other_gauss_amplitudes_err[384][4];
  Double_t other_gauss_sigmas[384][4];
  Double_t other_gauss_sigmas_err[384][4];
  Double_t landau_amplitudes[384];
  Double_t landau_amplitudes_err[384];
  Double_t landau_mpvs[384];
  Double_t landau_mpvs_err[384];
  Double_t landau_sigmas[384];
  Double_t landau_sigmas_err[384];

  Double_t actual_peak1_height[384];
  Double_t actual_peak2_height[384];

  for (int i = 0; i < 384; i++) {
    if (h_alladc[i]->GetEntries() < 1e2) {
      printf("rejected channel %i for too few entries\n", i);
      continue;
    }
    
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
    if (sp_fit->IsEmpty()) {
      printf("rejected channel %i for null or empty sp_fit\n", i);
      continue;
    }
    Double_t sp_gap = sp_fit->Parameter(3);
    sp_gaps[i] = sp_gap;
    sp_gaps_err[i] = sp_fit->ParError(3);
    chisqr_ndfs[i] = sp_fit->Chi2() / sp_fit->Ndf();
    first_gauss_mean[i] = sp_fit->Parameter(4);
    first_gauss_mean_err[i] = sp_fit->ParError(4);
    first_gauss_amplitude[i] = sp_fit->Parameter(5);
    first_gauss_amplitude_err[i] = sp_fit->ParError(5);
    first_gauss_sigma[i] = sp_fit->Parameter(6);
    first_gauss_sigma_err[i] = sp_fit->ParError(6);
    for (int j = 0; j < 4; j++) {
      other_gauss_amplitudes[i][j] = sp_fit->Parameter(2*j + 7);
      other_gauss_amplitudes_err[i][j] = sp_fit->ParError(2*j + 7);
      other_gauss_sigmas[i][j] = sp_fit->Parameter(2*j + 8);
      other_gauss_sigmas_err[i][j] = sp_fit->ParError(2*j + 8);
    }
    landau_amplitudes[i] = sp_fit->Parameter(0);
    landau_amplitudes_err[i] = sp_fit->ParError(0);
    landau_mpvs[i] = sp_fit->Parameter(1);
    landau_mpvs_err[i] = sp_fit->ParError(1);
    landau_sigmas[i] = sp_fit->Parameter(2);
    landau_sigmas_err[i] = sp_fit->ParError(2);

    TF1 *landau = new TF1("landau", "landau", 0.0, 200.0);
    landau->SetParameter(0, landau_amplitudes[i]);
    landau->SetParameter(1, landau_mpvs[i]);
    landau->SetParameter(2, landau_sigmas[i]);
    landau->SetLineStyle(1); landau->SetLineColor(kGreen); landau->SetLineWidth(1);
    
    //TLine* landau_ampl = new TLine(0.0, landau_amplitudes[i], 200.0, landau_amplitudes[i]);
    //landau_ampl->SetLineColor(kGreen);
    //landau_ampl->SetLineStyle(2);
    //landau_ampl->SetLineWidth(1);
    

    TF1 *gaus1 = new TF1("gaus1", "gaus", 0.0, 200.0);
    gaus1->SetParameter(0, first_gauss_amplitude[i]);
    gaus1->SetParameter(1, first_gauss_mean[i]);
    gaus1->SetParameter(2, first_gauss_sigma[i]);
    gaus1->SetLineStyle(1); gaus1->SetLineColor(kBlue); gaus1->SetLineWidth(1);

    //TLine* gauss1_ampl = new TLine(0.0, first_gauss_amplitude[i], 200.0, first_gauss_amplitude[i]);
    //gauss1_ampl->SetLineColor(kBlue);
    //gauss1_ampl->SetLineStyle(2);
    //gauss1_ampl->SetLineWidth(1);

    double p1, p2;
    sp_fit_get_peaks(f_singlepixels[i], landau->GetMaximumX(), first_gauss_mean[i], &p1, &p2);
    assert(0.0 < p1 && p1 < 50.0);
    assert(0.0 < p2 && p2 < 50.0);
    actual_peak1_height[i] = f_singlepixels[i]->Eval(p1);
    actual_peak2_height[i] = f_singlepixels[i]->Eval(p2);

    TF1 *gaus2 = new TF1("gaus2", "gaus", 0.0, 200.0);
    gaus2->SetParameter(0, other_gauss_amplitudes[i][0]);
    gaus2->SetParameter(1, first_gauss_mean[i] + sp_gap);
    gaus2->SetParameter(2, other_gauss_sigmas[i][0]);
    gaus2->SetLineStyle(1); gaus2->SetLineColor(kBlue); gaus2->SetLineWidth(1);

    TF1 *gaus3 = new TF1("gaus3", "gaus", 0.0, 200.0);
    gaus3->SetParameter(0, other_gauss_amplitudes[i][1]);
    gaus3->SetParameter(1, first_gauss_mean[i] + 2*sp_gap);
    gaus3->SetParameter(2, other_gauss_sigmas[i][1]);
    gaus3->SetLineStyle(1); gaus3->SetLineColor(kBlue); gaus3->SetLineWidth(1);

    TF1 *gaus4 = new TF1("gaus4", "gaus", 0.0, 200.0);
    gaus4->SetParameter(0, other_gauss_amplitudes[i][2]);
    gaus4->SetParameter(1, first_gauss_mean[i] + 3*sp_gap);
    gaus4->SetParameter(2, other_gauss_sigmas[i][2]);
    gaus4->SetLineStyle(1); gaus4->SetLineColor(kBlue); gaus4->SetLineWidth(1);
    
    TF1 *gaus5 = new TF1("gaus5", "gaus", 0.0, 200.0);
    gaus5->SetParameter(0, other_gauss_amplitudes[i][3]);
    gaus5->SetParameter(1, first_gauss_mean[i] + 4*sp_gap);
    gaus5->SetParameter(2, other_gauss_sigmas[i][3]);
    gaus5->SetLineStyle(1); gaus5->SetLineColor(kBlue); gaus5->SetLineWidth(1);


    
    TCanvas* c1 = new TCanvas(Form("c%i", i), "", 700, 500);
    gPad->SetLogy();
    h_alladc[i]->Draw();

    landau->Draw("SAME");
    //landau_ampl->Draw("SAME");
    gaus1->Draw("SAME");
    //gauss1_ampl->Draw("SAME");
    gaus2->Draw("SAME");
    gaus3->Draw("SAME");
    gaus4->Draw("SAME");
    gaus5->Draw("SAME");

    /*
    Double_t yndc = 1e5;
    TLine* p1_line = new TLine(p1, 0.0, p1, yndc);
    p1_line->SetLineColor(kRed);
    p1_line->SetLineStyle(2);
    p1_line->SetLineWidth(2);
    p1_line->Draw("SAME");
    TLine* p2_line = new TLine(p2, 0.0, p2, yndc);
    p2_line->SetLineColor(kRed);
    p2_line->SetLineStyle(2);
    p2_line->SetLineWidth(2);
    p2_line->Draw("SAME");
    */

    TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
    legend->AddEntry(f_landaugaus, Form("SP Gap: %.3f", sp_gap), "l");
    legend->AddEntry("", Form("ChiSqr/NDF: %.3f", chisqr_ndfs[i]));
    legend->Draw();

    c1->SaveAs(Form("%s/fit_channel_hists/channel_%i.png", file_prefix, i));
    n_channels++;
  }
  

  Double_t min_chisqr_ndf = chisqr_ndfs[0];
  int min_chisqr_idx = 0;
  Double_t max_chisqr_ndf = chisqr_ndfs[0];
  int max_chisqr_idx = 0;
  // now save all single pixel gaps to csv file
  FILE *fp = fopen(Form("%s/all_gaps.csv", file_prefix), "w+");
  fprintf(fp, "channel_num, sp_gap, chisqr/ndf, gauss1 mean, gauss1 ampl, gauss1 sigma, gauss2 ampl, gauss2 sigma, gauss3 ampl, gauss3 sigma, gauss4 ampl, gauss4 sigma, gauss5 ampl, gauss5 sigma, landau ampl, landau mpv, landau sigma, actual peak1 height, actual peak2 height"); // header
  for (int i = 0; i < 384; i++) {
    Double_t chisqr_ndf = chisqr_ndfs[i];
    fprintf(fp, "\n%i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f", i, sp_gaps[i], chisqr_ndf, first_gauss_mean[i], first_gauss_amplitude[i], first_gauss_sigma[i],
      other_gauss_amplitudes[i][0], other_gauss_sigmas[i][0], other_gauss_amplitudes[i][2], other_gauss_sigmas[i][2], other_gauss_amplitudes[i][2], other_gauss_sigmas[i][2], other_gauss_amplitudes[i][3], other_gauss_sigmas[i][3],
      landau_amplitudes[i], landau_mpvs[i], landau_sigmas[i], actual_peak1_height[i], actual_peak2_height[i]);
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
  c_chi->SaveAs(Form("%s/chisqr_by_channel.png", file_prefix));

  // compute the ratio of the 1st and 2nd peaks
  Double_t ratios[384];
  for (int i = 0; i < 384; i++) {
    ratios[i] = actual_peak1_height[i] / actual_peak2_height[i];
  }

  // compute average ratio by IB
  Double_t avg_ib_ratio[6];
  for (int ib = 0; ib < 6; ib++) {
    Double_t avg_ratio = 0.0;
    for (int i = ib*64; i < (ib + 1)*64; i++) {
      avg_ratio += ratios[i];
    }
    avg_ratio /= 64;
    avg_ib_ratio[ib] = avg_ratio;
  }
  Double_t ib_bins[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

  // plot ratio of 1st and 2nd peak as a fn of IB
  TCanvas *c_ib_ratio = new TCanvas("c_ib_ratio", "", 700, 500);
  TGraph *ib_ratio_graph = new TGraph(6, ib_bins, avg_ib_ratio);
  ib_ratio_graph->SetTitle("Peak Height Ratio by IB");
  ib_ratio_graph->GetXaxis()->SetTitle("IB");
  ib_ratio_graph->GetYaxis()->SetTitle("1st peak / 2nd peak");
  ib_ratio_graph->SetMarkerStyle(33);
  ib_ratio_graph->SetMarkerSize(2.0);
  ib_ratio_graph->GetXaxis()->SetLimits(0.0, 5.0);
  ib_ratio_graph->GetXaxis()->SetNdivisions(-5);
  ib_ratio_graph->Draw("AP");
  // canvas->Update();
  c_ib_ratio->SetGrid();
  c_ib_ratio->SaveAs(Form("%s/peak_ratio_by_ib.png", file_prefix));

  // plot ratio of 1st and 2nd peak as a fn of chisqr
  TCanvas *c_ratio_chi = new TCanvas("c_ratio_chi", "", 700, 500);
  TGraph *ratio_chi_graph = new TGraph(384, chisqr_ndfs, ratios);
  ratio_chi_graph->SetTitle("Peak Height Ratio by Chi Sqr");
  ratio_chi_graph->GetXaxis()->SetTitle("ChiSqr/NDF");
  ratio_chi_graph->GetYaxis()->SetTitle("1st peak / 2nd peak");
  ratio_chi_graph->SetMarkerStyle(33);
  ratio_chi_graph->SetMarkerSize(1.0);
  //ratio_chi_graph->GetXaxis()->SetLimits(0.0, 383.0);
  ratio_chi_graph->Draw("AP");
  // canvas->Update();
  c_ratio_chi->SetGrid();
  c_ratio_chi->SaveAs(Form("%s/peak_ratio_by_chisqr.png", file_prefix));



  //printf("max chisqr_ndf: %f (channel %i); min chisqr_ndf: %f (channel %i)\n", max_chisqr_ndf, max_chisqr_idx, min_chisqr_ndf, min_chisqr_idx);

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
  //printf("IBs 0-2: avg sp gap = %f; IBs 3-5 avg sp gap = %f\n", ib0_2_avg_gap, ib3_5_avg_gap);
  free(file_prefix);
  return n_channels;
}

void fit_all_runs() {
  const char *inDir = "./all_runs_test";
  
  char *dir = gSystem->ExpandPathName(inDir);
  void *dirp = gSystem->OpenDirectory(dir);

  char actualpath [PATH_MAX+1];
  char *path_to_fits = realpath("./all_run_fits", actualpath);

  const char *entry;
  const char *filename;
  Int_t n = 0;
  char const *start = "qa_output_000";
  TFile *hist_file;

  FILE *fp = fopen("./all_runs_stats.csv", "w+");
  fprintf(fp, "run_num, n_channels"); // header

  while ((entry = (char*) gSystem->GetDirEntry(dirp))) {
    if (!strncmp(entry, start, strlen(start))) {
      filename = gSystem->ConcatFileName(dir, entry);
      int run_num;
      struct stat statbuf;
      if (sscanf(entry, "qa_output_000%i", &run_num) == 1 && stat(filename, &statbuf) == 0 && S_ISDIR(statbuf.st_mode) && !strchr(entry, '.')
          && (hist_file = TFile::Open(Form("%s%s", filename, "/histograms.root")))) {
        printf("found run num %i at %s\n", run_num, filename);
        int n_channels = all_fits(run_num, path_to_fits, hist_file);
        fprintf(fp, "\n%i, %i", run_num, n_channels);
      }

    }
  }

  fclose(fp);
}