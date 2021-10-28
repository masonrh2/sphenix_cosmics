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
#include <TGraphErrors.h>

#include <cmath>
#include <stdexcept>
#include <sys/stat.h>
#include <pthread.h>

#define RUN_NUM 16945
#define SAVE_PLOTS true
#define BOUND_TOL 0.01

// TODO: multithreading

void sp_fit_get_peaks(TF1* sp_fit, double p1_0, double p2_0, double *p1, double *p2) {
  *p1 = sp_fit->GetMaximumX(p1_0 - 5, p1_0 + 5);
  *p2 = sp_fit->GetMaximumX(p2_0 - 5, p2_0 + 5);
}

/**
 * The function used for fitting
 */
Double_t fitf (Double_t *x, Double_t *par) {
  Double_t gaussians = 0;
  gaussians += par[2] * TMath::Exp(-TMath::Power(x[0]-par[3],2)/(2*par[4]*par[4]));
  gaussians += par[5] * TMath::Exp(-TMath::Power(x[0]-(par[0]+par[3]),2)/(2*par[6]*par[6]));
  for (int i = 0; i < 4; i++) {
    gaussians += par[2*i+7] * TMath::Exp(-TMath::Power(x[0]-((i+1)*par[1]+par[0]+par[3]),2)/(2*par[2*i+8]*par[2*i+8]));
  }
  // par[0] = first gap spacing
  // par[1] = other gap spacings
  // par[2] = first peak amplitude
  // par[3] = first peak mean
  // par[4] = first peak sigma
  // par[5,7,9,11,13] = other peak amplitudes
  // par[6,8,10,12,14] = other peak sigmas

  return gaussians;
}

Double_t first_gap_initial_guess = 22.0;
Double_t first_gap_bounds[2] = {18.0, 32.0};
Double_t other_gap_initial_guess = 30.0;
Double_t other_gap_bounds[2] = {20.0, 36.0};
Double_t first_two_peaks_initial_guess = 1e4;
Double_t first_two_peaks_bounds[2] = {1e3, 1e5};
Double_t other_peaks_initial_guess = 1e3;
Double_t other_peaks_bounds[2] = {10, 1e5};
Double_t sigma_initial_guess = 8.0;
Double_t sigma_bounds[2] = {2.0, 20.0};

void fit_no_pedestal(int run_num) {
  ROOT::EnableImplicitMT();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);   //make the plot list all the fit information
  
  TFile *hist_file = TFile::Open(Form("./all_runs/qa_output_nopedestal_000%i/histograms.root", run_num));
  if (!hist_file) {
    throw std::runtime_error("unable to open histogram file");
  }

  char *file_prefix = strdup(Form("./no_pedestal_hists_%i", run_num));

  // retrieve all histograms from the root file
  TH1D *h_alladc[384];
  for (int i = 0; i < 384;i++) {
    h_alladc[i] = (TH1D*) hist_file->Get(Form("h_alladc_%i", i));
  }
  
  TF1 *f_singlepixels[384];

  Double_t first_gap[384][2];
  Double_t other_gaps[384][2];
  Double_t first_ampl[384][2];
  Double_t first_mean[384][2];
  Double_t first_sigma[384][2];
  Double_t other_ampl[384][5][2];
  Double_t other_sigma[384][5][2];


  bool fit_success[384] = {false};

  for (int i = 0; i < 384; i++) {
    bool success = true;
    
    Double_t first_above = h_alladc[i]->FindFirstBinAbove(1e3);
    Double_t first_peak_guess = first_above + 16;
    
    Double_t mean_guesses[6];
    mean_guesses[0] = first_peak_guess;
    mean_guesses[1] = mean_guesses[0] + 22;
    for (int j = 2; j < 6; j++) {
      mean_guesses[j] = mean_guesses[j - 1] + 28;
    }

    Double_t left = mean_guesses[0] - 40.0;
    Double_t right = mean_guesses[5] + 14;

    // par[0] = first gap spacing
    // par[1] = other gap spacings
    // par[2] = first peak amplitude
    // par[3] = first peak mean
    // par[4] = first peak sigma
    // par[5,7,9,11,13] = other peak amplitudes
    // par[6,8,10,12,14] = other peak sigmas

    f_singlepixels[i] = new TF1(Form("f_%i", i), fitf, left, right, 15);

    f_singlepixels[i]->SetParameter(0, first_gap_initial_guess);
    f_singlepixels[i]->SetParLimits(0, first_gap_bounds[0], first_gap_bounds[1]);
    f_singlepixels[i]->SetParameter(1, other_gap_initial_guess);
    f_singlepixels[i]->SetParLimits(1, other_gap_bounds[0], other_gap_bounds[1]);
    f_singlepixels[i]->SetParameter(2, first_two_peaks_initial_guess);
    f_singlepixels[i]->SetParLimits(2, first_two_peaks_bounds[0], first_two_peaks_bounds[1]);
    f_singlepixels[i]->SetParameter(3, mean_guesses[0]);
    f_singlepixels[i]->SetParLimits(3, mean_guesses[0] - 5, mean_guesses[0] + 5);
    f_singlepixels[i]->SetParameter(4, sigma_initial_guess);
    f_singlepixels[i]->SetParLimits(4, sigma_bounds[0], sigma_bounds[1]);
    
    f_singlepixels[i]->SetParameter(5, first_two_peaks_initial_guess);
    f_singlepixels[i]->SetParLimits(5, first_two_peaks_bounds[0], first_two_peaks_bounds[1]);
    f_singlepixels[i]->SetParameter(6, sigma_initial_guess);
    f_singlepixels[i]->SetParLimits(6, sigma_bounds[0], sigma_bounds[1]);
    for (int j = 1; j < 5; j++) {
      f_singlepixels[i]->SetParameter(2*j+5, other_peaks_initial_guess);
      f_singlepixels[i]->SetParLimits(2*j+5, other_peaks_bounds[0], other_peaks_bounds[1]);
      f_singlepixels[i]->SetParameter(2*j+6, sigma_initial_guess);
      f_singlepixels[i]->SetParLimits(2*j+6, sigma_bounds[0], sigma_bounds[1]);
    }

    f_singlepixels[i]->SetLineWidth(2);

    TFitResultPtr sp_fit = h_alladc[i]->Fit(f_singlepixels[i], "Q S L M", "", left, right);
    if (sp_fit->IsEmpty()) {
      printf("channel %i: fit result empty!\n", i);
      success = false;
    } else if (!sp_fit->IsValid()) {
      printf("channel %i: fit did not converge!\n", i);
      success = false;
    }
    fit_success[i] = success;

    first_gap[i][0] = sp_fit->Parameter(0);
    first_gap[i][1] = sp_fit->ParError(0);
    other_gaps[i][0] = sp_fit->Parameter(1);
    other_gaps[i][1] = sp_fit->ParError(1);
    first_ampl[i][0] = sp_fit->Parameter(2);
    first_ampl[i][1] = sp_fit->ParError(2);
    first_mean[i][0] = sp_fit->Parameter(3);
    first_mean[i][1] = sp_fit->ParError(3);
    first_sigma[i][0] = sp_fit->Parameter(4);
    first_sigma[i][1] = sp_fit->ParError(4);
    for (int j = 0; j < 5; j++) {
      other_ampl[i][j][0] = sp_fit->Parameter(2*j+5);
      other_ampl[i][j][1] = sp_fit->ParError(2*j+5);
      other_sigma[i][j][0] = sp_fit->Parameter(2*j+6);
      other_sigma[i][j][1] = sp_fit->ParError(2*j+6);
    }

    // gaps
    if (first_gap[i][0] - first_gap_bounds[0] < BOUND_TOL) {
      printf("***channel %i: first gap (%f) at lower bound!\n", i, first_gap[i][0]);
    } else if (first_gap_bounds[1] - first_gap[i][0] < BOUND_TOL) {
      printf("***channel %i: first gap (%f) at upper bound!\n", i, first_gap[i][0]);
    }
    if (other_gaps[i][0] - other_gap_bounds[0] < BOUND_TOL) {
      printf("***channel %i: other gaps (%f) at lower bound!\n", i, other_gaps[i][0]);
    } else if (other_gap_bounds[1] - other_gaps[i][0] < BOUND_TOL) {
      printf("***channel %i: other gaps (%f) at upper bound!\n", i, other_gaps[i][0]);
    }
    // first peak
    if (first_ampl[i][0] - first_two_peaks_bounds[0] < BOUND_TOL) {
      printf("***channel %i: first ampl (%f) at lower bound!\n", i, first_ampl[i][0]);
    } else if (first_two_peaks_bounds[1] - first_ampl[i][0] < BOUND_TOL) {
      printf("***channel %i: first ampl (%f) at upper bound!\n", i, first_ampl[i][0]);
    }
    if (first_mean[i][0] - left < BOUND_TOL) {
      printf("***channel %i: first mean (%f) at lower bound!\n", i, first_mean[i][0]);
    } else if (right - first_mean[i][0] < BOUND_TOL) {
      printf("***channel %i: first mean (%f) at upper bound!\n", i, first_mean[i][0]);
    }
    if (first_sigma[i][0] - sigma_bounds[0] < BOUND_TOL) {
      printf("***channel %i: first sigma (%f) at lower bound!\n", i, first_sigma[i][0]);
    } else if (sigma_bounds[1] - first_sigma[i][0] < BOUND_TOL) {
      printf("***channel %i: first sigma (%f) at upper bound!\n", i, first_sigma[i][0]);
    }
    // second peak
    if (other_ampl[i][0][0] - first_two_peaks_bounds[0] < BOUND_TOL) {
      printf("***channel %i: second ampl (%f) at lower bound!\n", i, other_ampl[i][0][0]);
    } else if (first_two_peaks_bounds[1] - other_ampl[i][0][0] < BOUND_TOL) {
      printf("***channel %i: second ampl (%f) at upper bound!\n", i, other_ampl[i][0][0]);
    }
    if (other_sigma[i][0][0] - sigma_bounds[0] < BOUND_TOL) {
      printf("***channel %i: second sigma (%f) at lower bound!\n", i, other_sigma[i][0][0]);
    } else if (sigma_bounds[1] - other_sigma[i][0][0] < BOUND_TOL) {
      printf("***channel %i: second sigma (%f) at upper bound!\n", i, other_sigma[i][0][0]);
    }
    // other peaks
    for (int j = 1; j < 5; j++) {
      if (other_ampl[i][j][0] - other_peaks_bounds[0] < BOUND_TOL) {
        printf("***channel %i: %ith ampl (%f) at lower bound!\n", i, j+2, other_ampl[i][j][0]);
      } else if (other_peaks_bounds[1] - other_ampl[i][j][0] < BOUND_TOL) {
        printf("***channel %i: %ith ampl (%f) at upper bound!\n", i, j+2, other_ampl[i][j][0]);
      }
      if (other_sigma[i][j][0] - sigma_bounds[0] < BOUND_TOL) {
        printf("***channel %i: %ith sigma (%f) at lower bound!\n", i, j+2, other_sigma[i][j][0]);
      } else if (sigma_bounds[1] - other_sigma[i][j][0] < BOUND_TOL) {
        //printf("***channel %i: %ith sigma (%f) at upper bound!\n", i, j+2, other_sigma[i][j][0]);
      }
    }
    
    Double_t gaus_ampl[6];
    Double_t gaus_means[6];
    Double_t gaus_sigmas[6];
    gaus_ampl[0] = first_ampl[i][0];
    gaus_means[0] = first_mean[i][0];
    gaus_sigmas[0] = first_sigma[i][0];
    gaus_ampl[1] = other_ampl[i][0][0];
    gaus_means[1] = gaus_means[0] + first_gap[i][0];
    gaus_sigmas[1] = other_sigma[i][0][0];
    for (int j = 1; j < 5; j++) {
      gaus_ampl[j + 1] = other_ampl[i][j][0];
      gaus_means[j + 1] = gaus_means[j] + other_gaps[i][0];
      gaus_sigmas[j + 1] = other_sigma[i][j][0];
    }

    if (SAVE_PLOTS) {
      h_alladc[i]->GetXaxis()->SetRangeUser(left, right);    
      //h_alladc[i]->GetXaxis()->SetRangeUser(1200.0, 2000.0);    

      TCanvas* c1 = new TCanvas(Form("c%i", i), "", 700, 500);

      h_alladc[i]->Draw();

      gPad->SetLogy();
      gPad->Update();

      
      /*
      for (int j = 0; j < 6; j++) {
        TLine* line = new TLine(mean_guesses[j], pow(10.0, c1->GetUymin()), mean_guesses[j], pow(10.0, c1->GetUymax()));
        line->SetLineColor(kBlue);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("SAME");
      }
      */

      Double_t fit_means[6];
      fit_means[0] = first_mean[i][0];
      fit_means[1] = fit_means[0] + first_gap[i][0];
      for (int j = 2; j < 6; j++) {
        fit_means[j] = fit_means[j - 1] + other_gaps[i][0];
      }
      for (int j = 0; j < 6; j++) {
        TLine* line = new TLine(fit_means[j], pow(10.0, c1->GetUymin()), fit_means[j], pow(10.0, c1->GetUymax()));
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("SAME");
      }
      
      
      for (int j = 0; j < 6; j++) {
        TF1 *gaus = new TF1("gaus", "gaus", left, right);
        gaus->SetParameter(0, gaus_ampl[j]);
        gaus->SetParameter(1, gaus_means[j]);
        gaus->SetParameter(2, gaus_sigmas[j]);
        gaus->SetLineStyle(2);
        gaus->SetLineColor(kRed);
        gaus->SetLineWidth(1);
        gaus->SetRange(gaus_means[j] - 20, gaus_means[j] + 20);
        gaus->Draw("SAME");
      }
      

      /*
      TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
      legend->AddEntry(f_landaugaus, Form("SP Gap: %.3f", sp_gap), "l");
      legend->AddEntry("", Form("ChiSqr/NDF: %.3f", chisqr_ndfs[i]));
      legend->Draw();
      */

      //c1->SetGridy();
      c1->SaveAs(Form("%s/channel_%i.png", file_prefix, i));
    }
  }
  
  
  FILE *outfile = fopen(Form("%s/no_pedestal_params.csv", file_prefix), "w+");
  fprintf(outfile, "channel, first gap spacing, first gap spacing err, other gap spacing, other gap spacing err, 1st peak amplitude, 1st peak amplitude err, 1st peak mean, 1st peak mean err, 1st peak sigma, 1st peak sigma err, 2nd peak amplitudes, 2nd peak amplitudes err, 2nd peak sigmas, 2nd peak sigmas err, 3rd peak amplitudes, 3rd peak amplitudes err, 3rd peak sigmas, 3rd peak sigmas err, 4th peak amplitudes, 4th peak amplitudes err, 4th peak sigmas, 4th peak sigmas err, 5th peak amplitudes, 5th peak amplitudes err, 5th peak sigmas, 5th peak sigmas err, 6th peak amplitudes, 6th peak amplitudes err, 6th peak sigmas, 6th peak sigmas err");
  for (short i = 0; i < 384; i++) {
    fprintf(outfile, "\n%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", i, first_gap[i][0], first_gap[i][1], other_gaps[i][0], other_gaps[i][1], first_ampl[i][0], first_ampl[i][1], first_mean[i][0], first_mean[i][1], first_sigma[i][0], first_sigma[i][1]);
    for (short j = 0; j < 5; j++) {
      fprintf(outfile, ",%f,%f,%f,%f", other_ampl[i][j][0], other_ampl[i][j][1], other_sigma[i][j][0], other_sigma[i][j][1]);
    }
  }
  
  fclose(outfile);

  // compute block gaps
  Double_t blocks[96];
  Double_t blocks_err[96];
  Double_t block_gaps[96];
  Double_t block_gaps_err[96];
  for (int i = 0; i < 96; i++) {
    blocks[i] = i + 0.5;
    blocks_err[i] = 0;
    block_gaps[i] = 0.0;
    block_gaps_err[i] = 0.0;
    for (int j = 0; j < 4; j++) {
      block_gaps[i] += other_gaps[4*i + j][0];
      block_gaps_err[i] += other_gaps[4*i + j][1];
    }
    block_gaps[i] /= 4;
    block_gaps_err[i] /= 4;
  }
  TGraphErrors *gap_graph = new TGraphErrors(96, blocks, block_gaps, blocks_err, block_gaps_err);
  gap_graph->SetTitle(Form("Run %i Single Pixel Gaps (No Pedestal)", run_num));
  gap_graph->GetXaxis()->SetTitle("Block");
  gap_graph->GetYaxis()->SetTitle("Single Pixel Gap");
  gap_graph->GetXaxis()->SetRangeUser(0.0, 96.0);
  gap_graph->SetMarkerStyle(21);
  gap_graph->SetMarkerColor(kBlue);
  gap_graph->GetYaxis()->SetRangeUser(20.0, 40.0);
  TCanvas *c2 = new TCanvas("c2", "", 700, 500);
  gap_graph->Draw("AP");
  gPad->Update();
  for (int j = 1; j < 6; j++) {
    TLine* line = new TLine(j*16, c2->GetUymin(), j*16, c2->GetUymax());
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");
  }
  c2->SetGridy();
  c2->SaveAs(Form("%s/gaps.png", file_prefix));

  // compute difference between first gap and other gaps
  Double_t gap_diffs[96];
  Double_t gap_diffs_err[96];
  for (int i = 0; i < 96; i++) {
    gap_diffs[i] = 0.0;
    gap_diffs_err[i] = 0.0;
    for (int j = 0; j < 4; j++) {
      gap_diffs[i] += (first_gap[4*i + j][0] - other_gaps[4*i + j][0]);
      gap_diffs_err[i] += first_gap[4*i + j][1] + other_gaps[4*i + j][1];
    }
    gap_diffs[i] /= 4;
    gap_diffs_err[i] /= 4;
  }
  TGraphErrors *gap_diff_graph = new TGraphErrors(96, blocks, gap_diffs, blocks_err, gap_diffs_err);
  gap_diff_graph->SetTitle(Form("Run %i First Gap - Other Gaps (No Pedestal)", run_num));
  gap_diff_graph->GetXaxis()->SetTitle("Block");
  gap_diff_graph->GetYaxis()->SetTitle("First Gap - Other Gaps");
  gap_diff_graph->GetXaxis()->SetRangeUser(0.0, 96.0);
  gap_diff_graph->SetMarkerStyle(21);
  gap_diff_graph->SetMarkerColor(kBlue);
  gap_diff_graph->GetYaxis()->SetRangeUser(-20.0, 20.0);
  TCanvas *c3 = new TCanvas("c3", "", 700, 500);
  gap_diff_graph->Draw("AP");
  gPad->Update();
  for (int j = 1; j < 6; j++) {
    TLine* line = new TLine(j*16, c3->GetUymin(), j*16, c3->GetUymax());
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");
  }
  c3->SetGridy();
  c3->SaveAs(Form("%s/gap_diffs.png", file_prefix));

  // compute ratio between 1st and 2nd peak
  Double_t peak_ratios[96];
  for (int i = 0; i < 96; i++) {
    peak_ratios[i] = 0.0;
    for (int j = 0; j < 4; j++) {
      peak_ratios[i] += first_ampl[4*i + j][0] / other_ampl[4*i + j][0][0];
    }
    peak_ratios[i] /= 4;
  }
  Double_t peak_ratios_err[96];
  for (int i = 0; i < 96; i++) {
    peak_ratios_err[i] = 0.0;
    for (int j = 0; j < 4; j++) {
      peak_ratios_err[i] += peak_ratios[i]*pow(pow(first_ampl[4*i + j][1]/first_ampl[4*i + j][0], 2.0) + pow(other_ampl[4*i + j][0][1]/other_ampl[4*i + j][0][0], 2.0), 0.5);
    }
    peak_ratios_err[i] /= 4;
  }
  TGraphErrors *peak_ratio_graph = new TGraphErrors(96, blocks, peak_ratios, blocks_err, peak_ratios_err);
  peak_ratio_graph->SetTitle(Form("Run %i 1st Peak / 2nd Peak (No Pedestal)", run_num));
  peak_ratio_graph->GetXaxis()->SetTitle("Block");
  peak_ratio_graph->GetYaxis()->SetTitle("1st Peak / 2nd Peak");
  peak_ratio_graph->GetXaxis()->SetRangeUser(0.0, 96.0);
  peak_ratio_graph->SetMarkerStyle(21);
  peak_ratio_graph->SetMarkerColor(kBlue);
  peak_ratio_graph->GetYaxis()->SetRangeUser(0.0, 3.0);
  TCanvas *c4 = new TCanvas("c4", "", 700, 500);
  peak_ratio_graph->Draw("AP");
  gPad->Update();
  for (int j = 1; j < 6; j++) {
    TLine* line = new TLine(j*16, c4->GetUymin(), j*16, c4->GetUymax());
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");
  }
  c4->SetGridy();
  c4->SaveAs(Form("%s/peak_ratios.png", file_prefix));

  // compute ratio between 1st and 2nd peak
  Double_t chi2_block[96];
  Double_t chi2_block_err[96];
  Double_t max_chi2 = 0.0;
  for (int i = 0; i < 96; i++) {
    chi2_block[i] = 0.0;
    chi2_block_err[i] = 0.0;
    for (int j = 0; j < 4; j++) {
      chi2_block[i] += f_singlepixels[4*i + j]->GetChisquare() / f_singlepixels[4*i + j]->GetNDF();
    }
    chi2_block[i] /= 4;
    if (chi2_block[i] > max_chi2) {
      max_chi2 = chi2_block[i];
    }
  }
  TGraphErrors *chi2_block_graph = new TGraphErrors(96, blocks, chi2_block, blocks_err, chi2_block_err);
  chi2_block_graph->SetTitle(Form("Run %i Chi2/NDF (No Pedestal)", run_num));
  chi2_block_graph->GetXaxis()->SetTitle("Block");
  chi2_block_graph->GetYaxis()->SetTitle("Chi2/NDF");
  chi2_block_graph->GetXaxis()->SetRangeUser(0.0, 96.0);
  chi2_block_graph->SetMarkerStyle(21);
  chi2_block_graph->SetMarkerColor(kBlue);
  chi2_block_graph->GetYaxis()->SetRangeUser(0.0, max_chi2 * 1.2);
  TCanvas *c5 = new TCanvas("c5", "", 700, 500);
  chi2_block_graph->Draw("AP");
  gPad->Update();
  for (int j = 1; j < 6; j++) {
    TLine* line = new TLine(j*16, c5->GetUymin(), j*16, c5->GetUymax());
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");
  }
  c5->SetGridy();
  c5->SaveAs(Form("%s/chi2_block.png", file_prefix));

  int n_conv = 0;
  printf("<-- THESE FITS CONVERGED -->\n");
  for (int i = 0; i < 384; i++) {
    if (fit_success[i]) {
      printf("channel %3i: gap %2.2f ± %3.2f, ratio %4.3f, chi2 %8.2f\n", i, other_gaps[i][0], other_gaps[i][1], first_ampl[i][0]/other_ampl[i][0][0], f_singlepixels[i]->GetChisquare()/f_singlepixels[i]->GetNDF()); 
      n_conv++;
    }
  }
  printf("<-- THESE FITS DO NOT CONVERGE -->\n");
  for (int i = 0; i < 384; i++) {
    if (!fit_success[i]) {
      printf("channel %3i: gap %2.2f ± %3.2f, ratio %4.3f, chi2 %8.2f\n", i, other_gaps[i][0], other_gaps[i][1], first_ampl[i][0]/other_ampl[i][0][0], f_singlepixels[i]->GetChisquare()/f_singlepixels[i]->GetNDF()); 
    }
  }
  printf("%i/384 fits converged!\n", n_conv);
  //printf("max chisqr_ndf: %f (channel %i); min chisqr_ndf: %f (channel %i)\n", max_chisqr_ndf, max_chisqr_idx, min_chisqr_ndf, min_chisqr_idx);

  /*
  // compute averages for IBs 0-2 and 3-5
  Double_t avg_gap = 0;
  for (int i = 0; i < 192; i++) {
    avg_gap += other_gaps[i][0];
  }
  Double_t ib0_2_avg_gap = avg_gap / 192;
  avg_gap = 0;
  for (int i = 192; i < 384; i++) {
    avg_gap += other_gaps[i][0];
  }
  Double_t ib3_5_avg_gap = avg_gap / 192;
  printf("IBs 0-2: avg sp gap = %f; IBs 3-5 avg sp gap = %f\n", ib0_2_avg_gap, ib3_5_avg_gap);
  */

  free(file_prefix);
}