#include <TFile.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "csvFile.h"

const int mean_num_bins = 20;
const double mean_x_min = 0;
const double mean_x_max = 600;
const int sigma_num_bins = 20;
const double sigma_x_min = 0;
const double sigma_x_max = 100;
const int sigma_mean_num_bins = 20;
const double sigma_mean_x_min = 0;
const double sigma_mean_x_max = 0.4;
const int tower_num_bins = 20;
const double tower_x_min = 0;
const double tower_x_max = 600;
const std::vector<int> sectors {1, 2, 3, 4, 10, 16, 17};
const std::vector<int> interface_boards {0, 1, 2, 3, 4, 5};

void do_the_things() {
  csvfile csv("fit_data.csv");
  csv << "sector" << "ib" << "chi2" << "ndf" << "const" << "const_err" << "mean" << "mean_err" << "sigma" << "sigma_err" << csvfile::endrow;
  TH1D* mean_hist = new TH1D("mean_hist", "Distribution of Fit Mean;Mean (mpv);Count", mean_num_bins, mean_x_min, mean_x_max);
  TH1D* sigma_hist = new TH1D("sigma_hist", "Distribution of Fit Sigma;Sigma (mpv);Count", sigma_num_bins, sigma_x_min, sigma_x_max);
  TH1D* sigma_mean_hist = new TH1D("sigma_mean_hist", "Distribution of Fit Sigma/Mean;Sigma/Mean;Count", sigma_mean_num_bins, sigma_mean_x_min, sigma_mean_x_max);
  for (int sector : sectors) {
    std::stringstream sectorFileName;
    sectorFileName << "./sector_data/sector" << sector << ".root";
    TFile* sectorFile = new TFile(sectorFileName.str().c_str());
    for (int ib : interface_boards)  {
      TCanvas* canvas = new TCanvas();
      TH1D* data;
      std::stringstream ibFileName;
      ibFileName << "h_run" << sector << "_ib_" << ib << ";1";
      sectorFile->GetObject(ibFileName.str().c_str(), data);
      if(!data) { continue; }
      std::stringstream title;
      title << "Sector " << sector << " IB " << ib << ";mpv;ntowers";
      TH1D* hist = new TH1D("hist", title.str().c_str(), tower_num_bins, tower_x_min, tower_x_max);
      hist->GetSumw2();
      for (int i = 0; i < data->GetNcells(); i++) {
        double content = data->GetBinContent(i);
        // first check if this is data we are interested in...
        if (content > 0) {
          hist->Fill(content);
        }
      }
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(1);
      hist->Fit("gaus", "Q");
      hist->Draw("E0 SAME");
      std::stringstream pngFileName;
      pngFileName << "./graphs/png/sector" << sector << "_ib" << ib << ".png";
      canvas->SaveAs(pngFileName.str().c_str());
      std::stringstream svgFileName;
      svgFileName << "./graphs/svg/sector" << sector << "_ib" << ib << ".svg";
      canvas->SaveAs(svgFileName.str().c_str());
      std::stringstream pdfFileName;
      pdfFileName << "./graphs/pdf/sector" << sector << "_ib" << ib << ".pdf";
      canvas->SaveAs(pdfFileName.str().c_str());
      TF1* gaus = (TF1*) hist->GetListOfFunctions()->FindObject("gaus");
      double mean = gaus->GetParameter(1);
      double sigma = gaus->GetParameter(2);
      if (mean > 0 && abs(mean) < 1000) {
        mean_hist->Fill(mean);
      } else {
        std::cout << "rejected mean for sector " << sector << " ib " << ib << " (" << mean << ") for being out of range" << std::endl;
      }
      if (sigma < 500) {
        sigma_hist->Fill(sigma);
      } else {
        std::cout << "rejected sigma for sector " << sector << " ib " << ib << " (" << sigma << ") for being out of range" << std::endl;
      }
      if (sigma/mean > 0) {
        sigma_mean_hist->Fill(sigma/mean);
      } else {
        std::cout << "rejected sigma/mean for sector " << sector << " ib " << ib << " (" << sigma/mean << ") for being out of range" << std::endl;
      }
      csv << sector << ib << gaus->GetChisquare() << gaus->GetNDF() << gaus->GetParameter(0) << gaus->GetParError(0) << gaus->GetParameter(0) << gaus->GetParError(1) << gaus->GetParameter(2) << gaus->GetParError(2) << csvfile::endrow;
    }
  }
  // fit for mean and sigma histograms
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas* mean_canvas = new TCanvas();
  mean_hist->Fit("gaus");
  mean_hist->Draw("E0 SAME");
  mean_canvas->SaveAs("distr_fit_mean.svg");
  mean_canvas->SaveAs("distr_fit_mean.png");
  TCanvas* sigma_canvas = new TCanvas();
  sigma_hist->Fit("gaus");
  sigma_hist->Draw("E0 SAME");
  sigma_canvas->SaveAs("distr_fit_sigma.svg");
  sigma_canvas->SaveAs("distr_fit_sigma.png");
  TCanvas* sigma_mean_canvas = new TCanvas();
  sigma_mean_hist->Fit("gaus");
  sigma_mean_hist->Draw("E0 SAME");
  sigma_mean_canvas->SaveAs("distr_fit_sigma_over_mean.svg");
  sigma_mean_canvas->SaveAs("distr_fit_sigma_over_mean.png");
}