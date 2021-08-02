
#include <TFile.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "csvFile.h"

const int num_bins = 20;
const double x_min = 0;
const double x_max = 600;

void test() {
  std::stringstream sectorName;
  sectorName << "./sector_data/sector" << 1 << ".root";
  TFile* myFile = new TFile(sectorName.str().c_str());
  TCanvas* c1 = new TCanvas();
  TH1D* data;
  std::stringstream fileName;
  fileName << "h_run1_ib_" << 0 << ";1";
  myFile->GetObject(fileName.str().c_str(), data);
  TH1D* hist = new TH1D("h", "title;x axis;y axis", num_bins, x_min, x_max);
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
  hist->Fit("gaus");
  hist->Draw("E0 SAME");
  std::stringstream graphFileName;
  graphFileName << "./graphs/graph.svg";
  //c1->SaveAs(graphFileName.str().c_str());
  //std::cout << "mean: " << hist->GetMean() << std::endl;
  //std::cout << "std dev: " << hist->GetStdDev() << std::endl;
  TF1* gaus = (TF1*) hist->GetListOfFunctions()->FindObject("gaus");
  
  std::cout << "number of parameters: " << gaus->GetNumberFreeParameters() << std::endl;
  std::cout << gaus->GetParName(0) << ": " << gaus->GetParameter(0) << " +/- " << gaus->GetParError(0) << std::endl;
  std::cout << gaus->GetParName(1) << ": " << gaus->GetParameter(1) << " +/- " << gaus->GetParError(1) << std::endl;
  std::cout << gaus->GetParName(2) << ": " << gaus->GetParameter(2) << " +/- " << gaus->GetParError(2) << std::endl;
  std::cout << std::endl << std::stoi("") << std::endl;
}