#include <TROOT.h>
#include <TFile.h>
#include <TColor.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TMarker.h>
#include <TExec.h>
#include <TLatex.h>
#include <THStack.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "csvFile.h"
#include <exception>

// the root of all evil

const int num_bins = 30;
const double x_min = 0;
const double x_max = 600;

const std::vector<int> sectors {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 17};
const std::vector<int> interface_boards {0, 1, 2, 3, 4, 5};

std::vector<Color_t> hist_colors = {};
int hist_color_idx = 0;
std::vector<std::vector<double>> kelly_colors = {
  {255.0/255, 179.0/255, 0.0/255}, // vivid_yellow
  {128.0/255, 62.0/255, 117.0/255}, // strong_purple
  {255.0/255, 104.0/255, 0.0/255}, // vivid_orange
  {166.0/255, 189.0/255, 215.0/255}, // very_light_blue
  {193.0/255, 0.0/255, 32.0/255}, // vivid_red
  {206.0/255, 162.0/255, 98.0/255}, // grayish_yellow
  {129.0/255, 112.0/255, 102.0/255}, // medium_gray
  // these aren't good for people with defective color vision:
  {0.0/255, 125.0/255, 52.0/255}, // vivid_green
  {246.0/255, 118.0/255, 142.0/255}, // strong_purplish_pink
  {0.0/255, 83.0/255, 138.0/255}, // strong_blue
  {255.0/255, 122.0/255, 92.0/255}, // strong_yellowish_pink
  {83.0/255, 55.0/255, 122.0/255}, // strong_violet
  {255.0/255, 142.0/255, 0.0/255}, // vivid_orange_yellow
  {179.0/255, 40.0/255, 81.0/255}, // strong_purplish_red
  {244.0/255, 200.0/255, 0.0/255}, // vivid_greenish_yellow
  {127.0/255, 24.0/255, 13.0/255}, // strong_reddish_brown
  {147.0/255, 170.0/255, 0.0/255}, // vivid_yellowish_green
  {89.0/255, 51.0/255, 21.0/255}, // deep_yellowish_brown
  {241.0/255, 58.0/255, 19.0/255}, // vivid_reddish_orange
  {35.0/255, 44.0/255, 22.0/255} // dark_olive_green
};

std::vector<std::string> split_string(std::string s, std::string delimiter) {
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;
  while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
    token = s.substr (pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back (token);
  }
  res.push_back (s.substr (pos_start));
  return res;
}

std::vector<std::vector<double>> get_rainbow(int n) {
  std::vector<std::vector<double>> colors;
  for (int i = 0; i < n; i++) {
    double ratio = (double) i / (double) n;
    // std:: cout << "ratio=" << ratio << std::endl;
    double normalized = ratio * 256 * 6;
    // std:: cout << "norm=" << normalized << std::endl;
    //find the distance to the start of the closest region
    double x = normalized;
    while (x > 256) { x -= 256; }

    double red = 0, grn = 0, blu = 0;
    switch((int) normalized / 256) {
      case 0: red = 255;      grn = x;        blu = 0;       break;//red
      case 1: red = 255 - x;  grn = 255;      blu = 0;       break;//yellow
      case 2: red = 0;        grn = 255;      blu = x;       break;//green
      case 3: red = 0;        grn = 255 - x;  blu = 255;     break;//cyan
      case 4: red = x;        grn = 0;        blu = 255;     break;//blue
      case 5: red = 255;      grn = 0;        blu = 255 - x; break;//magenta
    }
    std::vector<double> rgb = {red/255, grn/255, blu/255};
    std::cout << "color: red=" << red << " green=" << grn << " blue=" << blu << std::endl;
    colors.push_back(rgb);
  }
  return colors;
}


void do_sipm() {
  gStyle->SetCanvasPreferGL(1);
  std::map<int, Color_t> sector_colors;
  for (int i = 0; i < sectors.size(); i++) {
    int sector = sectors[i];
    Color_t ci = TColor::GetFreeColorIndex();
    TColor *color = new TColor(ci, kelly_colors[i][0], kelly_colors[i][1], kelly_colors[i][2]);
    sector_colors[sector] = ci;
  }

  std::map<double, std::vector<double>> histogram_data;
  
  std::vector<double> hist_vop;
  std::vector<double> hist_mpv_err;
  std::vector<double> hist_mpv;

  std::fstream sector_map_file;
  sector_map_file.open("sector_maps.csv",std::ios::in);
  std::map<int, std::vector<std::string>> sector_map;
  std::vector<std::string> row;
  std::string line, word, temp;
  int line_num = 1;
  int sector_num = -1;
  int this_row_num = -1;
  int offset = -1;
  std::cout << "reading sector map" << std::endl;
  while (std::getline(sector_map_file, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    // read an entire row and
    // store it in a string variable 'line'
    // std::getline(sector_map_file, line);
    // used for breaking words
    std::stringstream s(line);
    // read every column data of a row and
    // store it in a string variable, 'word'
    while (std::getline(s, word, ',')) {
      // add all the column data
      // of a row to a vector
      row.push_back(word);
    }
    // now row has the data for this row
    if (line_num == 1) {
      // don't do anything with the header
      line_num++;
      continue;
    } else if (row[0].compare("") != 0) {
      // this is the first line of a new sector
      std::vector<std::string> split = split_string(row[0], "_");
      sector_num = std::stoi(split[1]);
      sector_map[sector_num] = std::vector<std::string>(96);
      for (int i = 0; i < 96; i++) {
        sector_map[sector_num][i] = "null";
      }
      offset = 0;
      // std::cout << "now sector is " << sector_num << std::endl;
    } // else this is a continuation of the previous sector
    sector_map[sector_num][offset] = row[1];
    sector_map[sector_num][offset + 1] = row[2];
    sector_map[sector_num][offset + 2] = row[3];
    sector_map[sector_num][offset + 3] = row[4];
    line_num++;
    offset += 4;
  }
  
  std::fstream sipm_file;
  sipm_file.open("new_vop.csv",std::ios::in);
  std::map<int, std::vector<double>> sipm_map;
  // std::vector<std::string> row;
  // std::string line, word, temp;
  line_num = 1;
  this_row_num = -1;
  offset = 0;
  int block_idx = -1;
  int possible_sectors = 0;
  std::vector<int> sector_has_data;
  std::map<int, int> sector_nums;
  std::map<int, double> sector_sums;
  std::cout << "reading sipm map" << std::endl;
  while (std::getline(sipm_file, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    // read an entire row and
    // store it in a string variable 'line'
    // std::getline(sipm_file, line);
    // used for breaking words
    std::stringstream s(line);
    // read every column data of a row and
    // store it in a string variable, 'word'
    while (std::getline(s, word, ',')) {
      // add all the column data
      // of a row to a vector
      row.push_back(word);
    }
    // now row has the data for this row
    if (line_num <= 1 || line_num == 3) {
      // don't do anything with the header
      line_num++;
      continue;
    } else if (line_num == 2) {
      for (int idx = 3; idx < row.size(); idx++) {
        std::vector<std::string> split = split_string(row[idx], " ");
        int sector = std::stoi(split[1]);
        sector_nums[idx] = sector;
        sector_sums[idx] = 0;
      }
      line_num++;
      continue;
    } else if (line_num == 4) {
      for (int idx = 3; idx < row.size(); idx++) {
        bool is_double = true;
        try {
          double x = std::stod(row[idx]);
        } catch (std::exception &e) {
          is_double = false;
        }
        if (is_double) {
          sector_has_data.push_back(idx);
          sipm_map[sector_nums[idx]] = std::vector<double>(96);
          std::cout << "sector " << sector_nums[idx] << " has data" << std::endl;
        }
      }
      block_idx = 3;
    }
    // std::cout << "block_idx = " << block_idx << std::endl;

    for (int idx : sector_has_data) {
      try {
        double channel_vop = std::stod(row[idx]);
        sector_sums[idx] += channel_vop;
      } catch (std::exception& e) {
        std::cout << "error:" << e.what() << " while converting [" << row[idx] << "] to double" << std::endl;
      }
    }

    offset++;
    if (offset == 4) {
      for (int idx : sector_has_data) {
        sipm_map[sector_nums[idx]][block_idx] = sector_sums[idx] / 4;
        sector_sums[idx] = 0;
      }
      if (block_idx % 4 == 0) {
        block_idx += 8;
      }
      offset = 0;
      block_idx--;
    }
    line_num++;
  }

  std::vector<double> all_vop;
  std::vector<double> all_mpv;
  std::map<std::string, double> dbn_mpv;

  std::map<int, std::map<double, std::vector<double>>> sector_vop_mpv;
  std::map<double, std::map<int, std::vector<double>>> vop_sector_mpv;

  for (int sector : sectors) {
    std::cout << "SECTOR " << sector << ":" << std::endl;
    std::stringstream sectorFileName;
    sectorFileName << "./sector_data/sector" << sector << ".root";
    TFile* sectorFile = new TFile(sectorFileName.str().c_str());
    TH1D* data;
    std::stringstream blockFileName;
    blockFileName << "h_run" << sector << "_block;1";
    // std::cout << "getting object" << std::endl;
    sectorFile->GetObject(blockFileName.str().c_str(), data);
    // std::cout << "got object" << std::endl;
    for (int i = 0; i < data->GetNcells(); i++) {
      // std::cout << "reading bin " << i << std::endl;
      double content = data->GetBinContent(i);
      int block_num = data->GetBinLowEdge(i);
      // first check if this is data we are interested in...
      if (sipm_map.find(sector) == sipm_map.end()) {
        std::cout << "unable to find sector " << sector << " in sipm_map" << std::endl;
        continue;
      } else if (block_num < 0 || block_num >= sipm_map[sector].size()) {
        std::cout << "block " << block_num << " was out of range for sipm_map sector " << sector << std::endl;
        continue;
      } else if (sector_map.find(sector) == sector_map.end()) {
        std::cout << "unable to find sector " << sector << " in sector_map" << std::endl;
        continue;
      } else if (block_num < 0 || block_num >= sector_map[sector].size()) {
        std::cout << "block " << block_num << " was out of range for sector_map sector " << sector << std::endl;
        continue;
      } else if (sector_map[sector][block_num][0] == 'F' || std::stoi(sector_map[sector][block_num]) >= 10000) {
        std::cout << "block " << block_num << ": rejected fudan block " << std::endl;
        continue;
      } else if (content <= 0 || content >= 1000) {
        std::cout << "block " << block_num << ": rejected bin content " << content << std::endl;
        continue;
      }
      // i guess this is valid data?
      double vop = sipm_map[sector][block_num];
      all_vop.push_back(vop);
      all_mpv.push_back(content);
      histogram_data[vop].push_back(content);
      dbn_mpv[sector_map[sector][block_num]] = content;
      sector_vop_mpv[sector][vop].push_back(content);
      vop_sector_mpv[vop][sector].push_back(content);
      std::cout << "block " << block_num << " (DBN " << std::stoi(sector_map[sector][block_num]) << "): good data (" << vop << ", " << content  << ")" << std::endl;
    }
  }

  int num_vops = histogram_data.size();
  std::cout << "there are " << num_vops << " vops" << std::endl;
  std::map<double, Color_t> vop_colors;
  int i = 0;
  for (const auto & p : histogram_data) {
    double vop = p.first;
    Color_t ci = TColor::GetFreeColorIndex();
    TColor *color = new TColor(ci, kelly_colors[i][0], kelly_colors[i][1], kelly_colors[i][2]);
    vop_colors[vop] = ci;
    i++;
  }

  THStack *hs = new THStack("hs", "distributions of mpv for each vop");
  csvfile vop_mpv("vop_mpv_fit_parameters.csv");
  vop_mpv << "vop" << "mean" << "mean err" << "sigma" << "sigma err" << csvfile::endrow;

  for (auto const &p : histogram_data) {
    TCanvas* canvas = new TCanvas();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    std::stringstream title;
    double vop = p.first;
    title << "vop " << vop << ";mpv;n_blocks";
    TH1D* hist = new TH1D("hist", title.str().c_str(), num_bins, x_min, x_max);
    hist->GetSumw2();
    for (double mpv : p.second) {
      hist->Fill(mpv);
    }
    hist->Fit("gaus", "Q");
    Color_t color = vop_colors[vop];
    hist_colors.push_back(color);
    hist->SetLineColor(color);
    hist->GetFunction("gaus")->SetLineColor(color);
    hist->Draw("E0 SAME");
    hs->Add(hist);

    std::stringstream fileName;
    fileName << "./vop_graphs/vop" << vop << ".svg";
    canvas->SaveAs(fileName.str().c_str());
    TF1* gaus = (TF1*) hist->GetListOfFunctions()->FindObject("gaus");
    double mean = gaus->GetParameter(1);
    double mean_err = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    double sigma_err = gaus->GetParError(2);
    vop_mpv << vop << mean << mean_err << sigma << sigma_err << csvfile::endrow;

    hist_vop.push_back(vop);
    hist_mpv_err.push_back(hist->GetStdDev());
    hist_mpv.push_back(hist->GetMean());
    hist_color_idx++;
  }
  TCanvas* canvas3 = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  hs->Draw("nostack");
  canvas3->BuildLegend();
  canvas3->SaveAs("./vop_graphs/vop_stacked.svg");

  std::cout << "creating vop by sector histograms" << std::endl;
  for (const auto & p : sector_vop_mpv) {
    int sector = p.first;
    std::map<double, std::vector<double>> vop_mpv_map = p.second;
    std::stringstream sector_title;
    sector_title << "sector " << sector << ";mpv;n_blocks";
    THStack *hs = new THStack("hs", sector_title.str().c_str());
    for (const auto & q : vop_mpv_map) {
      double vop = q.first;
      std::vector<double> mpvs = q.second;
      TCanvas* vop_canvas = new TCanvas("vop", "vop");
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(1);
      std::stringstream vop_title;
      vop_title << "sector " << sector << ", vop " << vop << ";mpv;n_blocks";
      TH1D* vop_hist = new TH1D("vop_hist", vop_title.str().c_str(), num_bins, x_min, x_max);
      vop_hist->GetSumw2();
      std::stringstream stack_title;
      stack_title << "vop " << vop;
      TH1D* stack_hist = new TH1D("stack_hist", stack_title.str().c_str(), num_bins, x_min, x_max);
      stack_hist->GetSumw2();
      for (double mpv : mpvs) {
        vop_hist->Fill(mpv);
        stack_hist->Fill(mpv);
      }
      vop_hist->SetLineColor(vop_colors[vop]);
      vop_hist->Fit("gaus", "Q");
      vop_hist->GetFunction("gaus")->SetLineColor(vop_colors[vop]);
      vop_hist->Draw("E0 SAME");
      std::stringstream fileName;
      fileName << "./vop_by_sector/sector-" << sector << "_vop-" << vop << ".svg";
      vop_canvas->SaveAs(fileName.str().c_str());

      stack_hist->SetFillColorAlpha(vop_colors[vop], 0.5);
      stack_hist->SetLineColorAlpha(vop_colors[vop], 0.0);
      stack_hist->Fit("gaus", "Q");
      stack_hist->GetFunction("gaus")->SetLineColor(vop_colors[vop]);
      hs->Add(stack_hist);
    }
    TCanvas* sector_canvas = new TCanvas("sector", "sector");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hs->Draw("nostackb");
    std::stringstream fileName;
    fileName << "./vop_by_sector/sector-" << sector << "_all.svg";
    sector_canvas->BuildLegend();
    sector_canvas->SaveAs(fileName.str().c_str());
  }

  std::cout << "creating the other histograms" << std::endl;
  for (const auto & p : vop_sector_mpv) {
    double vop = p.first;
    std::map<int, std::vector<double>> sector_mpv_map = p.second;
    std::stringstream vop_title;
    vop_title << "vop " << vop << ";mpv;n_blocks";
    THStack *hs = new THStack("hs", vop_title.str().c_str());
    for (const auto & q : sector_mpv_map) {
      int sector = q.first;
      std::vector<double> mpvs = q.second;
      TCanvas* sector_canvas = new TCanvas();
      // gStyle->SetOptStat(0);
      // gStyle->SetOptFit(1);
      std::stringstream sector_title;
      sector_title << "sector " << sector << ", vop " << vop << ";mpv;n_blocks";
      TH1D* vop_hist = new TH1D("vop_hist", sector_title.str().c_str(), num_bins, x_min, x_max);
      vop_hist->GetSumw2();
      std::stringstream stack_title;
      stack_title << "sector " << sector;
      TH1D* stack_hist = new TH1D("stack_hist", stack_title.str().c_str(), num_bins, x_min, x_max);
      stack_hist->GetSumw2();
      for (double mpv : mpvs) {
        vop_hist->Fill(mpv);
        stack_hist->Fill(mpv);
      }
      // vop_hist->SetLineColor(vop_colors[vop]);
      // vop_hist->Fit("gaus", "Q");
      // vop_hist->GetFunction("gaus")->SetLineColor(vop_colors[vop]);
      // vop_hist->Draw("E0 SAME");
      // std::stringstream fileName;
      // fileName << "./vop_by_sector/sector-" << sector << "_vop-" << vop << ".png";
      // sector_canvas->SaveAs(fileName.str().c_str());

      stack_hist->SetFillColorAlpha(sector_colors[sector], 0.5);
      stack_hist->SetLineColorAlpha(sector_colors[sector], 0.0);
      stack_hist->Fit("gaus", "Q");
      stack_hist->GetFunction("gaus")->SetLineColor(sector_colors[sector]);
      hs->Add(stack_hist);
    }
    TCanvas* vop_canvas = new TCanvas();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hs->Draw("nostackb");
    std::stringstream fileName;
    fileName << "./the_other_ones/vop-" << vop << "_all.svg";
    
    vop_canvas->BuildLegend();
    vop_canvas->SaveAs(fileName.str().c_str());
  }

  // scatter plot
  TCanvas* canvas = new TCanvas();
  TGraph* graph = new TGraph(all_vop.size(), &all_vop[0], &all_mpv[0]);
  graph->SetTitle("block mpv vs sipm vop");
  graph->GetXaxis()->SetTitle("vop");
  graph->GetYaxis()->SetTitle("mpv");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.5);
  graph->Draw("AP");
  // canvas->Update();
  canvas->SaveAs("./vop_graphs/vop_mpv_scatter.svg");

  std::vector<double> err_x;
  for (int i = 0; i < hist_vop.size(); i++) {
    err_x.push_back(0);
  }

  // discrete graph
  TCanvas* canvas2 = new TCanvas();
  TGraphErrors* graph2 = new TGraphErrors(hist_vop.size(), &hist_vop[0], &hist_mpv[0], &err_x[0], &hist_mpv_err[0]);
  graph2->SetTitle("mean block mpv vs sipm vop");
  graph2->GetXaxis()->SetTitle("vop");
  graph2->GetYaxis()->SetTitle("mpv");
  graph2->SetMarkerStyle(21);
  TExec *ex = new TExec("ex","DrawCol();");
  graph2->GetListOfFunctions()->Add(ex);
  graph2->Draw("AP");
  canvas2->SaveAs("./vop_graphs/vop_mpv_discrete.svg");

  csvfile csv("dbn_to_mpv.csv");
  csv << "dbn" << "mpv" << csvfile::endrow;
  for (auto const &p : dbn_mpv) {
    csv << std::stoi(p.first) << p.second << csvfile::endrow;
  }
}

void DrawCol()
{
   Int_t i,n;
   Double_t x,y;
   TLatex *l;

   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   n = g->GetN();
   TMarker *m;
   for (i=0; i<n; i++) {
      g->GetPoint(i,x,y);
      m = new TMarker(x,y,20);
      m->SetMarkerColor(hist_colors[i]);
      m->Paint();
   }
}
