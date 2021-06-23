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
#include <algorithm>
#include "csvFile.h"
#include <exception>

// the root of all evil

const int num_bins = 30;
const double x_min = 0;
const double x_max = 600;

const std::vector<int> sectors {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17};
const std::vector<int> interface_boards {0, 1, 2, 3, 4, 5};

std::map<double, std::string> vop_to_board_series = {
  {69.13, "A"},
  {69.18, "B"},
  {69.08, "C"},
  {69.23, "D"},
  {69.03, "E"},
  {69.28, "F"},
  {69.33, "G"},
  {68.93, "H"},
  {69.38, "I"},
  {68.88, "J"},
  {69.43, "K"},
  {68.83, "L"},
  {69.48, "M"},
  {68.78, "N"},
  {69.53, "O"},
  {68.58, "P"},
  {68.63, "Q"},
  //{, "R"},
  {68.68, "S"}
};

std::vector<Color_t> hist_colors;

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


void incl_fiber_batch() {
  
  
  
  gStyle->SetCanvasPreferGL(1);
  std::map<int, Color_t> sector_colors;
  for (int i = 0; i < sectors.size(); i++) {
    int sector = sectors[i];
    Color_t ci = TColor::GetFreeColorIndex();
    TColor *color = new TColor(ci, kelly_colors[i][0], kelly_colors[i][1], kelly_colors[i][2]);
    sector_colors[sector] = ci;
  }

  std::fstream sector_map_file;
  sector_map_file.open("sector_maps.csv", std::ios::in);
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
  
  std::map<std::string, int> dbn_to_fiber_batch;

  std::fstream blocks_1_12;
  blocks_1_12.open("blocks_1-12.csv",std::ios::in);
  // std::vector<std::string> row;
  // std::string line, word, temp;
  line_num = 1;
  this_row_num = -1;
  offset = 0;
  std::cout << "reading blocks_1-12" << std::endl;
  while (std::getline(blocks_1_12, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    // read an entire row and
    // store it in a string variable 'line'
    // std::getline(new_vop_file, line);
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
    if (line_num == 0) {
      // don't do anything with the header
      line_num++;
      continue;
    }
    // check if this is a block...
    bool is_int = true;
    int block_type;
    try {
      block_type = std::stoi(row[2]);
    } catch (std::exception &e) {
      is_int = false;
    }
    if (is_int) {
      // this row has block data
      std::string dbn = row[0];
      std::vector<std::string> split = split_string(row[9], "-");
      if (split.size() < 2) {
        std::cout << "skipped fiber batch for DBN " << dbn << " since it does not contain a '-'" << std::endl;
      } else {
        is_int = true;
        int fiber_batch = 0;
        try {
          fiber_batch = std::stoi(split[0]);
        } catch (std::exception &e) {
          is_int = false;
        }
        if (is_int) {
          dbn_to_fiber_batch[dbn] = fiber_batch;
        } else {
          std::cout << "error at DBN " << dbn << ": '" << split[0] << "' not castable to int" << std::endl; 
        }
      }
    }
    line_num++;
  }

  std::fstream blocks_13_64;
  blocks_13_64.open("blocks_13-64.csv",std::ios::in);
  // std::vector<std::string> row;
  // std::string line, word, temp;
  line_num = 1;
  this_row_num = -1;
  offset = 0;
  std::cout << "reading blocks_13-64" << std::endl;
  while (std::getline(blocks_13_64, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    // read an entire row and
    // store it in a string variable 'line'
    // std::getline(new_vop_file, line);
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
    if (line_num == 0) {
      // don't do anything with the header
      line_num++;
      continue;
    }
    // check if this is a block...
    bool is_int = true;
    int block_type;
    try {
      block_type = std::stoi(row[2]);
    } catch (std::exception &e) {
      is_int = false;
    }
    if (is_int) {
      // this row has block data
      std::string dbn = row[0];
      std::vector<std::string> split = split_string(row[9], "-");
      if (split.size() < 2) {
        std::cout << "skipped fiber batch for DBN " << dbn << " since it does not contain a '-'" << std::endl;
      } else {
        is_int = true;
        int fiber_batch = 0;
        try {
          fiber_batch = std::stoi(split[0]);
        } catch (std::exception &e) {
          is_int = false;
        }
        if (is_int) {
          dbn_to_fiber_batch[dbn] = fiber_batch;
        } else {
          std::cout << "error at DBN " << dbn << ": '" << split[0] << "' not castable to int" << std::endl; 
        }
      }
    }
    line_num++;
  }

  /*
  std::cout << "DBN to fiber batch:";
  for (const auto &p : dbn_to_fiber_batch) { 
    std::cout << "DBN " << p.first << ": " << p.second << "; ";
  }
  std::cout <<std::endl;
  */

 std::map<int, double> fiber_batch_to_scale_factor;
  std::fstream fiber_batch_map;
  fiber_batch_map.open("fiberbatchmeans.csv",std::ios::in);
  // std::vector<std::string> row;
  // std::string line, word, temp;
  line_num = 1;
  this_row_num = -1;
  offset = 0;
  std::cout << "reading fiber batch means" << std::endl;
  while (std::getline(fiber_batch_map, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    // read an entire row and
    // store it in a string variable 'line'
    // std::getline(new_vop_file, line);
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
    if (line_num == 0) {
      // don't do anything with the header
      line_num++;
      continue;
    }
    // check if this is a block...
    bool is_double = true;
    double scale_factor;
    try {
      scale_factor = std::stod(row[2]);
    } catch (std::exception &e) {
      is_double = false;
    }
    if (is_double) {
      int fiber_batch = std::stoi(row[0]);
      fiber_batch_to_scale_factor[fiber_batch] = scale_factor;
    }
    line_num++;
  }

  std::fstream new_vop_file;
  new_vop_file.open("new_vop.csv",std::ios::in);
  std::map<int, std::vector<double>> new_sipm_map;
  // std::vector<std::string> row;
  // std::string line, word, temp;
  line_num = 1;
  this_row_num = -1;
  offset = 0;
  int block_idx = -1;
  int possible_sectors = 0;
  std::vector<int> new_sector_has_data;
  std::map<int, int> new_sector_nums;
  std::map<int, double> new_sector_sums;
  std::cout << "reading new vop data" << std::endl;
  while (std::getline(new_vop_file, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    // read an entire row and
    // store it in a string variable 'line'
    // std::getline(new_vop_file, line);
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
        new_sector_nums[idx] = sector;
        new_sector_sums[idx] = 0;
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
          new_sector_has_data.push_back(idx);
          new_sipm_map[new_sector_nums[idx]] = std::vector<double>(96);
          std::cout << "sector " << new_sector_nums[idx] << " has new vop data" << std::endl;
        }
      }
      block_idx = 3;
    }
    // std::cout << "block_idx = " << block_idx << std::endl;

    for (int idx : new_sector_has_data) {
      try {
        double channel_vop = std::stod(row[idx]);
        new_sector_sums[idx] += channel_vop;
      } catch (std::exception& e) {
        std::cout << "error:" << e.what() << " while converting [" << row[idx] << "] to double" << std::endl;
      }
    }

    offset++;
    if (offset == 4) {
      for (int idx : new_sector_has_data) {
        new_sipm_map[new_sector_nums[idx]][block_idx] = new_sector_sums[idx] / 4;
        new_sector_sums[idx] = 0;
      }
      if (block_idx % 4 == 0) {
        block_idx += 8;
      }
      offset = 0;
      block_idx--;
    }
    line_num++;
  }

  // NEW VOP DATA
  std::vector<double> all_mpv;
  std::map<std::string, double> dbn_mpv;
  std::vector<double> new_all_vop;
  std::map<int, std::map<double, std::vector<double>>> new_sector_vop_mpv;
  std::map<double, std::map<int, std::vector<double>>> new_vop_sector_mpv;
  std::map<double, std::vector<double>> new_histogram_data;
  std::vector<double> new_hist_vop;
  std::vector<double> new_hist_mpv_err;
  std::vector<double> new_hist_mpv;

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
      if (new_sipm_map.find(sector) == new_sipm_map.end()) {
        std::cout << "unable to find sector " << sector << " in new_sipm_map" << std::endl;
        continue;
      } else if (block_num < 0 || block_num >= new_sipm_map[sector].size()) {
        std::cout << "block " << block_num << " was out of range for new_sipm_map sector " << sector << std::endl;
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
      double vop = new_sipm_map[sector][block_num];
      new_all_vop.push_back(vop);
      all_mpv.push_back(content);
      new_histogram_data[vop].push_back(content);
      dbn_mpv[sector_map[sector][block_num]] = content;
      new_sector_vop_mpv[sector][vop].push_back(content);
      new_vop_sector_mpv[vop][sector].push_back(content);
      //std::cout << "block " << block_num << " (DBN " << std::stoi(sector_map[sector][block_num]) << "): good data (" << vop << ", " << content  << ")" << std::endl;
    }
  }

  // got all mpv data...calculate 
  int num_dbns_with_mpv = 0;
  double sum_for_mpv_mean = 0;
  for (double mpv : all_mpv) {
    num_dbns_with_mpv++;
    sum_for_mpv_mean += mpv;
  }
  double mean_mpv = sum_for_mpv_mean / num_dbns_with_mpv;
  std::cout << "MEAN MPV: " << mean_mpv << std::endl;

  int new_num_vops = new_histogram_data.size();
  std::cout << "there are " << new_num_vops << " new vops" << std::endl;
  std::map<double, Color_t> new_vop_colors;
  int new_i = 0;
  for (const auto & p : new_histogram_data) {
    double vop = p.first;
    Color_t ci = TColor::GetFreeColorIndex();
    TColor *color = new TColor(ci, kelly_colors[new_i][0], kelly_colors[new_i][1], kelly_colors[new_i][2]);
    new_vop_colors[vop] = ci;
    new_i++;
  }

  THStack *new_hs = new THStack("hs", "distributions of mpv for each vop");
  csvfile new_vop_mpv("new_vop_mpv_fit_parameters.csv");
  new_vop_mpv << "vop" << "mean" << "mean err" << "sigma" << "sigma err" << csvfile::endrow;

  std::map<double, std::pair<double, double>> tim_y_axis_data;

  for (auto const &p : new_histogram_data) {
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
    Color_t color = new_vop_colors[vop];
    hist_colors.push_back(color);
    hist->SetLineColor(color);
    hist->GetFunction("gaus")->SetLineColor(color);
    hist->Draw("E0 SAME");
    new_hs->Add(hist);

    std::stringstream fileName;
    fileName << "./new_vop/vop_graphs/vop" << vop << ".svg";
    canvas->SaveAs(fileName.str().c_str());
    TF1* gaus = (TF1*) hist->GetListOfFunctions()->FindObject("gaus");
    double mean = gaus->GetParameter(1);
    double mean_err = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    double sigma_err = gaus->GetParError(2);
    new_vop_mpv << vop << mean << mean_err << sigma << sigma_err << csvfile::endrow;

    new_hist_vop.push_back(vop);
    new_hist_mpv_err.push_back(mean_err); //hist_mpv_err.push_back(hist->GetStdDev());
    new_hist_mpv.push_back(mean);
    tim_y_axis_data[vop] = std::pair<double, double>(mean, mean_err);
  }
  TCanvas* new_canvas3 = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  new_hs->Draw("nostack");
  new_canvas3->BuildLegend();
  new_canvas3->SaveAs("./new_vop/vop_graphs/vop_stacked.svg");

  std::cout << "creating vop by sector histograms" << std::endl;
  for (const auto & p : new_sector_vop_mpv) {
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
      vop_hist->SetLineColor(new_vop_colors[vop]);
      vop_hist->Fit("gaus", "Q");
      vop_hist->GetFunction("gaus")->SetLineColor(new_vop_colors[vop]);
      vop_hist->Draw("E0 SAME");
      std::stringstream fileName;
      fileName << "./new_vop/sector_by_vop/sector-" << sector << "_vop-" << vop << ".svg";
      vop_canvas->SaveAs(fileName.str().c_str());

      stack_hist->SetFillColorAlpha(new_vop_colors[vop], 0.5);
      stack_hist->SetLineColorAlpha(new_vop_colors[vop], 0.0);
      stack_hist->Fit("gaus", "Q");
      stack_hist->GetFunction("gaus")->SetLineColor(new_vop_colors[vop]);
      hs->Add(stack_hist);
    }
    TCanvas* sector_canvas = new TCanvas("sector", "sector");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hs->Draw("nostackb");
    std::stringstream fileName;
    fileName << "./new_vop/sector_by_vop/sector-" << sector << "_all.svg";
    sector_canvas->BuildLegend();
    sector_canvas->SaveAs(fileName.str().c_str());
  }

  std::cout << "creating the other histograms" << std::endl;
  for (const auto & p : new_vop_sector_mpv) {
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
    fileName << "./new_vop/vop_by_sector/vop-" << vop << "_all.svg";
    
    vop_canvas->BuildLegend();
    vop_canvas->SaveAs(fileName.str().c_str());
  }

  // scatter plot
  TCanvas* new_canvas = new TCanvas();
  TGraph* new_graph = new TGraph(new_all_vop.size(), &new_all_vop[0], &all_mpv[0]);
  new_graph->SetTitle("block mpv vs sipm vop");
  new_graph->GetXaxis()->SetTitle("vop");
  new_graph->GetYaxis()->SetTitle("mpv");
  new_graph->SetMarkerStyle(20);
  new_graph->SetMarkerSize(0.5);
  new_graph->Draw("AP");
  // canvas->Update();
  new_canvas->SaveAs("./new_vop/vop_graphs/vop_mpv_scatter.svg");

  std::vector<double> new_err_x;
  for (int i = 0; i < new_hist_vop.size(); i++) {
    new_err_x.push_back(0);
  }

  // discrete graph
  TCanvas* new_canvas2 = new TCanvas();
  TGraphErrors* new_graph2 = new TGraphErrors(new_hist_vop.size(), &new_hist_vop[0], &new_hist_mpv[0], &new_err_x[0], &new_hist_mpv_err[0]);
  new_graph2->SetTitle("Block MPV vs SiPM VOp (After)");
  new_graph2->GetXaxis()->SetTitle("vop");
  new_graph2->GetYaxis()->SetTitle("mpv");
  new_graph2->GetXaxis()->SetLimits(68.85, 69.4);
  new_graph2->GetYaxis()->SetLimits(220., 380.);
  new_graph2->SetMarkerStyle(21);
  //TExec *new_ex = new TExec("ex","DrawCol();");
  //new_graph2->GetListOfFunctions()->Add(new_ex);
  new_graph2->Draw("AP");
  new_canvas2->SetGrid();
  new_canvas2->SaveAs("./new_vop/vop_graphs/vop_mpv_discrete.svg");

  // tim's pcb series histogram
  //TFile* tim_outfile = new TFile("./new_vop/vop_graphs/tim_graphs.root","RECREATE");
  TCanvas* new_tim_canvas = new TCanvas();
  // first find all pcb series whcih have current data...
  std::vector<std::string> pcb_series_with_data;
  std::map<std::string, double> pcb_series_to_vop;
  for (int i = 0; i < new_hist_vop.size(); i++) {
    double vop = new_hist_vop[i];
    if (vop_to_board_series.find(vop) == vop_to_board_series.end()) {
      std::cout << "*unable to find vop " << vop << " in pcb series map" << std::endl;
    } else {
      //std::cout << "vop " << vop << " matches to board " << vop_to_board_series[vop] << std::endl;
      pcb_series_with_data.push_back(vop_to_board_series[vop]);
      pcb_series_to_vop[vop_to_board_series[vop]] = vop;
    }
  }
  int num_pcb_series_with_data = pcb_series_with_data.size();
  std::sort(pcb_series_with_data.begin(), pcb_series_with_data.end());

  TH1D* tim_empty_hist = new TH1D("h", "tim_empty", num_pcb_series_with_data, 0, num_pcb_series_with_data);
  tim_empty_hist->SetTitle("Mean Block MPV by SiPM PCB Series");
  tim_empty_hist->GetXaxis()->SetTitle("PCB Series");
  tim_empty_hist->GetYaxis()->SetTitle("MPV");
  //tim_empty_hist->GetYaxis()->SetRange(200, 400);
  tim_empty_hist->SetAxisRange(220., 380., "Y");
  for (int i = 0; i < num_pcb_series_with_data; i++) {
    tim_empty_hist->GetXaxis()->SetBinLabel(i + 1, pcb_series_with_data[i].c_str());
  }
  //tim_empty_hist->Write();
  tim_empty_hist->Draw();

  std::vector<double> tim_x_axis;
  for (int i = 0; i < num_pcb_series_with_data; i++) {
    tim_x_axis.push_back(i + 0.5);
  }
  std::vector<double> tim_x_axis_err;
  for (int i = 0; i < num_pcb_series_with_data; i++) { tim_x_axis_err.push_back(0); }
  std::vector<double> tim_y_axis;
  std::vector<double> tim_y_axis_err;
  for (std::string series : pcb_series_with_data) {
    double vop = pcb_series_to_vop[series];
    std::pair<double, double> p = tim_y_axis_data[vop];
    tim_y_axis.push_back(p.first);
    tim_y_axis_err.push_back(p.second);
    //std::cout << "y axis data for vop " << vop << ": (" << p.first << ", " << p.second << ")" << std::endl; 
  }
  TGraphErrors* tim_graph = new TGraphErrors(num_pcb_series_with_data, &tim_x_axis[0], &tim_y_axis[0], &tim_x_axis_err[0], &tim_y_axis_err[0]);

  tim_graph->GetYaxis()->SetTitle("MPV");
  tim_graph->GetYaxis()->SetLimits(220., 380.);
  tim_graph->GetYaxis()->SetTitle("MPV");
  tim_graph->SetMarkerStyle(21);
  //tim_graph->Write();
  tim_graph->Draw("P SAME");
  new_tim_canvas->SetGrid();
  new_tim_canvas->SaveAs("./new_vop/vop_graphs/vop_mpv_discrete_tim.svg");



  // ADJUSTED NEW VOP DATA
  std::vector<double> adjusted_all_mpv;
  std::map<std::string, double> adjusted_dbn_mpv;
  //std::vector<double> new_all_vop;
  std::map<int, std::map<double, std::vector<double>>> adjusted_sector_vop_mpv;
  std::map<double, std::map<int, std::vector<double>>> adjusted_vop_sector_mpv;
  std::map<double, std::vector<double>> adjusted_histogram_data;
  std::vector<double> adjusted_hist_vop;
  std::vector<double> adjusted_hist_mpv_err;
  std::vector<double> adjusted_hist_mpv;

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
      if (new_sipm_map.find(sector) == new_sipm_map.end()) {
        std::cout << "unable to find sector " << sector << " in new_sipm_map" << std::endl;
        continue;
      } else if (block_num < 0 || block_num >= new_sipm_map[sector].size()) {
        std::cout << "block " << block_num << " was out of range for new_sipm_map sector " << sector << std::endl;
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
      std::string dbn = sector_map[sector][block_num];
      int fiber_batch = dbn_to_fiber_batch[dbn];
      if (fiber_batch_to_scale_factor.find(fiber_batch) != fiber_batch_to_scale_factor.end()) {
        double correction_factor = fiber_batch_to_scale_factor[fiber_batch];
        double adjusted_content = content * correction_factor;
        //std::cout << "DBN " << std::stoi(dbn) << ": fiber batch " << fiber_batch << "; correction factor " << correction_factor
        //  << " (" << content << "->" << adjusted_content << ")" << std::endl;
        // i guess this is valid data?
        double vop = new_sipm_map[sector][block_num];
        //new_all_vop.push_back(vop);
        adjusted_all_mpv.push_back(content);
        adjusted_histogram_data[vop].push_back(content);
        adjusted_dbn_mpv[dbn] = content;
        adjusted_sector_vop_mpv[sector][vop].push_back(content);
        adjusted_vop_sector_mpv[vop][sector].push_back(content);
        //std::cout << "block " << block_num << " (DBN " << std::stoi(sector_map[sector][block_num]) << "): good data (" << vop << ", " << content  << ")" << std::endl;
      } else {
        std::cout << "** failed to add DBN " << std::stoi(dbn) << " since batch map does not contain batch " << fiber_batch << std::endl;
      }
    }
  }

  // got all mpv data...calculate 
  int adjusted_num_dbns_with_mpv = 0;
  double adjusted_sum_for_mpv_mean = 0;
  for (double mpv : adjusted_all_mpv) {
    adjusted_num_dbns_with_mpv++;
    adjusted_sum_for_mpv_mean += mpv;
  }
  double adjusted_mean_mpv = adjusted_sum_for_mpv_mean / adjusted_num_dbns_with_mpv;
  std::cout << "ADJUSTED MEAN MPV: " << adjusted_mean_mpv << std::endl;

  int adjusted_num_vops = adjusted_histogram_data.size();

  THStack *adjusted_hs = new THStack("hs", "distributions of mpv for each vop");
  csvfile adjusted_vop_mpv("adjusted_vop_mpv_fit_parameters.csv");
  adjusted_vop_mpv << "vop" << "mean" << "mean err" << "sigma" << "sigma err" << csvfile::endrow;

  std::map<double, std::pair<double, double>> adjusted_tim_y_axis_data;

  for (auto const &p : adjusted_histogram_data) {
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
    Color_t color = new_vop_colors[vop];
    hist_colors.push_back(color);
    hist->SetLineColor(color);
    hist->GetFunction("gaus")->SetLineColor(color);
    hist->Draw("E0 SAME");
    adjusted_hs->Add(hist);

    std::stringstream fileName;
    fileName << "./adjusted_vop/vop_graphs/vop" << vop << ".svg";
    canvas->SaveAs(fileName.str().c_str());
    TF1* gaus = (TF1*) hist->GetListOfFunctions()->FindObject("gaus");
    double mean = gaus->GetParameter(1);
    double mean_err = gaus->GetParError(1);
    double sigma = gaus->GetParameter(2);
    double sigma_err = gaus->GetParError(2);
    adjusted_vop_mpv << vop << mean << mean_err << sigma << sigma_err << csvfile::endrow;

    adjusted_hist_vop.push_back(vop);
    adjusted_hist_mpv_err.push_back(mean_err); //hist_mpv_err.push_back(hist->GetStdDev());
    adjusted_hist_mpv.push_back(mean);
    adjusted_tim_y_axis_data[vop] = std::pair<double, double>(mean, mean_err);
  }
  TCanvas* adjusted_canvas3 = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  adjusted_hs->Draw("nostack");
  adjusted_canvas3->BuildLegend();
  adjusted_canvas3->SaveAs("./adjusted_vop/vop_graphs/vop_stacked.svg");

  std::cout << "creating vop by sector histograms" << std::endl;
  for (const auto & p : adjusted_sector_vop_mpv) {
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
      vop_hist->SetLineColor(new_vop_colors[vop]);
      vop_hist->Fit("gaus", "Q");
      vop_hist->GetFunction("gaus")->SetLineColor(new_vop_colors[vop]);
      vop_hist->Draw("E0 SAME");
      std::stringstream fileName;
      fileName << "./adjusted_vop/sector_by_vop/sector-" << sector << "_vop-" << vop << ".svg";
      vop_canvas->SaveAs(fileName.str().c_str());

      stack_hist->SetFillColorAlpha(new_vop_colors[vop], 0.5);
      stack_hist->SetLineColorAlpha(new_vop_colors[vop], 0.0);
      stack_hist->Fit("gaus", "Q");
      stack_hist->GetFunction("gaus")->SetLineColor(new_vop_colors[vop]);
      hs->Add(stack_hist);
    }
    TCanvas* sector_canvas = new TCanvas("sector", "sector");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    hs->Draw("nostackb");
    std::stringstream fileName;
    fileName << "./adjusted_vop/sector_by_vop/sector-" << sector << "_all.svg";
    sector_canvas->BuildLegend();
    sector_canvas->SaveAs(fileName.str().c_str());
  }

  std::cout << "creating the other histograms" << std::endl;
  for (const auto & p : adjusted_vop_sector_mpv) {
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
    fileName << "./adjusted_vop/vop_by_sector/vop-" << vop << "_all.svg";
    
    vop_canvas->BuildLegend();
    vop_canvas->SaveAs(fileName.str().c_str());
  }

  // scatter plot
  TCanvas* adjusted_canvas = new TCanvas();
  TGraph* adjusted_graph = new TGraph(new_all_vop.size(), &new_all_vop[0], &adjusted_all_mpv[0]);
  adjusted_graph->SetTitle("block mpv vs sipm vop (adjusted)");
  adjusted_graph->GetXaxis()->SetTitle("vop");
  adjusted_graph->GetYaxis()->SetTitle("mpv");
  adjusted_graph->SetMarkerStyle(20);
  adjusted_graph->SetMarkerSize(0.5);
  adjusted_graph->Draw("AP");
  // canvas->Update();
  adjusted_canvas->SaveAs("./adjusted_vop/vop_graphs/vop_mpv_scatter.svg");

  std::vector<double> adjusted_err_x;
  for (int i = 0; i < adjusted_hist_vop.size(); i++) {
    adjusted_err_x.push_back(0);
  }

  // discrete graph
  TCanvas* adjusted_canvas2 = new TCanvas();
  TGraphErrors* adjusted_graph2 = new TGraphErrors(adjusted_hist_vop.size(), &adjusted_hist_vop[0], &adjusted_hist_mpv[0], &adjusted_err_x[0], &adjusted_hist_mpv_err[0]);
  adjusted_graph2->SetTitle("Block MPV vs SiPM VOp (Adjusted)");
  adjusted_graph2->GetXaxis()->SetTitle("vop");
  adjusted_graph2->GetYaxis()->SetTitle("mpv");
  adjusted_graph2->GetXaxis()->SetLimits(68.85, 69.4);
  adjusted_graph2->GetYaxis()->SetLimits(220., 380.);
  adjusted_graph2->SetMarkerStyle(21);
  //TExec *adjusted_ex = new TExec("ex","DrawCol();");
  //adjusted_graph2->GetListOfFunctions()->Add(adjusted_ex);
  adjusted_graph2->Draw("AP");
  adjusted_canvas2->SetGrid();
  adjusted_canvas2->SaveAs("./adjusted_vop/vop_graphs/vop_mpv_discrete.svg");

  // tim's pcb series histogram
  //TFile* tim_outfile = new TFile("./adjusted_vop/vop_graphs/tim_graphs.root","RECREATE");
  TCanvas* adjusted_tim_canvas = new TCanvas();
  // first find all pcb series whcih have current data...
  std::vector<std::string> adjusted_pcb_series_with_data;
  std::map<std::string, double> adjusted_pcb_series_to_vop;
  for (int i = 0; i < adjusted_hist_vop.size(); i++) {
    double vop = adjusted_hist_vop[i];
    if (vop_to_board_series.find(vop) == vop_to_board_series.end()) {
      std::cout << "*unable to find vop " << vop << " in pcb series map" << std::endl;
    } else {
      //std::cout << "vop " << vop << " matches to board " << vop_to_board_series[vop] << std::endl;
      adjusted_pcb_series_with_data.push_back(vop_to_board_series[vop]);
      adjusted_pcb_series_to_vop[vop_to_board_series[vop]] = vop;
    }
  }
  int adjusted_num_pcb_series_with_data = adjusted_pcb_series_with_data.size();
  std::sort(adjusted_pcb_series_with_data.begin(), adjusted_pcb_series_with_data.end());

  TH1D* adjusted_tim_empty_hist = new TH1D("h", "tim_empty", adjusted_num_pcb_series_with_data, 0, adjusted_num_pcb_series_with_data);
  adjusted_tim_empty_hist->SetTitle("Mean Block MPV by SiPM PCB Series (Adjusted)");
  adjusted_tim_empty_hist->GetXaxis()->SetTitle("PCB Series");
  adjusted_tim_empty_hist->GetYaxis()->SetTitle("MPV");
  //adjusted_tim_empty_hist->GetYaxis()->SetRange(200, 400);
  adjusted_tim_empty_hist->SetAxisRange(220., 380., "Y");
  for (int i = 0; i < adjusted_num_pcb_series_with_data; i++) {
    adjusted_tim_empty_hist->GetXaxis()->SetBinLabel(i + 1, adjusted_pcb_series_with_data[i].c_str());
  }
  //adjusted_tim_empty_hist->Write();
  adjusted_tim_empty_hist->Draw();

  std::vector<double> adjusted_tim_x_axis;
  for (int i = 0; i < adjusted_num_pcb_series_with_data; i++) {
    adjusted_tim_x_axis.push_back(i + 0.5);
  }
  std::vector<double> adjusted_tim_x_axis_err;
  for (int i = 0; i < adjusted_num_pcb_series_with_data; i++) { adjusted_tim_x_axis_err.push_back(0); }
  std::vector<double> adjusted_tim_y_axis;
  std::vector<double> adjusted_tim_y_axis_err;
  for (std::string series : adjusted_pcb_series_with_data) {
    double vop = adjusted_pcb_series_to_vop[series];
    std::pair<double, double> p = adjusted_tim_y_axis_data[vop];
    adjusted_tim_y_axis.push_back(p.first);
    adjusted_tim_y_axis_err.push_back(p.second);
    //std::cout << "y axis data for vop " << vop << ": (" << p.first << ", " << p.second << ")" << std::endl; 
  }
  TGraphErrors* adjusted_tim_graph = new TGraphErrors(adjusted_num_pcb_series_with_data, &adjusted_tim_x_axis[0], &adjusted_tim_y_axis[0], &adjusted_tim_x_axis_err[0], &adjusted_tim_y_axis_err[0]);

  adjusted_tim_graph->GetYaxis()->SetTitle("MPV");
  adjusted_tim_graph->GetYaxis()->SetLimits(220., 380.);
  adjusted_tim_graph->GetYaxis()->SetTitle("MPV");
  adjusted_tim_graph->SetMarkerStyle(21);
  //adjusted_tim_graph->Write();
  adjusted_tim_graph->Draw("P SAME");
  adjusted_tim_canvas->SetGrid();
  adjusted_tim_canvas->SaveAs("./adjusted_vop/vop_graphs/vop_mpv_discrete_tim.svg");
}