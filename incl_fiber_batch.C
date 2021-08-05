// Author: Mason Housenga

#include <TROOT.h>
#include <TFile.h>
#include <TColor.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TMarker.h>
#include <TExec.h>
#include <TLatex.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include "csvFile.h"
#include <exception>
#include <cmath>

// mason to-do:
// clean up the multiple mean calculations into one (get rid of old)
// clean up mess of code for making plots for 8/2 meeting
// look into "Error in <TAxis::SetBinLabel>: Illegal bin number: 20"

// PARAMETERS
const int num_bins = 30;
const double x_min = 0;
const double x_max = 600;

const std::vector<int> sectors {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17};
int max_sector = sectors[sectors.size() - 1];

const std::vector<int> interface_boards {0, 1, 2, 3, 4, 5};
const std::vector<Color_t> ib_colors {1, 2, 3, 4, 6, 7};

// 20 colors which are quite visually pleasing and maximally (or close to it) distinct
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

// C++ struct to hold information associated with a block
struct Block {
  bool has_mpv;
  int dbn;
  int fiber_batch;
  double mpv;
  double adj_mpv;
};

// HELPER FUNCTIONS
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

// helps get the new errors when dividing two errorful numbers (such as for sigma/mean)
std::pair<std::vector<double>, std::vector<double>> get_err_div(std::vector<double> x, std::vector<double> dx, std::vector<double> y, std::vector<double> dy) {
  if (x.size() != dx.size() || dx.size() != y.size() || y.size() != dy.size()) {
    return std::make_pair(std::vector<double>({-1}), std::vector<double>({-1}));
  }
  std::vector<double> z;
  std::vector<double> dz;
  for (int i = 0; i < x.size(); i++) {
    double x0 = x[i];
    double dx0 = dx[i];
    double y0 = y[i];
    double dy0 = dy[i];
    double z0 = x0/y0;
    double dz0 = z0*pow(pow(dx0/x0, 2.0) + pow(dy0/y0, 2.0), 0.5);
    z.push_back(z0);
    dz.push_back(dz0);
    //std::cout << "CALCULATING (" << i << "): " << x0 << " ± " << dx0 << " / " << y0 << " ± " << dy0 << " = " << z0 << " ± " << dz0 << std::endl;
  }
  return std::make_pair(z, dz);
}

// generates sector map from database "sectors" sheet
std::map<int, std::vector<std::string>> get_sector_map() {
  std::map<int, std::vector<std::string>> sector_map;
  std::set<std::string> dbns;
  for (int sector : sectors) {
    sector_map[sector] = std::vector<std::string>(96);
  }
  std::fstream sector_map_file;
  sector_map_file.open("db_sectors.csv", std::ios::in);
  std::vector<std::string> row;
  std::string line, word, temp;
  int line_num = 1; // 1-idx
  std::cout << "reading sector map" << std::endl;
  while (std::getline(sector_map_file, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    std::stringstream s(line);
    while (std::getline(s, word, ',')) {
      row.push_back(word);
    }
    if (line_num == 1) {
      // ignore first row
      line_num++;
      continue;
    } else if (line_num > 25) {
      // data has ended
      break;
    } else {
      int sector = 1;
      for (int start = 4; start <= 4 + 8 * 63; start += 8) {
        if (std::find(sectors.begin(), sectors.end(), sector) != sectors.end()) {
          int block_offset = (line_num - 2) * 4;
          if (dbns.find(row[start]) == dbns.end()) {
            sector_map[sector][block_offset] = row[start];
            dbns.insert(row[start]);
          } else {
            // found a duplicate dbn! panic!
            std::stringstream err_msg;
            err_msg << "tried to add dbn " << row[start] << " to sector_map, but map already contains this dbn!";
            throw std::logic_error(err_msg.str());
          }
          if (dbns.find(row[start + 1]) == dbns.end()) {
            sector_map[sector][block_offset + 1] = row[start + 1];
            dbns.insert(row[start + 1]);
          } else {
            // found a duplicate dbn! panic!
            std::stringstream err_msg;
            err_msg << "tried to add dbn " << row[start + 1] << " to sector_map, but map already contains this dbn!";
            throw std::logic_error(err_msg.str());
          }
          if (dbns.find(row[start + 2]) == dbns.end()) {
            sector_map[sector][block_offset + 2] = row[start + 2];
            dbns.insert(row[start + 2]);
          } else {
            // found a duplicate dbn! panic!
            std::stringstream err_msg;
            err_msg << "tried to add dbn " << row[start + 2] << " to sector_map, but map already contains this dbn!";
            throw std::logic_error(err_msg.str());
          }
          if (dbns.find(row[start + 3]) == dbns.end()) {
            sector_map[sector][block_offset + 3] = row[start + 3];
            dbns.insert(row[start + 3]);
          } else {
            // found a duplicate dbn! panic!
            std::stringstream err_msg;
            err_msg << "tried to add dbn " << row[start + 3] << " to sector_map, but map already contains this dbn!";
            throw std::logic_error(err_msg.str());
          }
        }
        sector++;
      }
    }
    line_num++;
  }
  return sector_map;
}

// & indicates "pass by reference" in C++, which is needed in this case
void add_fiber_batch_info(std::string file_name, std::map<int, int> &fb_map) {
  std::fstream blocks_1_12;
  blocks_1_12.open(file_name, std::ios::in);
  std::vector<std::string> row;
  std::string line, word, temp;
  int line_num = 1;
  int this_row_num = -1;
  int offset = 0;
  std::cout << "reading " << file_name << std::endl;
  while (std::getline(blocks_1_12, line)) {
    // std::cout << "...line " << line_num << std::endl;
    row.clear();
    std::stringstream s(line);
    while (std::getline(s, word, ',')) {
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
      bool is_int_dbn = true;
      int int_dbn = 0;
      try {
        int_dbn = std::stoi(dbn);
      } catch (std::exception &e) {
        is_int_dbn = false;
      }
      if (row[9] == "0") {
        // handles the special case of fiber batch 0 which does not contain a "-"
        if (!is_int_dbn) {
          std::cout << "error at DBN " << dbn << ": " << " dbn not castable to int!" << std::endl; 
        } else {
          if (fb_map.find(int_dbn) == fb_map.end()) {
            fb_map[int_dbn] = 0;
          } else {
            // found a duplicate dbn! panic!
            std::stringstream err_msg;
            err_msg << "tried to add dbn " << int_dbn << " to fiber batch map, but map already contains this dbn!";
            throw std::logic_error(err_msg.str());
          }
        }
      } else if (split.size() < 2) {
        std::cout << "skipped fiber batch [" << row[9] << "] for DBN " << dbn << " since it does not contain a '-'" << std::endl;
      } else {
        bool is_int_batch = true;
        int fiber_batch = 0;
        try {
          fiber_batch = std::stoi(split[0]);
        } catch (std::exception &e) {
          is_int_batch = false;
        }
        if (!is_int_batch) {
          std::cout << "error at DBN " << dbn << ": '" << split[0] << "' not castable to int!" << std::endl; 
        } else if (!is_int_dbn) {
          std::cout << "error at DBN " << dbn << ": " << " dbn not castable to int!" << std::endl; 
        } else {
          if (int_dbn == 1021) {
            std::cout << "FOUND THAT BITCH" << std::endl;
          }
          if (fb_map.find(int_dbn) == fb_map.end()) {
            fb_map[int_dbn] = fiber_batch;
          } else {
            // found a duplicate dbn! panic!
            std::stringstream err_msg;
            err_msg << "tried to add dbn " << int_dbn << " to fiber batch map, but map already contains this dbn!";
            throw std::logic_error(err_msg.str());
          }
        }
      }
    }
    line_num++;
  }
}


// THE MACRO
void incl_fiber_batch() {
  // READ SECTOR MAPS
  std::map<int, std::vector<std::string>> string_sector_map = get_sector_map();
  std::map<int, std::vector<int>> int_sector_map;
  for (auto const &p : string_sector_map) {
    int sector = p.first;
    std::vector<std::string> string_dbns = p.second;
    int_sector_map[sector] = std::vector<int>(96, -1);
    for (int i = 0; i < 96; i++) {
      std::string string_dbn = string_dbns[i];
      bool is_int_dbn = true;
      int int_dbn = 0;
      try {
        int_dbn = std::stoi(string_dbn);
      } catch (std::exception &e) {
        is_int_dbn = false;
      }
      if (is_int_dbn) {
        int_sector_map[sector][i] = int_dbn;
      }
    }
  }

  // IMPORT FIBER BATCH INFORMATION
  std::map<int, int> dbn_to_fiber_batch;
  add_fiber_batch_info("blocks_1-12.csv", dbn_to_fiber_batch);
  add_fiber_batch_info("blocks_13-64.csv", dbn_to_fiber_batch);

  // READ MPV DATA
  std::map<int, double> dbn_mpv;
  std::map<int, std::vector<Block>> block_data;
  for (int sector : sectors) {
    std::cout << "SECTOR " << sector << ":" << std::endl;
    std::stringstream sectorFileName;
    sectorFileName << "./sector_data/sector" << sector << ".root";
    TFile* sectorFile = new TFile(sectorFileName.str().c_str());
    TH1D* sector_data;
    std::stringstream blockFileName;
    blockFileName << "h_run" << sector << "_block;1";
    sectorFile->GetObject(blockFileName.str().c_str(), sector_data);
    block_data[sector] = std::vector<Block>(96);
    for (int ib : interface_boards) {
      TH1D* data;
      std::stringstream ibFileName;
      ibFileName << "h_run" << sector << "_block_ib_" << ib << ";1";
      // std::cout << "getting object" << std::endl;
      sectorFile->GetObject(ibFileName.str().c_str(), data);
      int min_block_num = ib * 16;
      int max_block_num = ib * 16 + 15;
      //std::cout << "ib " << ib << " has " <<  data->GetNcells() << " cells" << std::endl;
      for (int i = 1; i < data->GetNbinsX(); i++) {
        // std::cout << "reading bin " << i << std::endl;
        double content = data->GetBinContent(i);
        int block_num = i - 1;
        if (block_num < min_block_num || block_num > max_block_num) {
          // out of range for this ib, omit
          if (block_num < 0 || block_num > 95) {
            std::cout << "block " << block_num << " is out of range" << std::endl;
          }
          continue;
        }
        std::cout << "sector " << sector << ", block " << block_num << " (bin " << i << "); ";
        std::string dbn = string_sector_map[sector][block_num];
        if (dbn == "#N/A") {
          std::cout << "sector map missing value!" << std::endl;
          continue;
        } else if (dbn[0] == 'F') {
          // fudan block, omit (for now)
          std::cout << "rejected fudan block" << std::endl;
          continue;
        }
        // check for int dbn
        bool is_int_dbn = true;
        int int_dbn = 0;
        try {
          int_dbn = std::stoi(dbn);
        } catch (std::exception &e) {
          is_int_dbn = false;
        }
        if (!is_int_dbn) {
          // non-integer dbn, omit for now (casting to into solves issues with importing data from database)
          std::cout << "DBN " << dbn << "; not castable to int!";
          continue;
        }
        std::cout << "DBN " << std::stoi(dbn) << ": ";
        if (content <= 0 || content >= 1000) {
          // that mpv is probably not good data, omit
          std::cout << "rejected bin content " << content << std::endl;
          continue;
        }
        // check fiber batch information
        if (dbn_to_fiber_batch.find(int_dbn) != dbn_to_fiber_batch.end()) {
          int fiber_batch = dbn_to_fiber_batch[int_dbn];
          std::cout << "mpv " << content << std::endl;
          block_data[sector][block_num] = {true, int_dbn, fiber_batch, content, -1.0};
        } else {
          // unable to find fiber batch information for this block
          std::cout << "***DBN " << " does not have fiber batch information!" << std::endl;
          block_data[sector][block_num] = {true, int_dbn, -1, content, -1.0};
        }
        // do this regardless of whether we have fiber batch information or not
        dbn_mpv[int_dbn] = content;
      }
    }
  }

  // CALCULATE FIBER BATCH CORRECTION FACTORS
  int num_good_blocks = 0;
  double mean_mpv = 0;
  std::map<int, double> fiber_batch_avg_mpv;
  std::map<int, int> fiber_batch_counts;
  std::map<int, double> fiber_batch_correction_factors;
  for (auto const &p : block_data) {
    int sector = p.first;
    std::vector<Block> blocks = p.second;
    for (Block block : blocks) {
      if (block.has_mpv && block.fiber_batch >= 0) {
        fiber_batch_avg_mpv[block.fiber_batch] += block.mpv;
        fiber_batch_counts[block.fiber_batch]++;
        mean_mpv += block.mpv;
        num_good_blocks++;
      }
    }
  }
  mean_mpv /= num_good_blocks;
  for (auto const &p : fiber_batch_avg_mpv) {
    int batch = p.first;
    fiber_batch_avg_mpv[batch] = fiber_batch_avg_mpv[batch] / fiber_batch_counts[batch];
    fiber_batch_correction_factors[batch] = mean_mpv / fiber_batch_avg_mpv[batch];
    //std::cout << "fiber batch " << batch << ": avg mpv " << fiber_batch_avg_mpv[batch] << " => factor " << mean_mpv / fiber_batch_avg_mpv[batch] << std::endl;
  }

  // ADD ADJUSTED MPV TO BLOCKS
  for (auto &p : block_data) {
    int sector = p.first;
    std::vector<Block> *blocks = &p.second;
    for (int i = 0; i < blocks->size(); i++) {
      Block *block_ptr = &((*blocks)[i]);
      if (block_ptr->has_mpv && block_ptr->fiber_batch >= 0) {
        block_ptr->adj_mpv = fiber_batch_correction_factors[block_ptr->fiber_batch] * block_ptr->mpv;
      }
    }
  }

  double new_mean_before = 0.0;
  double new_mean_after = 0.0;
  int new_count = 0;
  std::map<int, double> new_sector_means_before;
  std::map<int, double> new_sector_means_after;
  std::map<int, int> new_sector_counts;
  std::map<int, double> fb_means_before;
  std::map<int, double> fb_means_after;
  std::map<int, int> fb_counts;
  std::map<int, std::vector<double>> fb_before_mpvs;
  for (auto const &p : block_data) {
    int sector = p.first;
    new_sector_means_before[sector] = 0.0;
    new_sector_means_after[sector] = 0.0;
    new_sector_counts[sector] = 0;
    std::vector<Block> blocks = p.second;
    for (Block block : blocks) {
      if (block.has_mpv && block.fiber_batch >= 0) {
        new_mean_before += block.mpv;
        new_mean_after += block.adj_mpv;
        new_sector_means_before[sector] += block.mpv;
        new_sector_means_after[sector] += block.adj_mpv;
        new_count++;
        new_sector_counts[sector]++;
        int fiber_batch = block.fiber_batch;
        if (fb_means_before.find(fiber_batch) == fb_means_before.end()) {
          fb_means_before[fiber_batch] = block.mpv;
          fb_means_after[fiber_batch] = block.adj_mpv;
          fb_counts[fiber_batch] = 1;
        } else {
          fb_means_before[fiber_batch] += block.mpv;
          fb_means_after[fiber_batch] += block.adj_mpv;
          fb_counts[fiber_batch]++;
        }
        fb_before_mpvs[fiber_batch].push_back(block.mpv);
      }
    }
    new_sector_means_before[sector] /= new_sector_counts[sector];
    new_sector_means_after[sector] /= new_sector_counts[sector];
  }
  new_mean_before /= new_count;
  new_mean_after /= new_count;
  //std::cout << "overall mean before: " << new_mean_before << "; overall mean after: " << new_mean_after << std::endl;
  for (int sector : sectors) {
    //std::cout << "sector " << sector << ": " << new_sector_means_before[sector] << " -> " << new_sector_means_after[sector] << " (" << new_sector_counts[sector] << " blocks)" << std::endl;
  }
  for (auto const &p : fb_counts) {
    int fiber_batch = p.first;
    int count = p.second;
    fb_means_before[fiber_batch] /= count;
    fb_means_after[fiber_batch] /= count;
    //std::cout << "fiber batch " << fiber_batch << ": " << fb_means_before[fiber_batch] << " -> " << fb_means_after[fiber_batch] << std::endl;
  }

  std::map<int, std::map<int, std::pair<double, int>>> sector_fb_before_means_counts;
  std::map<int, std::map<int, std::vector<double>>> sector_fb_before_mpvs;
  for (auto const &p : block_data) {
    int sector = p.first;
    sector_fb_before_means_counts[sector] = std::map<int, std::pair<double, int>>();
    std::vector<Block> blocks = p.second;
    for (Block block : blocks) {
      if (block.has_mpv && block.fiber_batch >= 0) {
        int fiber_batch = block.fiber_batch;
        if (sector_fb_before_means_counts[sector].find(fiber_batch) == sector_fb_before_means_counts[sector].end()) {
          sector_fb_before_means_counts[sector][fiber_batch] = std::make_pair(block.mpv, 1);
        } else {
          sector_fb_before_means_counts[sector][fiber_batch].first += block.mpv;
          sector_fb_before_means_counts[sector][fiber_batch].second++;
        }
        sector_fb_before_mpvs[sector][fiber_batch].push_back(block.mpv);
      }
    }
    std::cout << "SECTOR " << sector << std::endl;
    for (auto const &q : sector_fb_before_means_counts[sector]) {
      int fiber_batch = q.first;
      sector_fb_before_means_counts[sector][fiber_batch].first /= sector_fb_before_means_counts[sector][fiber_batch].second;
      std::cout << "-> fiber batch " << fiber_batch << ": " << sector_fb_before_means_counts[sector][fiber_batch].first << " [" << sector_fb_before_means_counts[sector][fiber_batch].second << " blocks] (overall " << fb_means_before[fiber_batch] << " [" << fb_counts[fiber_batch] << " blocks] => correction factor = " << new_mean_before / fb_means_before[fiber_batch] << ")" << std::endl;
    }
  }

  // BLOCK MPV HISTS FOR EACH SECTOR FOR EACH FIBER BATCH
  const int num_bins = 30;
  const double x_min = 0;
  const double x_max = 600;
  for (int sector : sectors) {
    for (auto const &p : sector_fb_before_mpvs[sector]) {
      int fiber_batch = p.first;
      std::vector<double> this_sector_fb_mpvs = p.second;
      std::vector<double> total_fb_mpvs = fb_before_mpvs[fiber_batch];
      std::stringstream hs_title;
      hs_title << "Sector " << sector << " FB " << fiber_batch << ";MPV;Num Blocks";
      THStack* hs = new THStack("hs", hs_title.str().c_str());
      std::stringstream this_sector_fb_hist_title;
      this_sector_fb_hist_title << "sector " << sector << " fb " << fiber_batch << ";MPV;Num Blocks";
      TH1D* sector_fb_hist = new TH1D("h", this_sector_fb_hist_title.str().c_str(), num_bins, x_min, x_max);
      std::stringstream total_fb_hist_title;
      total_fb_hist_title << "fb " << fiber_batch << ";MPV;Num Blocks";
      TH1D* total_fb_hist = new TH1D("h", total_fb_hist_title.str().c_str(), num_bins, x_min, x_max);
      for (double mpv : this_sector_fb_mpvs) {
        sector_fb_hist->Fill(mpv);
      }
      for (double mpv : total_fb_mpvs) {
        total_fb_hist->Fill(mpv);
      }
      sector_fb_hist->SetFillColorAlpha(kBlue, 0.7);
      total_fb_hist->SetFillColorAlpha(kBlack, 0.7);
      hs->Add(sector_fb_hist);
      hs->Add(total_fb_hist);
      std::stringstream hs_filename;
      hs_filename << "./plots/sector_fb/sector" << sector << "_fb" << fiber_batch << ".svg";
      TCanvas* hs_canvas = new TCanvas();
      //hs->SetMinimum(0.);
      //hs->SetMaximum(600.);
      hs->Draw("nostackb");
      hs_canvas->SetGrid();
      hs_canvas->Update();
      double sector_fb_mean = sector_fb_before_means_counts[sector][fiber_batch].first;
      TLine* sector_fb_mean_line = new TLine(sector_fb_mean, 0.0, sector_fb_mean, hs_canvas->GetUymax());
      sector_fb_mean_line->SetLineColor(kBlue);
      sector_fb_mean_line->SetLineStyle(2);
      sector_fb_mean_line->SetLineWidth(2);
      sector_fb_mean_line->Draw();
      double total_fb_mean = fb_means_before[fiber_batch];
      TLine* total_fb_mean_line = new TLine(total_fb_mean, 0.0, total_fb_mean, hs_canvas->GetUymax());
      std::cout << "ymax: " << hs_canvas->GetUymax() << std::endl;
      total_fb_mean_line->SetLineColor(kBlack);
      total_fb_mean_line->SetLineStyle(2);
      total_fb_mean_line->SetLineWidth(2);
      total_fb_mean_line->Draw();
      TLegend* hs_legend = new TLegend(0.7, 0.65, 0.9, 0.9);
      hs_legend->AddEntry(sector_fb_hist, "Sector", "f");
      hs_legend->AddEntry(total_fb_hist, "Batch", "f");
      hs_legend->AddEntry(sector_fb_mean_line, Form("Sector Avg (%.1f)", sector_fb_mean), "l");
      hs_legend->AddEntry(total_fb_mean_line, Form("Batch Avg (%.1f)", total_fb_mean), "l");
      hs_legend->Draw();
      hs_canvas->SaveAs(hs_filename.str().c_str());
    }
  }

  // BLOCK CALIBRATION FACTOR PLOTS FOR EACH SECTOR
  std::map<int, std::map<int, Color_t>> sector_fb_colors;
  std::map<int, std::vector<int>> sector_batches;
  std::map<int, std::vector<std::string>> sector_bin_labels;
  std::map<int, std::map<int, TH1D*>> sector_batch_hists;
  for (int sector : sectors) {
    // first determine the number of batches in this sector
    std::vector<int> batches;
    std::map<int, Color_t> fb_colors;
    for (int block_num = 0; block_num < 96; block_num++) {
      Block block = block_data[sector][block_num];
      if (block.has_mpv && block.fiber_batch >= 0) {
        if (std::find(batches.begin(), batches.end(), block.fiber_batch) == batches.end()) {
          // add this fb
          batches.push_back(block.fiber_batch);
          // create new color for this fb
          Color_t ci = TColor::GetFreeColorIndex();
          int nth_kelly = fb_colors.size();
          TColor* color = new TColor(ci, kelly_colors[nth_kelly][0], kelly_colors[nth_kelly][1], kelly_colors[nth_kelly][2]);
          fb_colors[block.fiber_batch] = ci;
          //std::cout << "sector " << sector << " fb " << block.fiber_batch << " has kelly color " << nth_kelly << std::endl;
        }
      }
    }
    //std::cout << "SECTOR " << sector << " HAS " << batches.size() << " BATCHES" << std::endl;
    std::sort(batches.begin(), batches.end());
    sector_batches[sector] = batches;
    sector_fb_colors[sector] = fb_colors;
    // generate bin labels
    std::vector<std::string> bin_labels;
    for (int batch : batches) {
      double correction_factor = new_mean_before / fb_means_before[batch];
      std::stringstream label;
      label << "FB " << batch << ": " << Form("%.3f", correction_factor);
      bin_labels.push_back(label.str());
    }
    sector_bin_labels[sector] = bin_labels;
    for (int batch : batches) {
      double correction_factor = new_mean_before / fb_means_before[batch];
      std::stringstream label;
      label << "FB " << batch << ": " << Form("%.3f", correction_factor);
      bin_labels.push_back(label.str());
    }
    std::stringstream hs_title;
    hs_title << "Sector " << sector << ";Fiber Batch;Num Blocks";
    THStack* hs = new THStack("hs", hs_title.str().c_str());
    std::map<int, TH1D*> batch_hists;
    for (int block_num = 0; block_num < 96; block_num++) {
      Block block = block_data[sector][block_num];
      if (block.has_mpv && block.fiber_batch >= 0) {
        if (batch_hists.find(block.fiber_batch) == batch_hists.end()) {
          // create hist for this fb
          std::stringstream batch_hist_title;
          batch_hist_title << "fb " << block.fiber_batch << ";Fiber Batch;Num Blocks";
          int num_batches = batches.size();
          TH1D* batch_hist = new TH1D("h", batch_hist_title.str().c_str(), num_batches, 0, num_batches);
          batch_hist->GetSumw2();
          batch_hist->SetLineColor(fb_colors[block.fiber_batch]);
          batch_hist->SetFillColor(fb_colors[block.fiber_batch]);
          //batch_hist->GetXaxis()->LabelsOption("v");
          for (int i = 0; i < bin_labels.size(); i++) {
            batch_hist->GetXaxis()->SetBinLabel(i + 1, std::to_string(batches[i]).c_str());
            //batch_hist->GetXaxis()->SetBinLabel(i + 1, bin_labels[i].c_str());
            //batch_hist->GetXaxis()->SetBinLabel(i + 1, "");
          }
          batch_hists[block.fiber_batch] = batch_hist;
        }
        batch_hists[block.fiber_batch]->Fill(std::find(batches.begin(), batches.end(), block.fiber_batch) - batches.begin());
      }
    }
    sector_batch_hists[sector] = batch_hists;
    for (auto const &p : batch_hists) {
      hs->Add(p.second, "B");
    }
    TCanvas* sector_canvas = new TCanvas();
    sector_canvas->SetGrid();
    //sector_canvas->SetRightMargin(0.3);
    hs->Draw();
    /*
    TLegend* legend = new TLegend(0.725, 0.025, 0.975, 0.975);
    for (int i = 0; i < batches.size(); i++) {
      legend->AddEntry(batch_hists[batches[i]], bin_labels[i].c_str());
    }
    legend->Draw();
    */
    std::stringstream hs_filename;
    hs_filename << "./plots/sector_fb_hists/sector" << sector << ".svg";
    sector_canvas->SaveAs(hs_filename.str().c_str());
  }

  //
  for (int sector : sectors) {
    std::stringstream hs_title;
    hs_title << "Sector " << sector << ";Block Number;Correction Factor";
    THStack* hs = new THStack("hs", hs_title.str().c_str());

    std::map<int, std::pair<TH1D*, TH1D*>> batch_hists;

    for (int block_num = 0; block_num < 96; block_num++) {
      Block block = block_data[sector][block_num];
      if (block.has_mpv && block.fiber_batch >= 0) {
        if (batch_hists.find(block.fiber_batch) == batch_hists.end()) {
          // create hist for this fb
          std::stringstream batch_hist_up_title;
          batch_hist_up_title << "fb " << block.fiber_batch << " up;Block Number;Correction Factor";
          TH1D* batch_hist_up = new TH1D("h_up", batch_hist_up_title.str().c_str(), 96, 0, 96);
          batch_hist_up->GetSumw2();
          batch_hist_up->SetMarkerColor(sector_fb_colors[sector][block.fiber_batch]);
          batch_hist_up->SetMarkerStyle(22); // up triangle
          // and the down hist
          std::stringstream batch_hist_down_title;
          batch_hist_down_title << "fb " << block.fiber_batch << " down;Block Number;Correction Factor";
          TH1D* batch_hist_down = (TH1D*) batch_hist_up->Clone("h_down");
          batch_hist_down->SetTitle(batch_hist_down_title.str().c_str());
          batch_hist_down->SetMarkerStyle(23); // down triangle
          // add to the map
          batch_hists[block.fiber_batch] = std::make_pair(batch_hist_up, batch_hist_down);
          //std::cout << "CREATED HISTS FOR FB " << fiber_batch << std::endl;
        }
        double correction_factor = new_mean_before / fb_means_before[block.fiber_batch];
        if (correction_factor >= 1) {
          batch_hists[block.fiber_batch].first->SetBinContent(block_num + 1, correction_factor);
          //batch_hists[fiber_batch].second->SetBinContent(block_num + 1, -9000);
          //std::cout << "added an up arrow at block " << block_num << " fb " << fiber_batch << " factor " << correction_factor << std::endl;
        } else {
          //batch_hists[fiber_batch].first->SetBinContent(block_num + 1, -9000);
          batch_hists[block.fiber_batch].second->SetBinContent(block_num + 1, correction_factor);
          //std::cout << "added a down arrow at block " << block_num << " fb " << fiber_batch << " factor " << correction_factor << std::endl;
        }
      }
    }
    for (auto const &p : batch_hists) {
      hs->Add(p.second.first, "P");
      hs->Add(p.second.second, "P");
    }
    TCanvas* sector_canvas = new TCanvas();
    sector_canvas->SetGrid();
    hs->SetMinimum(0.75);
    hs->SetMaximum(1.25);
    hs->Draw("NOSTACK");
    TLine* one_line = new TLine(0.0, 1.0, 96.0, 1.0);
    one_line->SetLineStyle(2);
    one_line->SetLineWidth(2);
    one_line->Draw();
    std::stringstream hs_filename;
    hs_filename << "./plots/sector_correction_factors/sector" << sector << ".svg";
    sector_canvas->SaveAs(hs_filename.str().c_str());
  }

  // make combined tim graphs, colored by fiber batch instead of interface board
  for (int sector : sectors) {
    std::stringstream hs_before_title;
    hs_before_title << "Sector " << sector << " (Before);Block Number;MPV";
    THStack* hs_before = new THStack("hs", hs_before_title.str().c_str());

    std::stringstream hs_after_title;
    hs_after_title << "Sector " << sector << " (After);Block Number;MPV";
    THStack* hs_after = new THStack("hs", hs_after_title.str().c_str());

    std::map<int, std::pair<TH1D*, TH1D*>> batch_hists;
    for (int block_num = 0; block_num < 96; block_num++) {
      Block block = block_data[sector][block_num];
      if (block.has_mpv && block.fiber_batch >= 0) {
        if (batch_hists.find(block.fiber_batch) == batch_hists.end()) {
          // original hist
          std::stringstream batch_hist_orig_title;
          batch_hist_orig_title << "fb " << block.fiber_batch << " orig up;Block Number;MPV";
          TH1D* batch_hist_orig = new TH1D("h_up", batch_hist_orig_title.str().c_str(), 96, 0, 96);
          batch_hist_orig->GetSumw2();
          batch_hist_orig->SetMarkerStyle(20); // solid circle
          //batch_hist_orig->SetMarkerColorAlpha(sector_fb_colors[sector][block.fiber_batch], 0.5);
          batch_hist_orig->SetMarkerColor(kBlack);
          
          // adj hist
          std::stringstream batch_hist_adj_title;
          batch_hist_adj_title << "fb " << block.fiber_batch << " adj;Block Number;MPV";
          TH1D* batch_hist_adj = (TH1D*) batch_hist_orig->Clone("h_adj");
          batch_hist_adj->SetTitle(batch_hist_adj_title.str().c_str());
          batch_hist_adj->SetMarkerStyle(20); // solid circle
          //batch_hist_adj->SetMarkerColor(sector_fb_colors[sector][block.fiber_batch]);
          //batch_hist_adj->SetMarkerColor(sector_fb_colors[sector][block.fiber_batch]);
          batch_hist_adj->SetMarkerColor(kBlack);
          // add to the map
          batch_hists[block.fiber_batch] = std::make_pair(batch_hist_orig, batch_hist_adj);
          //std::cout << "CREATED HISTS FOR FB " << block.fiber_batch << std::endl;
        }
        //std::cout << "ADJ MPV IS " << adj_mpv << std::endl;
        batch_hists[block.fiber_batch].first->SetBinContent(block_num + 1, block.mpv);
        batch_hists[block.fiber_batch].second->SetBinContent(block_num + 1, block.adj_mpv);
        //double correction_factor = new_mean_before / fb_means_before[block.fiber_batch];
      }
    }
    for (auto const &p : batch_hists) {
      hs_before->Add(p.second.first, "P");
      hs_after->Add(p.second.second, "P");
    }

    double sector_before_mean = new_sector_means_before[sector];
    double sector_after_mean = new_sector_means_after[sector];

    TCanvas* before_canvas = new TCanvas();
    hs_before->Draw();
    before_canvas->SetGrid();
    hs_before->SetMinimum(0.0);
    hs_before->SetMaximum(600.0);

    // draw sector mean lines
    TLine* bef_orig_avg_sector_mean_line = new TLine(0.0, sector_before_mean, 96.0, sector_before_mean);
    bef_orig_avg_sector_mean_line->SetLineColor(kRed);
    bef_orig_avg_sector_mean_line->SetLineStyle(2);
    bef_orig_avg_sector_mean_line->SetLineWidth(1);
    bef_orig_avg_sector_mean_line->Draw();
    TLine* bef_adj_avg_sector_mean_line = new TLine(0.0, new_mean_before, 96.0, new_mean_before);
    bef_adj_avg_sector_mean_line->SetLineColor(kBlue);
    bef_adj_avg_sector_mean_line->SetLineStyle(2);
    bef_adj_avg_sector_mean_line->SetLineWidth(1);
    bef_adj_avg_sector_mean_line->Draw();

    TLegend* bef_hs_legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    bef_hs_legend->AddEntry(bef_orig_avg_sector_mean_line, Form("Sector Avg (%.1f)", sector_before_mean), "l");
    bef_hs_legend->AddEntry(bef_adj_avg_sector_mean_line, Form("Overall Avg (%.1f)", new_mean_before), "l");
    bef_hs_legend->Draw();

    std::stringstream hs_before_filename;
    hs_before_filename << "./plots/sector_before_after/sector" << sector << "before.svg";
    hs_before_filename << "./plots/sector_before_after/sector" << sector << "before.png";
    before_canvas->SaveAs(hs_before_filename.str().c_str());

    TCanvas* after_canvas = new TCanvas();
    hs_after->Draw();
    after_canvas->SetGrid();
    hs_after->SetMinimum(0.0);
    hs_after->SetMaximum(600.0);

    // draw sector mean lines
    TLine* aft_orig_avg_sector_mean_line = new TLine(0.0, sector_after_mean, 96.0, sector_after_mean);
    aft_orig_avg_sector_mean_line->SetLineColor(kRed);
    aft_orig_avg_sector_mean_line->SetLineStyle(2);
    aft_orig_avg_sector_mean_line->SetLineWidth(1);
    aft_orig_avg_sector_mean_line->Draw();
    TLine* aft_adj_avg_sector_mean_line = new TLine(0.0, new_mean_before, 96.0, new_mean_before);
    aft_adj_avg_sector_mean_line->SetLineColor(kBlue);
    aft_adj_avg_sector_mean_line->SetLineStyle(2);
    aft_adj_avg_sector_mean_line->SetLineWidth(1);
    aft_adj_avg_sector_mean_line->Draw();

    TLegend* aft_hs_legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    aft_hs_legend->AddEntry(aft_orig_avg_sector_mean_line, Form("Sector Avg (%.1f)", sector_after_mean), "l");
    aft_hs_legend->AddEntry(aft_adj_avg_sector_mean_line, Form("Overall Avg (%.1f)", new_mean_before), "l");
    aft_hs_legend->Draw();

    std::stringstream hs_after_filename;
    hs_after_filename << "./plots/sector_before_after/sector" << sector << "after.svg";
    hs_after_filename << "./plots/sector_before_after/sector" << sector << "after.png";
    after_canvas->SaveAs(hs_after_filename.str().c_str());


    //TCanvas* sector_canvas = new TCanvas();
    //sector_canvas->SetGrid();
    //sector_canvas->SetRightMargin(0.3);
    //hs->SetMinimum(0.0);
    //hs->SetMaximum(600.0);
    //hs->Draw("NOSTACK");
    //std::stringstream hs_filename;
    //hs_filename << "./proof/combined/sector" << sector << ".svg";
    //sector_canvas->SaveAs(hs_filename.str().c_str());

  }

  // MAKE FIBER BATCH CORRECTION FACTOR LEGENDS FOR EACH SECTOR
  for (int sector : sectors) {
    TCanvas* correction_factor_canvas = new TCanvas();
    TLegend* legend = new TLegend(0.0, 0.0, 1.0, 1.0);
    for (int i = 0; i < sector_batches[sector].size(); i++) {
      legend->AddEntry(sector_batch_hists[sector][sector_batches[sector][i]], sector_bin_labels[sector][i].c_str());
    }
    legend->Draw();
    std::stringstream hs_filename;
    hs_filename << "./plots/sector_fb_legends/sector" << sector << ".svg";
    correction_factor_canvas->SaveAs(hs_filename.str().c_str());
  }
}