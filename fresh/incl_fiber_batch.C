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
  // <has_mpv, dbn, fiber_batch, mpv, adj_mpv>
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
    int num_blocks = 0;
    double original_mean_mpv = 0;
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
          num_blocks++;
          original_mean_mpv += content;
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
  for (auto const &p : block_data) {
    int sector = p.first;
    std::vector<Block> blocks = p.second;
    for (Block block : blocks) {
      if (block.has_mpv && block.fiber_batch >= 0) {
        block.adj_mpv = fiber_batch_correction_factors[block.fiber_batch] * block.mpv;
      }
    }
  }
}
