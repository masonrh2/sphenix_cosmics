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
#include <vector>
#include <map>
#include <algorithm>
#include "csvFile.h"
#include <exception>
#include <cmath>

// the root of all evil

const int num_bins = 30;
const double x_min = 0;
const double x_max = 600;

const std::vector<int> sectors {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17};
int max_sector = sectors[sectors.size() - 1];
const std::vector<int> interface_boards {0, 1, 2, 3, 4, 5};

const std::vector<Color_t> ib_colors {1, 2, 3, 4, 6, 7};

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

std::pair<std::vector<double>, std::vector<double>> get_err_div(std::vector<double> x, std::vector<double> dx, std::vector<double> y, std::vector<double> dy) {
  if (x.size() != dx.size() || dx.size() != y.size() || y.size() != dy.size()) {
    return std::make_pair(std::vector<double>({-1}), std::vector<double>({-1}));
  }
  std::vector<double> z;
  std::vector<double> dz;
  for (int i = 0; i < x.size(); i++) {
    double x0 = x[0];
    double dx0 = dx[0];
    double y0 = y[0];
    double dy0 = dy[0];
    double z0 = x0/y0;
    double dz0 = z0*pow(pow(dx0/x0, 2.0) + pow(dy0/y0, 2.0), 0.5);
    z.push_back(z0);
    dz.push_back(dz0);
  }
  return std::make_pair(z, dz);
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

  // READ SECTOR MAPS
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
  
  std::map<int, int> dbn_to_fiber_batch;

  // READ DATABASE
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
      if (dbn == "134") {
        std::cout << "***FOUND THE THING, FBN [" << row[9] << "]" << std::endl;
        bool is_zero = row[9] == "0";
        std::cout << "   AND the answer is " << is_zero << std::endl;
      }
      if (row[9] == "0") {
        // handles the special case of fiber batch 0 which does not contain a "-"
        bool is_int_dbn = true;
        int int_dbn = 0;
        try {
          int_dbn = std::stoi(dbn);
        } catch (std::exception &e) {
          is_int_dbn = false;
        }
        if (!is_int_dbn) {
          std::cout << "error at DBN " << dbn << ": " << " dbn not castable to int!" << std::endl; 
        } else {
          dbn_to_fiber_batch[int_dbn] = 0;
          //std::cout << "good at DBN [" << std::stoi(dbn) << "]: FB 0" << std::endl;
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
        bool is_int_dbn = true;
        int int_dbn = 0;
        try {
          int_dbn = std::stoi(dbn);
        } catch (std::exception &e) {
          is_int_dbn = false;
        }
        if (!is_int_batch) {
          std::cout << "error at DBN " << dbn << ": '" << split[0] << "' not castable to int!" << std::endl; 
        } else if (!is_int_dbn) {
          std::cout << "error at DBN " << dbn << ": " << " dbn not castable to int!" << std::endl; 
        } else {
          dbn_to_fiber_batch[int_dbn] = fiber_batch;
          //std::cout << "good at DBN [" << std::stoi(dbn) << "]: FB " << fiber_batch << std::endl;
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
        std::cout << "skipped fiber batch [" << row[9] << "] for DBN " << dbn << " since it does not contain a '-'" << std::endl;
      } else {
        bool is_int_batch = true;
        int fiber_batch = 0;
        try {
          fiber_batch = std::stoi(split[0]);
        } catch (std::exception &e) {
          is_int_batch = false;
        }
        bool is_int_dbn = true;
        int int_dbn = 0;
        try {
          int_dbn = std::stoi(dbn);
        } catch (std::exception &e) {
          is_int_dbn = false;
        }
        if (!is_int_batch) {
          std::cout << "error at DBN " << dbn << ": '" << split[0] << "' not castable to int!" << std::endl; 
        } else if (!is_int_dbn) {
          std::cout << "error at DBN " << dbn << ": " << " dbn not castable to int!" << std::endl; 
        } else {
          dbn_to_fiber_batch[int_dbn] = fiber_batch;
          //std::cout << "good at DBN [" << std::stoi(dbn) << "]: FB " << fiber_batch << std::endl;
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

  // READ FIBER BATCH CORRECTION FACTORS
  /*
  std::map<int, double> fiber_batch_to_scale_factor;
  std::fstream fiber_batch_map;
  fiber_batch_map.open("fiberbatchmeans.csv", std::ios::in);
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
  */

  // READ NEW VOP DATA
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
          //std::cout << "sector " << new_sector_nums[idx] << " has new vop data" << std::endl;
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

  csvfile sector_csv("./new_vop/sector/fit_parameters.csv");
  sector_csv << "sector" << "mean" << "mean err" << "sigma" << "sigma err" << csvfile::endrow;
  csvfile adj_sector_csv("./adjusted_vop/sector/fit_parameters.csv");
  adj_sector_csv << "sector" << "mean" << "mean err" << "sigma" << "sigma err" << csvfile::endrow;

  // READ MPV DATA AND MAKE ALL THE PLOTS
  std::vector<double> all_mpv;
  std::map<int, double> dbn_mpv;
  std::vector<double> new_all_vop;
  std::map<int, std::map<double, std::vector<double>>> new_sector_vop_mpv;
  std::map<double, std::map<int, std::vector<double>>> new_vop_sector_mpv;
  std::map<double, std::vector<double>> new_histogram_data;
  std::vector<double> new_hist_vop;
  std::vector<double> new_hist_mpv_err;
  std::vector<double> new_hist_mpv;
  
  std::vector<double> adjusted_all_mpv;
  std::map<int, double> adjusted_dbn_mpv;
  //std::vector<double> new_all_vop;
  std::map<int, std::map<double, std::vector<double>>> adjusted_sector_vop_mpv;
  std::map<double, std::map<int, std::vector<double>>> adjusted_vop_sector_mpv;
  std::map<double, std::vector<double>> adjusted_histogram_data;
  std::vector<double> adjusted_hist_vop;
  std::vector<double> adjusted_hist_mpv_err;
  std::vector<double> adjusted_hist_mpv;

  // create sector diagrams
  std::vector<double> original_means;
  std::vector<double> original_mean_errs;
  std::vector<double> original_sigmas;
  std::vector<double> original_sigma_errs;
  std::vector<double> adj_means;
  std::vector<double> adj_mean_errs;
  std::vector<double> adj_sigmas;
  std::vector<double> adj_sigma_errs;

  int total_num_blocks = 0;
  double total_original_mean = 0;

  std::map<int, std::vector<std::tuple<bool, int, int, double, double>>> tuple_data;
  std::map<int, double> original_sector_means;

  gStyle->SetOptStat(1000000001); // just header
  gStyle->SetOptFit(1); // default

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

    tuple_data[sector] = std::vector<std::tuple<bool, int, int, double, double>>(96);

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
        std::string dbn = sector_map[sector][block_num];
        if (sector_map[sector][block_num][0] == 'F' || std::stoi(sector_map[sector][block_num]) >= 10000) {
          // fudan block, omit for now
          std::cout << "rejected fudan block" << std::endl;
          continue;
        }

        bool is_int_dbn = true;
        int int_dbn = 0;
        try {
          int_dbn = std::stoi(dbn);
        } catch (std::exception &e) {
          is_int_dbn = false;
        }
        if (!is_int_dbn) {
          // non-integer dbn, omit for now (casting to into solves issues with importing data from database)
          std::cout << "DBN " << std::stoi(dbn) << "; not castable to int!";
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
          //double correction_factor = fiber_batch_to_scale_factor[fiber_batch];
          //double adjusted_content = content * correction_factor;
          //std::cout << "DBN " << std::stoi(dbn) << ": fiber batch " << fiber_batch << "; correction factor " << correction_factor
          //  << " (" << content << "->" << adjusted_content << ")" << std::endl;
          // i guess this is valid data?
          int fiber_batch = dbn_to_fiber_batch[int_dbn];
          double vop = new_sipm_map[sector][block_num];

          num_blocks++;
          original_mean_mpv += content;

          total_num_blocks++;
          total_original_mean += content;

          std::cout << "mpv " << content << std::endl;
          //std::cout << "block " << block_num << " (DBN " << std::stoi(sector_map[sector][block_num]) << "): good data (" << vop << ", " << content  << ")" << std::endl;

          tuple_data[sector][block_num] = std::make_tuple(true, int_dbn, fiber_batch, content, -1.0);
        } else {
          // unable to find fiber batch information for this block
          std::cout << "***DBN " << " does not have fiber batch information!" << std::endl;
          tuple_data[sector][block_num] = std::make_tuple(true, int_dbn, -1, content, -1.0);
        }

        // do this regardless of whether we have fiber batch information or not
        double vop = new_sipm_map[sector][block_num];
        new_all_vop.push_back(vop);
        all_mpv.push_back(content);
        new_histogram_data[vop].push_back(content);
        dbn_mpv[std::stoi(dbn)] = content;
        new_sector_vop_mpv[sector][vop].push_back(content);
        new_vop_sector_mpv[vop][sector].push_back(content);
      }
    }
    original_mean_mpv = original_mean_mpv / num_blocks;
    original_sector_means[sector] = original_mean_mpv;
  }

  // GOT ALL MPV DATA
  total_original_mean = total_original_mean / total_num_blocks;

  // calculate MY fiber batch corrections...
  std::map<int, double> fiber_batch_avg_mpv;
  std::map<int, int> fiber_batch_counts;
  for (auto const &p : dbn_to_fiber_batch) {
    int dbn = p.first;
    int batch = p.second;
    if (dbn_mpv.find(dbn) == dbn_mpv.end()) {
      // unable to find that dbn in dbn_mpv
      //std::cout << "unable to find dbn " << dbn << " in dbn_mpv" << std::endl;
    } else {
      double mpv = dbn_mpv[dbn];
      fiber_batch_avg_mpv[batch] += mpv;
      fiber_batch_counts[batch]++;
      //std::cout << "dbn " << dbn << " in fb " << batch << " has mpv " << mpv << std::endl;
    }
  }
  for (auto const &p : fiber_batch_avg_mpv) {
    int batch = p.first;
    fiber_batch_avg_mpv[batch] = fiber_batch_avg_mpv[batch] / fiber_batch_counts[batch];
    //std::cout << "fiber batch " << batch << ": avg mpv " << fiber_batch_avg_mpv[batch] << " => factor " << total_original_mean / fiber_batch_avg_mpv[batch] << std::endl;
  }

  // make modifications to tims colorful histograms
  std::map<int, double> adjusted_sector_means;
  double max_diff = 0;
  double max_ratio_diff = 0;
  double total_adjusted_mean = 0;
  for (int sector : sectors) {
    std::stringstream sectorFileName;
    sectorFileName << "./sector_data/sector" << sector << ".root";
    TFile* sectorFile = new TFile(sectorFileName.str().c_str());

    TH1D* sector_data;
    std::stringstream blockFileName;
    blockFileName << "h_run" << sector << "_block;1";
    sectorFile->GetObject(blockFileName.str().c_str(), sector_data);

    std::stringstream sector_title;
    sector_title << "Sector " << sector << " (Original);Block Number;MPV";
    THStack* hs = new THStack("hs", sector_title.str().c_str());
    std::stringstream adj_sector_title;
    adj_sector_title << "Sector " << sector << " (Adjusted);Block Number;MPV";
    THStack* adj_hs = new THStack("adj_hs", adj_sector_title.str().c_str());
    std::stringstream combined_sector_title;
    combined_sector_title << "Sector " << sector << " (Combined);Block Number;MPV";
    THStack* combined_hs = new THStack("combined_hs", combined_sector_title.str().c_str());
    std::stringstream diff_sector_title;
    diff_sector_title << "Sector " << sector << " (Difference);Block Number;Δ MPV";
    THStack* diff_hs = new THStack("diff_hs", diff_sector_title.str().c_str());

    // gaussian fit histograms
    std::stringstream sector_hist_title;
    sector_hist_title << "sector " << sector << ";MPV;Num Blocks";
    TH1D* sector_hist = new TH1D("h", sector_hist_title.str().c_str(), num_bins, x_min, x_max);
    sector_hist->GetSumw2();
    std::stringstream adj_sector_hist_title;
    adj_sector_hist_title << "sector " << sector << " (adjusted);MPV;Num Blocks";
    TH1D* adj_sector_hist = new TH1D("adj_h", adj_sector_hist_title.str().c_str(), num_bins, x_min, x_max);
    adj_sector_hist->GetSumw2();

    int num_blocks = 0;
    double adj_mean_mpv = 0;

    for (int ib : interface_boards) {
      TH1D* data;
      std::stringstream ibFileName;
      ibFileName << "h_run" << sector << "_block_ib_" << ib << ";1";
      sectorFile->GetObject(ibFileName.str().c_str(), data);
      TH1D* adj_data = (TH1D*) data->Clone("adj_data");
      TH1D* diff_data = (TH1D*) data->Clone("diff_data");
      int min_block_num = ib * 16;
      int max_block_num = ib * 16 + 15;
      for (int i = min_block_num; i <= max_block_num; i++) {
        auto block_info = tuple_data[sector][i];
        bool has_mpv_data = std::get<0>(block_info);
        int dbn = std::get<1>(block_info);
        int fiber_batch = std::get<2>(block_info);
        double mpv = std::get<3>(block_info);
        
        if (has_mpv_data && fiber_batch >= 0) {
          // has fiber batch correction data
          if (fiber_batch_avg_mpv[fiber_batch] == 0.0) { std::cout << "***PANIC, fiber batch " << fiber_batch << " has zero avg mpv (DBN " << dbn << ")" << std::endl; }
          double correction_factor = total_original_mean / fiber_batch_avg_mpv[fiber_batch];
          //double adjusted_content = content * correction_factor;
          double adj_mpv = mpv * correction_factor;
          //std::get<4>(block_info) = adj_mpv; // DOES NOT WORK
          std::get<4>(tuple_data[sector][i]) = adj_mpv;
          //std::cout << "*****TRIED TO SET ADJ MPV FOR DBN " << dbn << " TO " << adj_mpv << ": " << std::get<4>(block_info) << std::endl;
          adj_mean_mpv += adj_mpv;
          total_adjusted_mean += adj_mpv;

          data->SetBinContent(i, mpv);
          adj_data->SetBinContent(i, adj_mpv);
          diff_data->SetBinContent(i, adj_mpv - mpv);
          sector_hist->Fill(mpv);
          adj_sector_hist->Fill(adj_mpv);

          double vop = new_sipm_map[sector][i];

          adjusted_all_mpv.push_back(adj_mpv);
          adjusted_histogram_data[vop].push_back(adj_mpv);
          adjusted_dbn_mpv[dbn] = adj_mpv;
          adjusted_sector_vop_mpv[sector][vop].push_back(adj_mpv);
          adjusted_vop_sector_mpv[vop][sector].push_back(adj_mpv);

          if (std::abs(adj_mpv - mpv) > max_diff) { max_diff = std::abs(adj_mpv - mpv); }
          if (std::abs(correction_factor - 1) > max_ratio_diff) { max_ratio_diff = std::abs(correction_factor - 1); }
        } else if (has_mpv_data) {
          std::cout << "***WARNING, dbn " << dbn << " has MPV data but no fiber batch!" << std::endl;
        } else {
          // omit this data
          data->SetBinContent(i, -9000.);
          adj_data->SetBinContent(i, -9000.);
          diff_data->SetBinContent(i, -9000.);
        }
      }
      data->GetXaxis()->SetLimits(0.0, 95.0);
      adj_data->GetXaxis()->SetLimits(0.0, 95.0);
      diff_data->GetXaxis()->SetLimits(0.0, 95.0);

      hs->Add(data);
      adj_hs->Add(adj_data);
      diff_hs->Add(diff_data);

      data->SetMarkerStyle(21); // sqruares for tim
      adj_data->SetMarkerStyle(21); // sqruares for tim
      combined_hs->Add(data);
      combined_hs->Add(adj_data);
    }

    adj_mean_mpv = adj_mean_mpv / num_blocks;
    adjusted_sector_means[sector] = adj_mean_mpv;

    std::stringstream sector_filename;
    TCanvas* sector_canvas = new TCanvas();
    sector_filename << "./sector_diagrams/original/sector" << sector << ".svg";
    hs->SetMinimum(0.);
    hs->SetMaximum(600.);
    hs->Draw("nostack");
    sector_canvas->SetGrid();
    double original_mean = original_sector_means[sector];
    TLine* original_mean_line = new TLine(0.0, original_mean, 95.0, original_mean);
    //original_mean_line->SetLineColor(kRed);
    original_mean_line->SetLineStyle(2);
    original_mean_line->SetLineWidth(2);
    original_mean_line->Draw();
    TPaveText* original_pave_text = new TPaveText(0.7, 0.8, 0.9, 0.9, "NDC");
    std::stringstream original_mean_text;
    original_mean_text << "Mean MPV: " << Form("%.1f", original_mean);
    original_pave_text->AddText(original_mean_text.str().c_str());
    original_pave_text->Draw();
    sector_canvas->SaveAs(sector_filename.str().c_str());

    std::stringstream adj_sector_filename;
    TCanvas* adj_sector_canvas = new TCanvas();
    adj_sector_filename << "./sector_diagrams/adjusted/sector" << sector << ".svg";
    adj_hs->SetMinimum(0.);
    adj_hs->SetMaximum(600.);
    adj_hs->Draw("nostack");
    adj_sector_canvas->SetGrid();
    double adjusted_mean = adjusted_sector_means[sector];
    TLine* adj_mean_line = new TLine(0.0, adjusted_mean, 95.0, adjusted_mean);
    //adj_mean_line->SetLineColor(kRed);
    adj_mean_line->SetLineStyle(2);
    adj_mean_line->SetLineWidth(2);
    adj_mean_line->Draw();
    TPaveText* adj_pave_text = new TPaveText(0.7, 0.8, 0.9, 0.9, "NDC");
    std::stringstream adj_mean_text;
    adj_mean_text << "Mean MPV: " << Form("%.1f", adjusted_mean);
    adj_pave_text->AddText(adj_mean_text.str().c_str());
    adj_pave_text->Draw();
    adj_sector_canvas->SaveAs(adj_sector_filename.str().c_str());

    std::stringstream combined_sector_filename;
    TCanvas* combined_sector_canvas = new TCanvas();
    combined_sector_filename << "./sector_diagrams/combined/sector" << sector << ".svg";
    combined_hs->SetMinimum(0.);
    combined_hs->SetMaximum(600.);
    combined_hs->Draw("nostack");
    combined_sector_canvas->SetGrid();
    combined_sector_canvas->SaveAs(combined_sector_filename.str().c_str());

    std::stringstream diff_sector_filename;
    TCanvas* diff_sector_canvas = new TCanvas();
    diff_sector_filename << "./sector_diagrams/diff/sector" << sector << ".svg";
    diff_hs->SetMinimum(-100.);
    diff_hs->SetMaximum(100.);
    diff_hs->Draw("nostack");
    diff_sector_canvas->SetGrid();
    diff_sector_canvas->SaveAs(diff_sector_filename.str().c_str());

    // create sector hists
    TCanvas* sector_hist_canvas = new TCanvas();
    //gStyle->SetOptFit(1);
    TFitResultPtr gaus_fit = sector_hist->Fit("gaus", "Q");
    TF1* gaus = (TF1*) sector_hist->GetListOfFunctions()->FindObject("gaus");
    double mean = gaus->GetParameter(1); original_means.push_back(mean);
    double mean_err = gaus->GetParError(1); original_mean_errs.push_back(mean_err);
    double sigma = gaus->GetParameter(2); original_sigmas.push_back(sigma);
    double sigma_err = gaus->GetParError(2); original_sigma_errs.push_back(sigma_err);
    sector_csv << sector << mean << mean_err << sigma << sigma_err << csvfile::endrow;
    sector_hist->Draw("E0 SAME");
    std::stringstream sector_hist_filename;
    sector_hist_filename << "./new_vop/sector/sector" << sector << ".svg";
    sector_hist_canvas->SaveAs(sector_hist_filename.str().c_str());

    TCanvas* adj_sector_hist_canvas = new TCanvas();
    //gStyle->SetOptFit(1);
    TFitResultPtr adj_gaus_fit = adj_sector_hist->Fit("gaus", "Q");
    TF1* adj_gaus = (TF1*) adj_sector_hist->GetListOfFunctions()->FindObject("gaus");
    double adj_mean = adj_gaus->GetParameter(1); adj_means.push_back(adj_mean);
    double adj_mean_err = adj_gaus->GetParError(1); adj_mean_errs.push_back(adj_mean_err);
    double adj_sigma = adj_gaus->GetParameter(2); adj_sigmas.push_back(adj_sigma);
    double adj_sigma_err = adj_gaus->GetParError(2); adj_sigma_errs.push_back(adj_sigma_err);
    adj_sector_csv << sector << adj_mean << adj_mean_err << adj_sigma << adj_sigma_err << csvfile::endrow;
    adj_sector_hist->Draw("E0 SAME");
    std::stringstream adj_sector_hist_filename;
    adj_sector_hist_filename << "./adjusted_vop/sector/sector" << sector << ".svg";
    adj_sector_hist_canvas->SaveAs(adj_sector_hist_filename.str().c_str());

    TCanvas* combined_sector_hist_canvas = new TCanvas();
    
    std::stringstream combined_hist_title;
    combined_hist_title << "Sector " << sector << " Before and After FB Adjustment;MPV;Num Blocks";
    THStack* combined_hist = new THStack("hs", combined_hist_title.str().c_str());
    sector_hist->SetLineColorAlpha(kRed, 0.5);
    //sector_hist->SetFillColorAlpha(kRed, 0.5);
    sector_hist->GetFunction("gaus")->SetLineColor(kRed);
    adj_sector_hist->SetLineColorAlpha(kBlue, 0.5);
    //adj_sector_hist->SetFillColorAlpha(kBlue, 0.5);
    adj_sector_hist->GetFunction("gaus")->SetLineColor(kBlue);
    combined_hist->Add(sector_hist);
    combined_hist->Add(adj_sector_hist);
    combined_hist->Draw("E0 NOSTACK");

    //the following lines will force the stats for h[1] and h[2]
    //to be drawn at a different position to avoid overlaps
    //gPad->Modified();
    //gPad->Update();
    combined_sector_hist_canvas->Update(); //to for the generation of the 'stat" boxes
    TPaveStats* st1 = (TPaveStats*) sector_hist->GetListOfFunctions()->FindObject("stats");
    TPaveStats* st2 = (TPaveStats*) adj_sector_hist->GetListOfFunctions()->FindObject("stats");
    st1->SetX1(0.8); st1->SetX2(1.0); st1->SetY1(0.45); st1->SetY2(0.75);
    st2->SetX1(0.8); st2->SetX2(1.0); st2->SetY1(0.15); st2->SetY2(0.45);
    //combined_sector_hist_canvas->Modified();
    
    TLegend* legend = new TLegend(0.8, 0.75, 1.0, 0.85);
    legend->AddEntry(sector_hist, "original");
    legend->AddEntry(adj_sector_hist, "adjusted");
    legend->Draw();

    std::stringstream combined_sector_hist_filename;
    combined_sector_hist_filename << "./new_vop/sector/sector" << sector << "_combined.svg";
    combined_sector_hist_canvas->SaveAs(combined_sector_hist_filename.str().c_str());
  }
  total_adjusted_mean = total_adjusted_mean / total_num_blocks;

  // do the old things
  csvfile csv("dbn_to_mpv.csv");
  csv << "dbn" << "mpv" << csvfile::endrow;
  for (auto const &p : dbn_mpv) {
    csv << p.first << p.second << csvfile::endrow;
  }
  // old plots for original mpv data
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
    TFitResultPtr gaus_fit = hist->Fit("gaus", "Q");
    Color_t color = new_vop_colors[vop];
    hist_colors.push_back(color);
    hist->SetLineColor(color);
    hist->GetFunction("gaus")->SetLineColor(color);
    hist->Draw("E0 SAME");
    new_hs->Add(hist);

    std::stringstream fileName;
    fileName << "./new_vop/vop_graphs/vop" << vop << ".svg";
    canvas->SaveAs(fileName.str().c_str());
    TF1* gaus = (TF1*) hist->GetFunction("gaus");
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
  new_graph2->SetTitle("Block MPV vs SiPM VOp");
  new_graph2->GetXaxis()->SetTitle("vop");
  new_graph2->GetYaxis()->SetTitle("mpv");
  new_graph2->GetXaxis()->SetLimits(68.85, 69.4);
  new_graph2->GetHistogram()->SetMinimum(220.);
  new_graph2->GetHistogram()->SetMaximum(380.);
  //new_graph2->GetYaxis()->SetLimits(220., 380.);
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

  // do the old things for adjusted mpv
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
  adjusted_graph2->GetHistogram()->SetMinimum(220.);
  adjusted_graph2->GetHistogram()->SetMaximum(380.);
  //adjusted_graph2->GetYaxis()->SetLimits(220., 380.);
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

  // calculate overall means and per-sector means
  double overall_orig_mean_mpv = total_original_mean;
  double overall_adj_mean_mpv = total_adjusted_mean;
  std::map<int, double> sector_orig_mean_mpv;
  std::map<int, double> sector_adj_mean_mpv;
  for (int sector : sectors) {
    auto tuples = tuple_data[sector];
    double sector_orig_avg = 0;
    double sector_adj_avg = 0;
    int num_blocks = 0;
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      double mpv = std::get<3>(block_info);
      double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) {
        sector_orig_avg += mpv;
        sector_adj_avg += adj_mpv;
        num_blocks++;
      }
    }
    sector_orig_avg = sector_orig_avg / num_blocks;
    sector_adj_avg = sector_adj_avg / num_blocks;
    sector_orig_mean_mpv[sector] = sector_orig_avg;
    sector_adj_mean_mpv[sector] = sector_adj_avg;
  }

  std::vector<double> sectors_as_double;
  std::vector<double> sector_err;
  for (int i = 0; i < sectors.size(); i++) {
    sectors_as_double.push_back(sectors[i]);
    sector_err.push_back(0.0);
  }

  // create sigma, mean, and sigma/mean plots
  auto orig_get_err_div_result = get_err_div(original_sigmas, original_sigma_errs, original_means, original_mean_errs);
  std::vector<double> orig_sigma_over_means = orig_get_err_div_result.first;
  std::vector<double> orig_sigma_over_mean_errs = orig_get_err_div_result.second;
  auto adj_get_err_div_result = get_err_div(adj_sigmas, adj_sigma_errs, adj_means, adj_mean_errs);
  std::vector<double> adj_sigma_over_means = adj_get_err_div_result.first;
  std::vector<double> adj_sigma_over_mean_errs = adj_get_err_div_result.second;
  
  // original means
  TCanvas* original_mean_canvas = new TCanvas();
  TGraphErrors* original_mean_graph = new TGraphErrors(sectors.size(), &sectors_as_double[0], &original_means[0], &sector_err[0], &original_mean_errs[0]);
  original_mean_graph->SetTitle("Sector MPV Means (Original)");
  original_mean_graph->GetXaxis()->SetTitle("Sector");
  original_mean_graph->GetYaxis()->SetTitle("Mean (MPV)");
  original_mean_graph->GetXaxis()->SetLimits(0.0, max_sector + 1);
  original_mean_graph->GetHistogram()->SetMinimum(0.0);
  original_mean_graph->GetHistogram()->SetMaximum(600.0);
  //original_mean_graph->GetYaxis()->SetLimits(220., 380.);
  original_mean_graph->SetMarkerStyle(33); // diamond
  original_mean_graph->Draw("AP");
  original_mean_canvas->SetGrid();
  original_mean_canvas->SaveAs("./sector_diagrams/original_means.svg");
  
  // adjusted means
  TCanvas* adj_mean_canvas = new TCanvas();
  TGraphErrors* adj_mean_graph = new TGraphErrors(sectors.size(), &sectors_as_double[0], &adj_means[0], &sector_err[0], &adj_mean_errs[0]);
  adj_mean_graph->SetTitle("Sector MPV Means (Adjusted)");
  adj_mean_graph->GetXaxis()->SetTitle("Sector");
  adj_mean_graph->GetYaxis()->SetTitle("Mean (MPV)");
  adj_mean_graph->GetXaxis()->SetLimits(0.0, max_sector + 1);
  adj_mean_graph->GetHistogram()->SetMinimum(0.0);
  adj_mean_graph->GetHistogram()->SetMaximum(600.0);
  //adj_mean_graph->GetYaxis()->SetLimits(220., 380.);
  adj_mean_graph->SetMarkerStyle(33); // diamond
  adj_mean_graph->Draw("AP");
  adj_mean_canvas->SetGrid();
  adj_mean_canvas->SaveAs("./sector_diagrams/adj_means.svg");

  // original sigmas
  TCanvas* original_sigma_canvas = new TCanvas();
  TGraphErrors* original_sigma_graph = new TGraphErrors(sectors.size(), &sectors_as_double[0], &original_sigmas[0], &sector_err[0], &original_sigma_errs[0]);
  original_sigma_graph->SetTitle("Sector MPV Sigmas (Original)");
  original_sigma_graph->GetXaxis()->SetTitle("Sector");
  original_sigma_graph->GetYaxis()->SetTitle("Sigma (MPV)");
  original_sigma_graph->GetXaxis()->SetLimits(0.0, max_sector + 1);
  original_sigma_graph->GetHistogram()->SetMinimum(0.0);
  original_sigma_graph->GetHistogram()->SetMaximum(100.0);
  //original_sigma_graph->GetYaxis()->SetLimits(220., 380.);
  original_sigma_graph->SetMarkerStyle(33); // diamond
  original_sigma_graph->Draw("AP");
  original_sigma_canvas->SetGrid();
  original_sigma_canvas->SaveAs("./sector_diagrams/original_sigmas.svg");

  // adjsuted sigmas
  TCanvas* adj_sigma_canvas = new TCanvas();
  TGraphErrors* adj_sigma_graph = new TGraphErrors(sectors.size(), &sectors_as_double[0], &adj_sigmas[0], &sector_err[0], &adj_sigma_errs[0]);
  adj_sigma_graph->SetTitle("Sector MPV Sigmas (Adjusted)");
  adj_sigma_graph->GetXaxis()->SetTitle("Sector");
  adj_sigma_graph->GetYaxis()->SetTitle("Sigma (MPV)");
  adj_sigma_graph->GetXaxis()->SetLimits(0.0, max_sector + 1);
  adj_sigma_graph->GetHistogram()->SetMinimum(0.0);
  adj_sigma_graph->GetHistogram()->SetMaximum(100.0);
  //adj_sigma_graph->GetYaxis()->SetLimits(220., 380.);
  adj_sigma_graph->SetMarkerStyle(33); // diamond
  adj_sigma_graph->Draw("AP");
  adj_sigma_canvas->SetGrid();
  adj_sigma_canvas->SaveAs("./sector_diagrams/adj_sigmas.svg");

  // original sigma_over_means
  TCanvas* original_sigma_over_mean_canvas = new TCanvas();
  TGraphErrors* original_sigma_over_mean_graph = new TGraphErrors(sectors.size(), &sectors_as_double[0], &orig_sigma_over_means[0], &sector_err[0], &orig_sigma_over_mean_errs[0]);
  original_sigma_over_mean_graph->SetTitle("Sector MPV Sigma/Mean (Original)");
  original_sigma_over_mean_graph->GetXaxis()->SetTitle("Sector");
  original_sigma_over_mean_graph->GetYaxis()->SetTitle("Sigma/Mean");
  original_sigma_over_mean_graph->GetXaxis()->SetLimits(0.0, max_sector + 1);
  original_sigma_over_mean_graph->GetHistogram()->SetMinimum(0.0);
  original_sigma_over_mean_graph->GetHistogram()->SetMaximum(1.0);
  //original_sigma_over_mean_graph->GetYaxis()->SetLimits(220., 380.);
  original_sigma_over_mean_graph->SetMarkerStyle(33); // diamond
  original_sigma_over_mean_graph->Draw("AP");
  original_sigma_over_mean_canvas->SetGrid();
  original_sigma_over_mean_canvas->SaveAs("./sector_diagrams/original_sigma_over_means.svg");

  // adjsuted sigma_over_means
  TCanvas* adj_sigma_over_mean_canvas = new TCanvas();
  TGraphErrors* adj_sigma_over_mean_graph = new TGraphErrors(sectors.size(), &sectors_as_double[0], &adj_sigma_over_means[0], &sector_err[0], &adj_sigma_over_mean_errs[0]);
  adj_sigma_over_mean_graph->SetTitle("Sector MPV Sigma/Mean (Adjusted)");
  adj_sigma_over_mean_graph->GetXaxis()->SetTitle("Sector");
  adj_sigma_over_mean_graph->GetYaxis()->SetTitle("Sigma/Mean");
  adj_sigma_over_mean_graph->GetXaxis()->SetLimits(0.0, max_sector + 1);
  adj_sigma_over_mean_graph->GetHistogram()->SetMinimum(0.0);
  adj_sigma_over_mean_graph->GetHistogram()->SetMaximum(1.0);
  //adj_sigma_over_mean_graph->GetYaxis()->SetLimits(220., 380.);
  adj_sigma_over_mean_graph->SetMarkerStyle(33); // diamond
  adj_sigma_over_mean_graph->Draw("AP");
  adj_sigma_over_mean_canvas->SetGrid();
  adj_sigma_over_mean_canvas->SaveAs("./sector_diagrams/adj_sigma_over_means.svg");

  // PROOF THAT MY CODE IS CORRECT @ANNE

  std::map<int, std::map<int, int>> sector_fb_counts;
  for (int sector : sectors) {
    auto tuples = tuple_data[sector];
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      double mpv = std::get<3>(block_info);
      double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) {
        sector_fb_counts[sector][fiber_batch]++;
      }
    }
  }
  
  // make fb histograms for each sector
  std::map<int, std::map<int, Color_t>> sector_fb_colors;
  std::map<int, std::vector<int>> sector_batches;
  for (int sector : sectors) {
    auto tuples = tuple_data[sector];
    // first determine the number of batches in this sector
    std::vector<int> batches;
    std::map<int, Color_t> fb_colors;
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      //int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      double mpv = std::get<3>(block_info);
      //double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) {
        if (std::find(batches.begin(), batches.end(), fiber_batch) == batches.end()) {
          // add this fb
          batches.push_back(fiber_batch);
          // create new color for this fb
          Color_t ci = TColor::GetFreeColorIndex();
          int nth_kelly = fb_colors.size();
          TColor* color = new TColor(ci, kelly_colors[nth_kelly][0], kelly_colors[nth_kelly][1], kelly_colors[nth_kelly][2]);
          fb_colors[fiber_batch] = ci;
          //std::cout << "sector " << sector << " fb " << fiber_batch << " has kelly color " << nth_kelly << std::endl;
        }
      }
    }
    std::cout << "SECTOR " << sector << " HAS " << batches.size() << " BATCHES" << std::endl;
    std::sort(batches.begin(), batches.end());
    sector_batches[sector] = batches;
    sector_fb_colors[sector] = fb_colors;
    // generate bin labels
    std::vector<std::string> bin_labels;
    for (int batch : batches) {
      double correction_factor = total_original_mean / fiber_batch_avg_mpv[batch];
      std::stringstream label;
      label << "FB " << batch << ": " << Form("%.3f", correction_factor);
      bin_labels.push_back(label.str());
    }
    std::stringstream hs_title;
    hs_title << "Sector " << sector << ";Fiber Batch;Num Blocks";
    THStack* hs = new THStack("hs", hs_title.str().c_str());
    std::map<int, TH1D*> batch_hists;
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      //int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      double mpv = std::get<3>(block_info);
      //double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) {
        if (batch_hists.find(fiber_batch) == batch_hists.end()) {
          // create hist for this fb
          std::stringstream batch_hist_title;
          batch_hist_title << "fb " << fiber_batch << ";Fiber Batch;Num Blocks";
          int num_batches = batches.size();
          TH1D* batch_hist = new TH1D("h", batch_hist_title.str().c_str(), num_batches, 0, num_batches);
          batch_hist->GetSumw2();
          batch_hist->SetLineColor(fb_colors[fiber_batch]);
          batch_hist->SetFillColor(fb_colors[fiber_batch]);
          //batch_hist->GetXaxis()->LabelsOption("v");
          for (int i = 0; i < bin_labels.size(); i++) {
            batch_hist->GetXaxis()->SetBinLabel(i + 1, std::to_string(batches[i]).c_str());
            //batch_hist->GetXaxis()->SetBinLabel(i + 1, bin_labels[i].c_str());
            //batch_hist->GetXaxis()->SetBinLabel(i + 1, "");
          }
          batch_hists[fiber_batch] = batch_hist;
        }
        batch_hists[fiber_batch]->Fill(std::find(batches.begin(), batches.end(), fiber_batch) - batches.begin());
      }
    }
    for (auto const &p : batch_hists) {
      hs->Add(p.second, "B");
    }
    TCanvas* sector_canvas = new TCanvas();
    sector_canvas->SetGrid();
    sector_canvas->SetRightMargin(0.3);
    hs->Draw();
    TLegend* legend = new TLegend(0.725, 0.025, 0.975, 0.975);
    for (int i = 0; i < batches.size(); i++) {
      legend->AddEntry(batch_hists[batches[i]], bin_labels[i].c_str());
    }
    legend->Draw();
    std::stringstream hs_filename;
    hs_filename << "./proof/batches/sector" << sector << ".svg";
    sector_canvas->SaveAs(hs_filename.str().c_str());
  }

  std::cout << "FIBER BATCH COUNTS:" << std::endl;
  for (auto const &p : sector_fb_counts) {
    int sector = p.first;
    std::cout << "sector " << sector << std:: endl;
    for (auto const &q : p.second) {
      int batch = q.first;
      int count = q.second;
      std::cout << "  fb " << batch << ": " << count << std::endl;
    }
  }

  // make calibration factor plots for each sector
  for (int sector : sectors) {
    std::cout << "SECTOR " << sector << ":" << std::endl;
    auto tuples = tuple_data[sector];
    std::stringstream hs_title;
    hs_title << "Sector " << sector << ";Block Number;Correction Factor";
    THStack* hs = new THStack("hs", hs_title.str().c_str());
    std::map<int, std::pair<TH1D*, TH1D*>> batch_hists;
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      //int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      double mpv = std::get<3>(block_info);
      //double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) {
        if (batch_hists.find(fiber_batch) == batch_hists.end()) {
          // create hist for this fb
          std::stringstream batch_hist_up_title;
          batch_hist_up_title << "fb " << fiber_batch << " up;Block Number;Correction Factor";
          TH1D* batch_hist_up = new TH1D("h_up", batch_hist_up_title.str().c_str(), 96, 0, 96);
          batch_hist_up->GetSumw2();
          batch_hist_up->SetMarkerColor(sector_fb_colors[sector][fiber_batch]);
          batch_hist_up->SetMarkerStyle(22); // up triangle
          // and the down hist
          std::stringstream batch_hist_down_title;
          batch_hist_down_title << "fb " << fiber_batch << " down;Block Number;Correction Factor";
          TH1D* batch_hist_down = (TH1D*) batch_hist_up->Clone("h_down");
          batch_hist_down->SetTitle(batch_hist_down_title.str().c_str());
          batch_hist_down->SetMarkerStyle(23); // down triangle
          // add to the map
          batch_hists[fiber_batch] = std::make_pair(batch_hist_up, batch_hist_down);
          //std::cout << "CREATED HISTS FOR FB " << fiber_batch << std::endl;
        }
        double correction_factor = total_original_mean / fiber_batch_avg_mpv[fiber_batch];
        if (correction_factor >= 1) {
          batch_hists[fiber_batch].first->SetBinContent(block_num + 1, correction_factor);
          //batch_hists[fiber_batch].second->SetBinContent(block_num + 1, -9000);
          //std::cout << "added an up arrow at block " << block_num << " fb " << fiber_batch << " factor " << correction_factor << std::endl;
        } else {
          //batch_hists[fiber_batch].first->SetBinContent(block_num + 1, -9000);
          batch_hists[fiber_batch].second->SetBinContent(block_num + 1, correction_factor);
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
    //sector_canvas->SetRightMargin(0.3);
    hs->SetMinimum(0.75);
    hs->SetMaximum(1.25);
    hs->Draw("NOSTACK");
    /*
    TLegend* legend = new TLegend(0.725, 0.025, 0.975, 0.975);
    for (int i = 0; i < batches.size(); i++) {
      legend->AddEntry(batch_hists[batches[i]], bin_labels[i].c_str());
    }
    legend->Draw();
    */
    TLine* one_line = new TLine(0.0, 1.0, 96.0, 1.0);
    //one_line->SetLineColor(kRed);
    one_line->SetLineStyle(2);
    one_line->SetLineWidth(2);
    one_line->Draw();
    std::stringstream hs_filename;
    hs_filename << "./proof/tim/sector" << sector << ".svg";
    sector_canvas->SaveAs(hs_filename.str().c_str());
  }

  // make combined tim graphs, colored by fiber batch instead of interface board
  for (int sector : sectors) {
    std::cout << "SECTOR " << sector << ":" << std::endl;
    auto tuples = tuple_data[sector];
    std::stringstream hs_title;
    hs_title << "Sector " << sector << ";Block Number;MPV";
    THStack* hs = new THStack("hs", hs_title.str().c_str());
    std::map<int, std::pair<std::pair<TH1D*, TH1D*>, TH1D*>> batch_hists;
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      //int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      double mpv = std::get<3>(block_info);
      double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) {
        if (batch_hists.find(fiber_batch) == batch_hists.end()) {
          // original hist for adj up
          std::stringstream batch_hist_orig_up_title;
          batch_hist_orig_up_title << "fb " << fiber_batch << " orig up;Block Number;MPV";
          TH1D* batch_hist_orig_up = new TH1D("h_up", batch_hist_orig_up_title.str().c_str(), 96, 0, 96);
          batch_hist_orig_up->GetSumw2();
          batch_hist_orig_up->SetMarkerStyle(22); // solid up arrow
          batch_hist_orig_up->SetMarkerColorAlpha(sector_fb_colors[sector][fiber_batch], 0.5);
          // original hist for adj down
          std::stringstream batch_hist_orig_down_title;
          batch_hist_orig_down_title << "fb " << fiber_batch << " orig down;Block Number;MPV";
          TH1D* batch_hist_orig_down = (TH1D*) batch_hist_orig_up->Clone("h_down");
          batch_hist_orig_down->SetTitle(batch_hist_orig_down_title.str().c_str());
          batch_hist_orig_down->SetMarkerStyle(23); // solid down arrow
          batch_hist_orig_down->SetMarkerColorAlpha(sector_fb_colors[sector][fiber_batch], 0.5);
          // adj hist
          std::stringstream batch_hist_adj_title;
          batch_hist_adj_title << "fb " << fiber_batch << " adj;Block Number;MPV";
          TH1D* batch_hist_adj = (TH1D*) batch_hist_orig_up->Clone("h_adj");
          batch_hist_adj->SetTitle(batch_hist_adj_title.str().c_str());
          batch_hist_adj->SetMarkerStyle(20); // solid circle
          //batch_hist_adj->SetMarkerColor(sector_fb_colors[sector][fiber_batch]);
          batch_hist_adj->SetMarkerColor(sector_fb_colors[sector][fiber_batch]);
          // add to the map
          batch_hists[fiber_batch] = std::make_pair(std::make_pair(batch_hist_orig_up, batch_hist_orig_down), batch_hist_adj);
          //std::cout << "CREATED HISTS FOR FB " << fiber_batch << std::endl;
        }
        std::cout << "ADJ MPV IS " << adj_mpv << std::endl;
        batch_hists[fiber_batch].second->SetBinContent(block_num + 1, adj_mpv);
        double correction_factor = total_original_mean / fiber_batch_avg_mpv[fiber_batch];
        if (correction_factor >= 1) {
          batch_hists[fiber_batch].first.first->SetBinContent(block_num + 1, mpv);
        } else {
          batch_hists[fiber_batch].first.second->SetBinContent(block_num + 1, mpv);
        }
      }
    }
    for (auto const &p : batch_hists) {
      hs->Add(p.second.first.first, "P");
      hs->Add(p.second.first.second, "P");
      hs->Add(p.second.second, "P");
    }
    TCanvas* sector_canvas = new TCanvas();
    sector_canvas->SetGrid();
    //sector_canvas->SetRightMargin(0.3);
    hs->SetMinimum(0.0);
    hs->SetMaximum(600.0);
    hs->Draw("NOSTACK");
    /*
    TLegend* legend = new TLegend(0.725, 0.025, 0.975, 0.975);
    for (int i = 0; i < batches.size(); i++) {
      legend->AddEntry(batch_hists[batches[i]], bin_labels[i].c_str());
    }
    legend->Draw();
    */
    std::stringstream hs_filename;
    hs_filename << "./proof/combined/sector" << sector << ".svg";
    sector_canvas->SaveAs(hs_filename.str().c_str());
  }



  for (int i = 0; i < sectors.size(); i++) {
    std::cout << "sector " << sectors[i] << " mean " << original_means[i] << " -> " << adj_means[i] << ", sigma " << original_sigmas[i] << " -> " << adj_sigmas[i] << std::endl;
  }
  std::cout << "max diff: " << max_diff << std::endl;
  std::cout << "max ratio diff: " << max_ratio_diff << std::endl;
  std::cout << "mean over all sectors " << total_original_mean << " -> " << total_adjusted_mean << std::endl;

  for (int sector : sectors) {
    int count = 0;
    auto tuples = tuple_data[sector];
    for (int block_num = 0; block_num < 96; block_num++) {
      auto block_info = tuples[block_num];
      bool has_mpv_data = std::get<0>(block_info);
      //int dbn = std::get<1>(block_info);
      int fiber_batch = std::get<2>(block_info);
      //double mpv = std::get<3>(block_info);
      //double adj_mpv = std::get<4>(block_info);
      if (has_mpv_data && fiber_batch >= 0) { count++; }
    }
    if (count == 96) {
      //std::cout << "sector " << sector << " looks good :)" << std::endl;
    } else {
      //std::cout << "sector " << sector << " has " << count << " != 96 blocks with proper data >:(" << std::endl;
    }
  }

  // print sector means and their deviations from their respective means
  double orig_avg_sector_mean = 0;
  double adj_avg_sector_mean = 0;
  for (int sector : sectors) {
    orig_avg_sector_mean += sector_orig_mean_mpv[sector];
    adj_avg_sector_mean += sector_adj_mean_mpv[sector];
    double orig_dev_mean = sector_orig_mean_mpv[sector] - overall_orig_mean_mpv;
    double adj_dev_mean = sector_adj_mean_mpv[sector] - overall_adj_mean_mpv;
    std::cout << "sector " << sector << ": " << sector_orig_mean_mpv[sector] << " (" << orig_dev_mean << ") -> " << sector_adj_mean_mpv[sector] << " (" << adj_dev_mean << ")" << std::endl;
  }
  orig_avg_sector_mean = orig_avg_sector_mean / sectors.size();
  adj_avg_sector_mean = adj_avg_sector_mean / sectors.size();
  std::cout << "avg of sector means: " << orig_avg_sector_mean << " (" << orig_avg_sector_mean - overall_orig_mean_mpv << ") -> " << adj_avg_sector_mean << " (" << adj_avg_sector_mean - overall_adj_mean_mpv << ")" << std::endl;
}
