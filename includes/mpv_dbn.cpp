#include "mpv_dbn.h"
#include <TMath.h>

constexpr bool debug = false;

constexpr char PRINT_BOLD[] = "\x1b[1m";
constexpr char PRINT_GREY[] = "\x1b[0;90m";
constexpr char PRINT_RED[] = "\x1b[1;31m";
constexpr char PRINT_YELLOW[] = "\x1b[1;33m";
constexpr char PRINT_GREEN[] = "\x1b[1;32m";
constexpr char PRINT_BLUE[] = "\x1b[1;34m";
constexpr char PRINT_END[] = "\x1b[0m";

/**
 * @brief The lower nonphysical threshold for MPV. MPVs <= this value will be saved as -1.
 */
constexpr double MPV_CUTOFF_LOW = 0.0;
/**
 * @brief The upper nonphysical threshold for MPV. MPVs >= this value will be saved as -1.
 */
constexpr double MPV_CUTOFF_HIGH = 1000.0;

/**
 * @brief Get the channel numbers associated with a block.
 * 
 * @param block_num block number (/96, 0-based, corresponding to h_allblocks).
 * @return std::vector<int> channel numbers (0-based, corresponding to h_allchannels).
 */
std::vector<int> block_to_channel(int block_num) {
  int t_x_off = block_num%4;
  int t_y_off = block_num/4;
  int f_x_off = t_y_off;
  int f_y_off = 3 - t_x_off;
  int ib = f_x_off/4;
  int ib_block_idx = 3 - f_x_off%4 + 4*f_y_off + 16*ib;
  int chnl = ib_block_idx*4;
  return {chnl, chnl + 1, chnl + 2, chnl + 3};
}

/**
 * @brief Get the block number of a channel.
 * 
 * @param channel channel number (0-based, corresponding to h_allchannels).
 * @return int block number (/96, 0-based, corresponding to h_allblocks).
 */
int channel_to_block(int channel) {
  int i = channel/4;
  int ib = i/16;
  int f_x_off = 3 - ((i%16)%4) + 4*ib;
  int f_y_off = (i%16)/4;
  int t_x_off = 3 - f_y_off;
  int t_y_off = f_x_off;
  return t_x_off + 4*t_y_off;
}

/**
 * @brief Get the channel numbers on the perimeter of a sector.
 * 
 * @param drop_low_rap_edge whether to count low rapidity edge as part of the perimeter.
 * @return std::set<int> set of channel numbers (0-based) which are on the perimeter of a sector. 
 */
std::set<int> perimeter_channels(bool drop_low_rap_edge = true) {
  std::set<int> perimeter;
  std::vector<std::vector<int>> channels(48);
  for (auto &vec : channels) {
    vec = std::vector<int>(8, -1);
  }

  for (int ib = 0; ib < 6; ib++) {
    int x_ib_off = 4*2*ib;
    for (int i_block = 0; i_block < 16; i_block++) {
      int x_block_off = 2*(3 - (i_block%4));
      int y_block_off = 2*(i_block/4);
      for (int channel = 0; channel < 4; channel++) {
        int x_channel_off = channel%2;
        int y_channel_off = channel/2;
        int x_off = x_ib_off + x_block_off + x_channel_off;
        int y_off = y_block_off + y_channel_off;
        int channel_num = 64*ib + 4*i_block + channel;
        // printf("(%2d, %2d)\n", x_ib_off + x_block_off + x_channel_off, y_block_off + y_channel_off);
        channels[x_off][y_off] = channel_num;
        if (!drop_low_rap_edge && (x_off == 47 || y_off == 0 || y_off == 7)) {          
          perimeter.insert(channel_num);
        } else if (drop_low_rap_edge && (x_off == 0 || x_off == 47 || y_off == 0 || y_off == 7)) {
          perimeter.insert(channel_num);
        }
      }
    }
  }
  // for (auto &vec : channels) {
  //   for (auto &chnl : vec) {
  //     printf("%3d ", chnl);
  //   }
  //   printf("\n");
  // }
  return perimeter;
}

/**
 * @brief Get the channel numbers which contribute to a block's block mpv.
 * 
 * @param perimeter set of channels to be excluded from block mpv calculation (edges).
 * @param block block number (/96, 0-based).
 * @return std::vector<int> channel numbers (0-based) which contribute to block's mpv calculation.
 */
std::vector<int> block_contributing_channels(const std::set<int> &perimeter, int block) {
  std::vector<int> all_chnls = block_to_channel(block);
  std::vector<int> contributing_channels;
  for (int &chnl : all_chnls) {
    if (perimeter.find(chnl) == perimeter.end()) {
      contributing_channels.push_back(chnl);
    }
  }
  return contributing_channels;
}

/**
 * @brief Get run numbers for each sector from csv. Run numbers are 1-based.

 * NOTE: depends on files/physics_runs.csv.
 * 
 * @return std::map<int, int> (sector -> run number).
 */
std::map<int, int> read_physics_runs() {
  std::map<int, int> new_sector_runs;
  std::fstream runs_file;
  runs_file.open("files/physics_runs.csv", std::ios::in);
  std::pair<int, int> *row;
  std::string line, word;
  int line_num = 1;
  while (std::getline(runs_file, line) && line_num <= 65) {
    if (line_num == 1) {
      line_num++;
      continue;
    }
    std::stringstream s(line);
    int sector;
    int run;
    if (sscanf(line.c_str(), "%i, %i", &sector, &run) != 2) {
      throw std::runtime_error("failed to parse 'int, int' in physics_runs.csv");
    }
    // else we have read ints sector, run
    new_sector_runs[sector] = run;
    line_num++;
  }
  return new_sector_runs;
}

/**
 * @brief Copy Tim's run folders that don't exist locally.
 * 
 * NOTE: must be run on SDCC.
 */
void get_physics_runs() {
  std::map<int, int> sector_runs = read_physics_runs();
  for (int sector = 1; sector <= 64; sector++) {
    int run = sector_runs[sector];
    if (run > 0) {
      // check if local folder exists
      char *local_folder_name = Form("physics_runs/qa_output_000%05d", run);
      if (opendir(local_folder_name)) {
        // folder exists locally!
        continue;
      } else if (errno == ENOENT) {
        // folder does not exist locally
        char *server_folder_name = Form("/gpfs/mnt/gpfs02/sphenix/user/trinn/sPHENIX_emcal_cosmics_sector0/macros/qa_output_000%05d/", run);
        // check if server folder exists!
        if (opendir(server_folder_name)) {
          // it worked!
        } else if (errno == ENOENT) {
          throw std::runtime_error(Form("%s does not exist!", server_folder_name));
        } else {
          throw std::runtime_error(Form("opendir failed"));
        }
        system(Form("cp -r %s physics_runs", server_folder_name));
        // check if it worked...
        if (opendir(local_folder_name)) {
          // it worked!
          printf("got new physics fun for sector %2d\n", sector);
          continue;
        } else if (errno == ENOENT) {
          throw std::runtime_error(Form("failed to cp %s -> physics_runs", server_folder_name));
        } else {
          throw std::runtime_error(Form("opendir failed"));
        }
      } else {
        throw std::runtime_error(Form("opendir failed"));
      }
    } // else we don't have the info OR it's in sector_data already
  }
}

/**
 * @brief Get DBNs for each block for each sector.
 * 
 * NOTE: depends on files/Blocks database - Sectors.csv.
 * 
 * @return std::vector<std::vector<std::string>> [sector][block number] -> dbn, or "" if none.
 */
std::vector<std::vector<std::string>> get_dbns() {
  std::fstream sector_map_file;
  sector_map_file.open("files/Blocks database - Sectors.csv", std::ios::in);
  std::vector<std::vector<std::string>> all_data(24);
  std::vector<std::string> *row;
  std::string line, word;
  int line_num = 1;
  while (std::getline(sector_map_file, line) && line_num < 26) {
    if (line_num == 1) {
      line_num++;
      continue;
    }
    row = &all_data[line_num - 2];
    std::stringstream s(line);
    size_t n = 0;
    while (std::getline(s, word, ',')) {
      if (word[word.length() - 1] == '\r') {
        word.erase(word.end() - 1);
      }
      row->push_back(word);
      n++;
    }
    if (n != 64*8) {
      throw std::runtime_error("expected row to have 64*8 cells");
    }
    line_num++;
  }

  std::vector<std::vector<std::string>> dbns(64);
  for (short i = 0; i < 64; i++) {
    size_t data_start = 4 + 8 * i;
    for (short j = 0; j < 24; j++) {
      for (short k = 0; k < 4; k++) {
        std::string dbn = all_data[j][data_start + k];
        if (dbn == "#N/A") {
          dbns[i].push_back("");
        } else {
          int d = 0;
          if (dbn[0] == 'F') {
            // expect that the rest is castable to int
            int res = sscanf(dbn.c_str(), "F%i", &d);
            if (res != 1) {
              throw std::runtime_error("expected dbn to be F[int]");
            }
          } else if (dbn[0] == 'C') {
            // expect that the rest is castable to int
            int res = sscanf(dbn.c_str(), "C%i", &d);
            if (res != 1) {
              throw std::runtime_error("expected dbn to be C[int]");
            }
          } else {
            // expect dbn is castable to int
            int res = sscanf(dbn.c_str(), "%i", &d);
            if (res != 1) {
              throw std::runtime_error("expected dbn to be [int]");
            }
          }
          dbns[i].push_back(dbn);
        } 
      }
    }
  }
  for (short sector = 0; sector < 64; sector++) {
    int n_blocks = 0;
    for (short block = 0; block < 96; block++) {
      std::string dbn = dbns[sector][block];
      if (dbn != "") {
        n_blocks++;
      }
    }
    if (debug) {
      if (n_blocks == 96) {
        printf("sector %2d: got dbn for all blocks!\n", sector + 1);
      } else if (n_blocks > 0) {
        printf("sector %2d: got dbn for %2d/96 blocks\n", sector + 1, n_blocks);
      }
    }
  }
  return dbns;
}

/**
 * @brief Get MPVs for each channel for each sector from h_allchannels.
 * 
 * @return std::vector<std::vector<double>> [sector][channel number] -> mpv 
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> get_chnl_mpv_with_err() {
  std::vector<std::vector<double>> chnl_mpv(64);
  std::vector<std::vector<double>> chnl_mpv_err(64);
  for (auto &vec : chnl_mpv) {
    vec = std::vector<double>(384, -1.0);
  }
  for (auto &vec : chnl_mpv_err) {
    vec = std::vector<double>(384, -1.0);
  }

  TFile *hist_file;
  char filename[64];
  std::map<int, int> physics_runs = read_physics_runs();
  for (std::pair<int, int> p : physics_runs) {
    int sector = p.first;
    int run_num = p.second;
    if (run_num > 0) {
      // get data from the run number
      sprintf(filename, "physics_runs/qa_output_000%i/histograms.root ", run_num);
      if ((hist_file = TFile::Open(filename))) {
        if (debug) {
          printf("found run file for sector %i: physics_runs/qa_output_000%i/histograms.root\n", sector, run_num);
        }
        TH1D* data;
        hist_file->GetObject("h_allchannels;1", data);
        if (!data) {
          std::cerr << "  unable to get histogram" << std::endl;
          continue;       
        }
        for (int chnl = 0; chnl < 384; chnl++) {
          double this_err = data->GetBinError(chnl + 1);
          double this_mpv = data->GetBinContent(chnl + 1);
          chnl_mpv[sector - 1][chnl] = this_mpv;
          chnl_mpv_err[sector - 1][chnl] = this_err;
        }
      } else {
        printf("FAILED to find run file for sector %i: qa_output_000%i/histograms.root\n", sector, run_num);
      }
    }
  }
  return std::make_pair(chnl_mpv, chnl_mpv_err);
}

/**
 * @brief Calculate block MPVs and error for each block for each sector from channel mpv and error.
 * 
 * @param chnl_mpv_with_err channel mpv and error.
 * @param perimeter set of channels to be excluded from block mpv calculation (edges).
 * @return std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> calculate_block_mpv_with_err(const std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> &chnl_mpv_with_err, const std::set<int> &perimeter) {
  std::vector<std::vector<double>> block_mpvs(64);
  std::vector<std::vector<double>> block_mpv_errs(64);
  for (auto &vec : block_mpvs) {
    vec = std::vector<double>(96, -1.0);
  }
  for (auto &vec : block_mpv_errs) {
    vec = std::vector<double>(96, -1.0);
  }
  for (int sector = 0; sector < 64; sector++) {
    for (int block = 0; block < 96; block++) {
      std::vector<int> chnls = block_contributing_channels(perimeter, block);
      double block_mpv = 0.0;
      double block_mpv_err = 0.0;
      for (int &chnl : chnls) {
        block_mpv += chnl_mpv_with_err.first[sector][chnl];
        block_mpv_err += chnl_mpv_with_err.second[sector][chnl];
      }
      block_mpv /= chnls.size();
      block_mpv_err /= chnls.size();
      block_mpvs[sector][block] = block_mpv;
      block_mpv_errs[sector][block] = block_mpv_err;
    }
  }
  return std::make_pair(block_mpvs, block_mpv_errs);
}

/**
 * @brief Get single pixel gaps for each block for each sector (average of block's four towers).
 * 
 * @param write_ib write IB mean and sigma to csv file.
 * @return std::vector<std::vector<double>> [sector][block number] -> gap.
 */
std::vector<std::vector<double>> get_sp_gaps(bool write_ib = false) {
  std::vector<std::vector<double>> sp_gaps(64);
  for (auto &vec : sp_gaps) {
    vec = std::vector<double>(96, -1.0);
  }

  TFile *hist_file;
  char filename[64];
  std::map<int, int> physics_runs = read_physics_runs();
  for (std::pair<int, int> p : physics_runs) {
    int sector = p.first;
    int run_num = p.second;
    if (run_num > 0) {
      // get data from the run number
      sprintf(filename, "physics_runs/qa_output_000%i/histograms.root", run_num);
      if ((hist_file = TFile::Open(filename))) {
        if (debug) {
          printf("found run file for sector %i: physics_runs/qa_output_000%i/histograms.root\n", sector, run_num);
        }
        TH1D* data;
        hist_file->GetObject("h_sp_perchnl;1", data);
        if (!data) {
          std::cerr << "  unable to get sp histogram" << std::endl;
          continue;       
        }
        for (int block_num = 0; block_num < 96; block_num++) {
          double avg_sp_gap = 0;
          for (int &chnl : block_to_channel(block_num)) {
            double content = data->GetBinContent(chnl + 1);
            if (content <= 0) {
              printf("complaint at sector %i channel %i: sp_gap <= 0 (%f)\n", sector, chnl, content);
            }
            avg_sp_gap += content;
          }
          sp_gaps[sector - 1][block_num] = avg_sp_gap/4;
          // printf("sector %2d block %2d: sp gap = %f\n", sector, block_num + 1, sp_gaps[sector - 1][block_num]);
        }
      } else {
        printf("FAILED to find run file for sector %i: qa_output_000%i/histograms.root\n", sector, run_num);
      }
    }
  }
  if (write_ib) {
    FILE *mean_outfile = fopen("files/sp_gap_avg_ib_mean.csv", "w+");
    FILE *sigma_outfile = fopen("files/sp_gap_avg_ib_sigma.csv", "w+");
    fprintf(mean_outfile, "SECTOR, RUN, IB0, IB1, IB2, IB3, IB4, IB5");
    fprintf(sigma_outfile, "SECTOR, RUN, IB0, IB1, IB2, IB3, IB4, IB5");
    for (short sector = 0; sector < 64; sector++) {
      fprintf(mean_outfile, "\n%i, %i", sector + 1, physics_runs[sector + 1]);
      fprintf(sigma_outfile, "\n%i, %i", sector + 1, physics_runs[sector + 1]);
      for (short ib = 0; ib < 6; ib++) {
        std::vector<double> entries;
        for (short block = 0; block < 16; block++) {
          double content = sp_gaps[sector][ib*16 + block];
          if (content > 0) {
            entries.push_back(content);
          }
        }
        if (entries.size() > 0) {
          fprintf(mean_outfile, ", %f", TMath::Mean(entries.begin(), entries.end()));
          fprintf(sigma_outfile, ", %f", TMath::StdDev(entries.begin(), entries.end()));
        } else {
          fprintf(mean_outfile, ", ");
          fprintf(sigma_outfile, ", ");
        }
      }
    }
    fclose(mean_outfile);
    fclose(sigma_outfile);
  }

  return sp_gaps;
}

// TODO: how to handle channels with zero mpv, just a for loop to skip the edges, probably

/**
 * @brief Write "database" to file as csv.
 * 
 * @param drop_low_rap_edge whether to drop low rapidity edge like all other edges (TRUE, better for calibration)
 *    or keep it (FALSE, default behavior of h_allblocks).
 */
void write_map_to_file(bool drop_low_rap_edge) {
  FILE *outfile = fopen("files/dbn_mpv.csv", "w+");
  fprintf(outfile, "sector, block, dbn, mpv, mpv_err, ch0_mpv, ch0_mpv_err, ch1_mpv, ch1_mpv_err, ch2_mpv, ch2_mpv_err, ch3_mpv, ch3_mpv_err");
  auto dbns = get_dbns();
  // auto mpvs = get_mpvs();
  // auto mpv_errs = get_mpv_errs(mpvs);
  std::set<int> perimeter = perimeter_channels(drop_low_rap_edge);
  auto chnl_mpv_and_err = get_chnl_mpv_with_err();
  auto block_mpv_and_err = calculate_block_mpv_with_err(chnl_mpv_and_err, perimeter);
  auto chnl_mpvs = chnl_mpv_and_err.first;
  auto chnl_mpv_errs = chnl_mpv_and_err.second;
  auto block_mpvs = block_mpv_and_err.first;
  auto block_mpv_errs = block_mpv_and_err.second;
  for (int sector = 0; sector < 64; sector++) {
    int n_blocks = 0;
    for (int block = 0; block < 96; block++) {
      std::string dbn = dbns[sector][block];
      double mpv = block_mpvs[sector][block];
      double mpv_err = block_mpv_errs[sector][block];
      auto all_chnls = block_to_channel(block);
      // double ch0_mpv = chnl_mpvs[sector][chnls[0]];
      // double ch0_mpv_err = chnl_mpv_errs[sector][chnls[0]];
      // double ch1_mpv = chnl_mpvs[sector][chnls[1]];
      // double ch1_mpv_err = chnl_mpv_errs[sector][chnls[1]];
      // double ch2_mpv = chnl_mpvs[sector][chnls[2]];
      // double ch2_mpv_err = chnl_mpv_errs[sector][chnls[2]];
      // double ch3_mpv = chnl_mpvs[sector][chnls[3]];
      // double ch3_mpv_err = chnl_mpv_errs[sector][chnls[3]];
      if (dbn != "") {
        if (mpv > 0) {
          fprintf(outfile, "\n%i, %i, %s, %f, %f", sector + 1, block + 1, dbn.c_str(), mpv, mpv_err);
          for (int &chnl : all_chnls) { // order : 0, 1, 2, 3
            if (perimeter.find(chnl) == perimeter.end()) {
              // contributes to block mpv
              fprintf(outfile, ", %f, %f", chnl_mpvs[sector][chnl], chnl_mpv_errs[sector][chnl]);
            } else {
              // does not contribute to block mpv
              fprintf(outfile, ", , ");
            }
          }
          n_blocks++;
        } else {
          fprintf(outfile, "\n%i, %i, %s, , , , , , , , , , ", sector + 1, block + 1, dbn.c_str());
        }
      }
    }
    if (debug) {
      if (n_blocks == 96) {
        printf("sector %2d: got (dbn, mpv) for all blocks!\n", sector + 1);
      } else if (n_blocks > 0) {
        printf("sector %2d: got (dbn, mpv) for %2d/96 blocks\n", sector + 1, n_blocks);
      }
    }
  }
  fclose(outfile);
}