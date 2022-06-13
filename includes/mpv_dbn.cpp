#include "mpv_dbn.h"
#include <TMath.h>

constexpr bool debug = false;

/**
 * @brief Get run numbers for each sector from csv.

 * NOTE: depends on files/physics_runs.csv.
 * 
 * @return std::map<int, int> (sector, run number).
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
 * @return std::vector<std::vector<std::string>> sector[block number] -> dbn, or "" if none.
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
 * @brief Get MPVs for each block for each sector.
 * 
 * NOTE: rejects MPVs <= 0 and >= 1000.
 * 
 * @param write_ib write IB mean and sigma to csv file.
 * @return std::vector<std::vector<double>> sector[block number] -> mpv, or -1 if none.
 */
std::vector<std::vector<double>> get_mpvs(bool write_ib = false) {
  std::vector<std::vector<double>> mpvs(64);
  for (auto &vec : mpvs) {
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
      sprintf(filename, "physics_runs/qa_output_000%i/histograms.root ", run_num);
      if ((hist_file = TFile::Open(filename))) {
        if (debug) {
          printf("found run file for sector %i: physics_runs/qa_output_000%i/histograms.root\n", sector, run_num);
        }
        TH1D* data;
        hist_file->GetObject("h_allblocks;1", data);
        if (!data) {
          std::cerr << "  unable to get histogram" << std::endl;
          continue;       
        }
        for (int block_num = 1; block_num <= 96; block_num++) {
          double content = data->GetBinContent(block_num);
          // first check if this is data we are interested in...
          if (content <= 0 || content >= 1000) {
            if (debug) std::cout << "  block " << block_num << ": rejected bin content " << content << std::endl;
            continue;
          } else {
            // if (debug) std::cout << "  block " << block_num << ": mpv " << content << std::endl;
          }
          // std::cout << "block " << block_num << " (DBN " << dbns[sector - 1][block_num] << "): good mpv (" << content  << ")" << std::endl;
          if (mpvs[sector - 1][block_num - 1] != -1) {
            throw std::runtime_error(Form("would overwrite mpv data at sector %i block %i (is there a repeated sector?)", sector, block_num));
          }
          mpvs[sector - 1][block_num - 1] = content;
        }
      } else {
        printf("FAILED to find run file for sector %i: qa_output_000%i/histograms.root\n", sector, run_num);
      }
    }
  }
  if (write_ib) {
    FILE *mean_outfile = fopen("files/mpv_avg_ib_mean.csv", "w+");
    FILE *sigma_outfile = fopen("files/mpv_avg_ib_sigma.csv", "w+");
    fprintf(mean_outfile, "SECTOR, RUN, IB0, IB1, IB2, IB3, IB4, IB5");
    fprintf(sigma_outfile, "SECTOR, RUN, IB0, IB1, IB2, IB3, IB4, IB5");
    for (short sector = 0; sector < 64; sector++) {
      fprintf(mean_outfile, "\n%i, %i", sector + 1, physics_runs[sector + 1]);
      fprintf(sigma_outfile, "\n%i, %i", sector + 1, physics_runs[sector + 1]);
      for (short ib = 0; ib < 6; ib++) {
        std::vector<double> entries;
        for (short block = 0; block < 16; block++) {
          double mpv = mpvs[sector][ib*16 + block];
          if (mpv != -1) {
            entries.push_back(mpv);
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
  for (short sector = 0; sector < 64; sector++) {
    int n_blocks = 0;
    for (short block = 0; block < 96; block++) {
      double mpv = mpvs[sector][block];
      if (mpv > 0) {
        n_blocks++;
      }
    }
    if (debug) {
      if (n_blocks == 96) {
        printf("sector %2d: got mpv for all blocks!\n", sector + 1);
      } else if (n_blocks > 0) {
        printf("sector %2d: got mpv for %2d/96 blocks\n", sector + 1, n_blocks);
      }
    }
  }
  return mpvs;
}

/**
 * @brief Get MPV errors for each block for each sector.
 * 
 * @return std::vector<std::vector<double>> sector[block number] -> mpv error.
 */
std::vector<std::vector<double>> get_mpv_errs() {
  std::vector<std::vector<double>> mpv_errs(64);
  for (auto &vec : mpv_errs) {
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
        for (int block_num = 1; block_num <= 96; block_num++) {
          double err = 0;
          for (int tower = 0; tower < 4; tower++) {
            err += data->GetBinError(4*(block_num - 1) + tower + 1);
          }
          if (mpv_errs[sector - 1][block_num - 1] != -1) {
            throw std::runtime_error(Form("would overwrite mpv_err data at sector %i block %i (is there a repeated sector?)", sector, block_num));
          }
          printf("sector %2d, block %2d, mpv_err = %f\n", sector, block_num, err);
          mpv_errs[sector - 1][block_num - 1] = err;
        }
      } else {
        printf("FAILED to find run file for sector %i: qa_output_000%i/histograms.root\n", sector, run_num);
      }
    }
  }
  // for (short sector = 0; sector < 64; sector++) {
  //   int n_blocks = 0;
  //   for (short block = 0; block < 96; block++) {
  //     double mpv = mpv_errs[sector][block];
  //     if (mpv > 0) {
  //       n_blocks++;
  //     }
  //   }
  //   if (debug) {
  //     if (n_blocks == 96) {
  //       printf("sector %2d: got mpv_err for all blocks!\n", sector + 1);
  //     } else if (n_blocks > 0) {
  //       printf("sector %2d: got mpv_err for %2d/96 blocks\n", sector + 1, n_blocks);
  //     }
  //   }
  // }
  return mpv_errs;
}

/**
 * @brief Get single pixel gaps for each block for each sector (average of block's four towers).
 * 
 * @param write_ib write IB mean and sigma to csv file.
 * @return std::vector<std::vector<double>> sector[block number] -> gap.
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
        for (int block_num = 1; block_num <= 96; block_num++) {
          double avg_sp_gap = 0;
          for (int i = 1; i <= 4; i++) {
            double content = data->GetBinContent((block_num - 1)*4 + i);
            if (content <= 0) {
              printf("complaint at sector %i channel %i: sp_gap < 0 (%f)\n", sector, (block_num - 1)*4 + i, content);
            }
            avg_sp_gap += content;
          }
          if (sp_gaps[sector - 1][block_num - 1] != -1) {
            throw std::runtime_error(Form("would overwrite sp_gap data at sector %i block %i (is there a repeated sector?)", sector, block_num));
          }
          sp_gaps[sector - 1][block_num - 1] = avg_sp_gap/4;
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

/**
 * @brief Write "database" to file as csv. 
 */
void write_map_to_file() {
  FILE *outfile = fopen("files/dbn_mpv.csv", "w+");
  fprintf(outfile, "SECTOR, BLOCK, DBN, MPV");
  auto dbns = get_dbns();
  auto mpvs = get_mpvs();
  for (short sector = 0; sector < 64; sector++) {
    int n_blocks = 0;
    for (short block = 0; block < 96; block++) {
      std::string dbn = dbns[sector][block];
      double mpv = mpvs[sector][block];
      if (dbn != "") {
        if (mpv > 0) {
          fprintf(outfile, "\n%i, %i, %s, %f", sector + 1, block + 1, dbn.c_str(), mpv);
          n_blocks++;
        } else {
          fprintf(outfile, "\n%i, %i, %s, ", sector + 1, block + 1, dbn.c_str());
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