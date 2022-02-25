#include "mpv_dbn.h"

#define DEBUG 0

std::vector<std::pair<int, int>> read_physics_runs() {
  std::vector<std::pair<int, int>> new_sector_runs(64);
  std::fstream runs_file;
  runs_file.open("files/physics_runs.csv", std::ios::in);
  std::pair<int, int> *row;
  std::string line, word;
  int line_num = 1;
  std::cout << "reading sector map" << std::endl;
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
    new_sector_runs[line_num - 2] = std::make_pair(sector, run);
    line_num++;
  }
  return new_sector_runs;
}

void get_physics_runs() {
  std::vector<std::pair<int, int>> sector_runs = read_physics_runs();
  for (int i = 0; i < 64; i++) {
    int sector = sector_runs[i].first;
    int run = sector_runs[i].second;
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

std::vector<std::vector<std::string>> get_dbns() {
  std::fstream sector_map_file;
  sector_map_file.open("files/Blocks database - Sectors.csv", std::ios::in);
  std::vector<std::vector<std::string>> all_data(24);
  std::vector<std::string> *row;
  std::string line, word;
  int line_num = 1;
  std::cout << "reading sector map" << std::endl;
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
  return dbns;
}

std::vector<std::vector<double>> get_mpvs() {
  std::vector<std::vector<double>> mpvs(64);
  for (auto &vec : mpvs) {
    vec = std::vector<double>(96, -1.0);
  }

  TFile *hist_file;
  char filename[64];
  for (std::pair<int, int> p : read_physics_runs()) {
    int sector = p.first;
    int run_num = p.second;
    if (run_num < 0) {
      // use tim's old data
      sprintf(filename, "./sector_data/sector%i.root", sector);
      hist_file = TFile::Open(filename);
      if ((hist_file = TFile::Open(filename))) {
        printf("found file for sector %i: sector_data/sector%i.root\n", sector, sector);
        TH1D* data;
        hist_file->GetObject(Form("h_run%i_block;1", sector), data);
        if (!data) {
          std::cerr << "  unable to get histogram" << std::endl;
          continue;       
        }
        for (int block_num = 1; block_num <= 96; block_num++) {
          double content = data->GetBinContent(block_num);
          // first check if this is data we are interested in...
          if (content <= 0 || content >= 1000) {
            if (DEBUG) std::cout << "  block " << block_num << ": rejected bin content " << content << std::endl;
            continue;
          } else {
            if (DEBUG) std::cout << "  block " << block_num << ": mpv " << content << std::endl;
          }
          //std::cout << "block " << block_num << " (DBN " << dbns[sector - 1][block_num] << "): good mpv (" << content  << ")" << std::endl;
          if (mpvs[sector - 1][block_num - 1] != -1) {
            throw std::runtime_error(Form("would overwrite mpv data at sector %i block %i (is there a repeated sector?)", sector, block_num));
          }
          mpvs[sector - 1][block_num - 1] = content;
        }
      } else {
        printf("FAILED to find sector %d file (%s)\n", sector, filename);
      }
    } else if (run_num > 0) {
      // get data from the run number
      sprintf(filename, "physics_runs/qa_output_000%i/histograms.root", run_num);
      if ((hist_file = TFile::Open(filename))) {
        printf("found run file for sector %i: physics_runs/qa_output_000%i/histograms.root\n", sector, run_num);
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
            if (DEBUG) std::cout << "  block " << block_num << ": rejected bin content " << content << std::endl;
            continue;
          } else {
            if (DEBUG) std::cout << "  block " << block_num << ": mpv " << content << std::endl;
          }
          //std::cout << "block " << block_num << " (DBN " << dbns[sector - 1][block_num] << "): good mpv (" << content  << ")" << std::endl;
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
  return mpvs;
}

std::map<std::string, double> get_map() {
  std::map<std::string, double> dbn_mpv_map;
  auto dbns = get_dbns();
  auto mpvs = get_mpvs();
  for (short sector = 0; sector < 64; sector++) {
    for (short block = 0; block < 96; block++) {
      std::string dbn = dbns[sector][block];
      double mpv = mpvs[sector][block];
      if (dbn != "" && mpv > 0) {
        if (dbn_mpv_map.find(dbn) != dbn_mpv_map.end()) {
          throw std::runtime_error(Form("duplicate dbn: %s (sector %i block %i)", dbn.c_str(), sector + 1, block + 1));
        } else {
         dbn_mpv_map[dbn] = mpv; 
        }
      }
    }
  }
  return dbn_mpv_map;
}

void write_map_to_file() {
  FILE *outfile = fopen("files/dbn_mpv.csv", "w+");
  fprintf(outfile, "DBN, MPV");
  auto dbns = get_dbns();
  auto mpvs = get_mpvs();
  for (short sector = 0; sector < 64; sector++) {
    int n_blocks = 0;
    for (short block = 0; block < 96; block++) {
      std::string dbn = dbns[sector][block];
      double mpv = mpvs[sector][block];
      if (dbn != "" && mpv > 0) {
        fprintf(outfile, "\n%s, %f", dbn.c_str(), mpv);
        n_blocks++;
      }
    }
    if (n_blocks > 0) {
      printf("sector %2d: got dbn, mpv for %2d/96 blocks\n", sector, n_blocks);
    }
  }
  fclose(outfile);
}