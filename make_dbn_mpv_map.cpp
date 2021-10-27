#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <filesystem>

#include <assert.h>
#include <string.h>
#include <sys/stat.h>

#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>


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

void make_dbn_mpv_map () {
  std::fstream sector_map_file;
  sector_map_file.open("Blocks database - Sectors.csv",std::ios::in);
  std::map<int, std::vector<std::string>> sector_map;
  std::vector<std::vector<std::string>> all_data(24);
  std::vector<std::string> *row;
  std::string line, word;
  int line_num = 1;
  int sector_num = -1;
  int this_row_num = -1;
  int offset = -1;
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
              throw std::runtime_error("expected dbn to be F(int)");
            }
          } else if (dbn[0] == 'C') {
            // expect that the rest is castable to int
            int res = sscanf(dbn.c_str(), "C%i", &d);
            if (res != 1) {
              throw std::runtime_error("expected dbn to be C(int)");
            }
          } else {
            // expect dbn is castable to int
            int res = sscanf(dbn.c_str(), "%i", &d);
            if (res != 1) {
              throw std::runtime_error("expected dbn to be (int)");
            }
          }
          
          
          dbns[i].push_back(dbn);
        } 
      }
    }
  }




  // now get mpvs
  std::vector<std::vector<double>> mpvs(64);
  for (auto &vec : mpvs) {
    vec = std::vector<double>(96, -1);
  }

  const char *inDir = "./sector_data";
  char *dir = gSystem->ExpandPathName(inDir);
  void *dirp = gSystem->OpenDirectory(dir);

  const char *entry;
  const char *filename;
  Int_t n = 0;
  char const *start = "sector";
  TFile *hist_file;

  

  while ((entry = (char*) gSystem->GetDirEntry(dirp))) {
    filename = gSystem->ConcatFileName(dir, entry);
    TFile *hist_file;
    int sector = 0;
    char str[64];
    int res = sscanf(entry, "sector%i%s", &sector, str);
    if (res == 2 && !strcmp(str, ".root") && (hist_file = TFile::Open(filename))) {
      printf("found sector file %s ...\n", entry);
      TH1D* data;
      hist_file->GetObject(Form("h_run%i_block;1", sector), data);
      if (!data) {
        std::cerr << "  unable to get histogram" << std::endl;
        continue;       
      }
      for (int i = 0; i < data->GetNcells(); i++) {
        double content = data->GetBinContent(i);
        int block_num = data->GetBinLowEdge(i);
        // first check if this is data we are interested in...
        if (block_num < 0 || block_num >= 96) {
          std::cout << "  block " << block_num << " was out of range" << std::endl;
          continue;
        } else if (content <= 0 || content >= 1000) {
          std::cout << "  block " << block_num << ": rejected bin content " << content << std::endl;
          continue;
        }
        //std::cout << "block " << block_num << " (DBN " << dbns[sector - 1][block_num] << "): good mpv (" << content  << ")" << std::endl;
        mpvs[sector - 1][block_num] = content;
      }
    }
  }

  // and write to csv file
  FILE *outfile = fopen("dbn_mpv_mpv.csv", "w+");
  fprintf(outfile, "DBN, MPV");
  for (short sector = 0; sector < 64; sector++) {
    for (short block = 0; block < 96; block++) {
      std::string dbn = dbns[sector][block];
      double mpv = mpvs[sector][block];
      if (dbn != "" && mpv > 0) {
        fprintf(outfile, "\n%s, %f", dbn.c_str(), mpv);
      }
    }
  }
  fclose(outfile);
}