#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <filesystem>
#include <stdexcept>

#include <TSystem.h>
#include <TH1D.h>
#include <TFile.h>

/**
 * @brief gives the run number (5-digit) for each sector, or -1
 * if data is in the old sector files from Tim
 * 
 */
std::vector<std::pair<int, int>> new_sector_runs = {
  {1, -1},
  {2, -1},
  {3, -1},
  {4, -1},
  {5, -1},
  {6, -1},
  {7, -1},
  {8, -1},
  {9, -1},
  {10, 16758},
  {11, -1},
  {12, -1},
  {13, 16798},
  {14, -1},
  {15, -1},
  {16, 16878},
  {17, -1},
  {18, -1},
  {19, 16945},
  {20, 16163},
  {21, 16231},
  {22, 16366},
  {23, 16720},
  {24, 16535},
  {25, 17063},
  {26, 17179},
  {27, 17262},
  {28, 17398},
  {29, 17465},
  {30, 17537},
  {31, 17667},
  {32, 17804},
  {32, 17867}
};

std::vector<std::vector<std::string>> get_dbns();
std::vector<std::vector<double>> get_mpvs();
std::map<std::string, double> get_map();
void write_map_to_file();

#include "mpv_dbn.cpp"