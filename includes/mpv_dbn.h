#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
//#include <filesystem>
#include <stdexcept>

#include <dirent.h>

#include <TSystem.h>
#include <TH1D.h>
#include <TFile.h>

std::map<int, int> read_physics_runs();
void get_physics_runs();
std::vector<std::vector<std::string>> get_dbns();
std::vector<std::vector<double>> get_mpvs(bool write_ib);
std::vector<std::vector<double>> get_mpv_errs();
std::vector<std::vector<double>> get_sp_gaps(bool write_ib);
void write_map_to_file();

#include "mpv_dbn.cpp"