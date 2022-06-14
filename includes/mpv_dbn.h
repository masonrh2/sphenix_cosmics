#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
//#include <filesystem>
#include <stdexcept>

#include <dirent.h>

#include <TSystem.h>
#include <TH1D.h>
#include <TFile.h>

std::set<int> perimeter_channels(bool drop_low_rap_edge);
std::map<int, int> read_physics_runs();
void get_physics_runs();
std::vector<std::vector<std::string>> get_dbns();
std::vector<std::vector<double>> get_mpvs(bool write_ib);
std::vector<std::vector<double>> get_mpv_errs(std::vector<std::vector<double>> mpvs);
std::vector<std::vector<std::pair<double, double>>> get_mpv_with_err(bool drop_low_rap_edge);
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> get_chnl_mpv();
std::vector<std::vector<double>> get_sp_gaps(bool write_ib);
void write_map_to_file();

#include "mpv_dbn.cpp"