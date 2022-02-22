#include "includes/mpv_dbn.h"

// C++ struct to hold information associated with a block
struct Block {
  bool has_mpv;
  int dbn;
  int fiber_batch;
  double mpv;
  double adj_mpv;
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

// & indicates "pass by reference" in C++, which is needed in this case
void add_fiber_batch_info(std::string file_name, std::map<int, int> &fb_map) {
  std::fstream database_sheet_file;
  database_sheet_file.open(file_name, std::ios::in);
  if (!database_sheet_file) {
    throw std::runtime_error(Form("unable to open file '%s'", file_name.c_str()));
  }
  std::vector<std::string> row;
  std::string line, word, temp;
  int line_num = 1;
  int this_row_num = -1;
  int offset = 0;
  std::cout << "reading " << file_name << std::endl;
  while (std::getline(database_sheet_file, line)) {
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

void fiber_batch_correction() {
  std::map<int, int> dbn_to_fiber_batch;
  add_fiber_batch_info("files/blocks_1-12.csv", dbn_to_fiber_batch);
  add_fiber_batch_info("files/blocks_13-64.csv", dbn_to_fiber_batch);
  auto dbns = get_dbns();
  auto mpvs = get_mpvs();

  std::map<int, double> batch_averages;
  for (auto vec : dbns) {
    for (std::string dbn : vec) {
      int batch = dbn_to_fiber_batch[dbn];
      printf("%s: %i\n");
    }
  }
  // DESIGN CHOICE: use strings or ints as DBNs (map keys)??
}