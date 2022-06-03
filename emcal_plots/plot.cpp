#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

typedef struct Block {
  unsigned int sector;
  unsigned int block_number;
  std::string dbn;
  double mpv;
  double density;
  double fiber_count;
  double fiber_t1_count;
  double fiber_t2_count;
  double fiber_t3_count;
  double fiber_t4_count;
  double scint_ratio;
} Block;

void plot() {
  std::vector<Block> all_blocks;
  std::fstream database;
  database.open("emcal_plots/sPHENIX_EMCal_blocks - dbn_mpv.csv", std::ios::in);
  std::string line, word;
  int line_num = 1;
  std::cout << "reading 'new database'" << std::endl;
  while (std::getline(database, line)) {
    if (line_num == 1) {
      line_num++;
      continue;
    }
    std::stringstream csvStream(line);
    unsigned int sector;
    unsigned int block_number;
    std::string dbn;
    double mpv;
    double density;
    double fiber_count;
    double fiber_t1_count;
    double fiber_t2_count;
    double fiber_t3_count;
    double fiber_t4_count;
    double scint_ratio;
    
    std::string s, tmp;

    // std::getline(csvStream, s, ',');
    

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> sector;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> block_number;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> dbn;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> mpv;

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      density = -1;
    } else {
      std::istringstream(tmp) >> density;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      fiber_count = -1;
    } else {
      std::istringstream(tmp) >> fiber_count;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      fiber_t1_count = -1;
    } else {
      std::istringstream(tmp) >> fiber_t1_count;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      fiber_t2_count = -1;
    } else {
      std::istringstream(tmp) >> fiber_t2_count;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      fiber_t3_count = -1;
    } else {
      std::istringstream(tmp) >> fiber_t3_count;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      fiber_t4_count = -1;
    } else {
      std::istringstream(tmp) >> fiber_t4_count;
    }

    std::getline(csvStream, tmp);
    if (tmp == "" || tmp == "\r") {
      scint_ratio = -1;
    } else {
      std::istringstream(tmp) >> scint_ratio;
    }



    all_blocks.push_back({sector, block_number, dbn, mpv, density, fiber_count, fiber_t1_count, fiber_t2_count, fiber_t3_count, fiber_t4_count, scint_ratio});
    line_num++;
  }


  TH2D* h = new TH2D("", "", 128, 0, 128, 48, 0, 48);

  for (const Block& block : all_blocks) {
    bool is_north = block.sector % 2 == 0;
    unsigned int x_offset = block.block_number % 4;
    if (!is_north) {
      x_offset = 3 - x_offset;
    }
    x_offset += 4*((block.sector - 1) / 2);
    unsigned int y_offset;
    unsigned int y_idx = block.block_number / 4;
    if (is_north) {
      y_offset = 23 - y_idx;
    } else {
      y_offset = 24 + y_idx;
    }
    if (block.mpv > 0) {
      h->SetBinContent(x_offset, y_offset, block.mpv);
    }
  }

  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(0);
  // h->Draw("COLZ");
  h->Draw("LEGO2Z");
  c->SaveAs("emcal_plots/plot.pdf");
}
