#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>

#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>

#include "../includes/mpv_dbn.h"

/**
 * @brief Struct that packages all information associated with a block.
 */
typedef struct Block {
  unsigned int sector;
  unsigned int block_number;
  std::string dbn;
  double mpv;
  double mpv_err;
  double density;
  double fiber_count;
  double fiber_t1_count;
  double fiber_t2_count;
  double fiber_t3_count;
  double fiber_t4_count;
  double scint_ratio;
} Block;

typedef std::function<std::pair<bool, double>(Block)> value_getter;

/**
 * @brief Struct that packages all details of a single plot (e.g., density).
 */
typedef struct PlotConfig {
  std::string file_name;
  std::string title;
  std::string units;
  value_getter get_value;
  double plot_min;
} PlotConfig;

/**
 * @brief Get the (x, y) location of a block within the EMCal plot.
 * 
 * @param block 
 * @return std::pair<unsigned int, unsigned int> (x, y), zero-based.
 */
std::pair<unsigned int, unsigned int> get_plot_indices(Block block) {
  bool is_north = block.sector % 2 == 0;
  unsigned int x_offset = (block.block_number - 1) % 4;
  if (!is_north) {
    x_offset = 3 - x_offset;
  }
  x_offset += 4*((block.sector - 1) / 2);
  unsigned int y_offset;
  unsigned int y_idx = (block.block_number - 1) / 4;
  if (is_north) {
    y_offset = 23 - y_idx;
  } else {
    y_offset = 24 + y_idx;
  }
  return std::make_pair(x_offset, y_offset);
}

/**
 * @brief Plot sector numbers (1 - 64) and block type (1 - 24) onto the current canvas.
 */
void plot_sector_and_block_labels(bool channel_lvl = false) {
  double sector_box_width = 3;
  double sector_box_height = 1.5;
  // odd sectors
  for (unsigned int i = 0; i < 32; i++) {
    double x_center = 2 + 4*i;
    double y_center = 48 + 1.5;
    TPaveText *text;
    if (!channel_lvl) {
      text = new TPaveText(x_center - sector_box_width/2, y_center - sector_box_height/2, x_center + sector_box_width/2, y_center + sector_box_height/2, "NB");
    } else {
      text = new TPaveText(2*(x_center - sector_box_width/2), 2*(y_center - sector_box_height/2), 2*(x_center + sector_box_width/2), 2*(y_center + sector_box_height/2), "NB");
    }
    text->SetTextFont(42);
    text->SetFillColorAlpha(0, 0);
    text->AddText(Form("%d", 2*i + 1));
    text->SetTextAlign(22);
    text->Draw();
  }
  // even sectors
  for (unsigned int i = 0; i < 32; i++) {
    double x_center = 2 + 4*i;
    double y_center = -1.5;
    TPaveText *text;
    if (!channel_lvl) {
      text = new TPaveText(x_center - sector_box_width/2, y_center - sector_box_height/2, x_center + sector_box_width/2, y_center + sector_box_height/2, "NB");
    } else {
      text = new TPaveText(2*(x_center - sector_box_width/2), 2*(y_center - sector_box_height/2), 2*(x_center + sector_box_width/2), 2*(y_center + sector_box_height/2), "NB");
    }
    text->SetTextFont(42);
    text->SetFillColorAlpha(0, 0);
    text->AddText(Form("%d", 2*i + 2));
    text->SetTextAlign(22);
    text->Draw();
  }

  double block_box_width = 3;
  double block_box_height = 1;
  for (unsigned int i = 0; i < 24; i++) {
    double x_center1 = -2;
    double y_center1 = 23.5 - i;
    TPaveText *text1;
    if (!channel_lvl) {
      text1 = new TPaveText(x_center1 - block_box_width/2, y_center1 - block_box_height/2, x_center1 + block_box_width/2, y_center1 + block_box_height/2, "NB");
    } else {
      text1 = new TPaveText(2*(x_center1 - block_box_width/2), 2*(y_center1 - block_box_height/2), 2*(x_center1 + block_box_width/2), 2*(y_center1 + block_box_height/2), "NB");
    }
    
    text1->SetTextFont(42);
    text1->SetFillColorAlpha(0, 0);
    text1->AddText(Form("%d", i + 1));
    text1->SetTextAlign(32);
    text1->Draw();

    double x_center2 = -2;
    double y_center2 = 24.5 + i;
    TPaveText *text2;
    if (!channel_lvl) {
      text2 = new TPaveText(x_center2 - block_box_width/2, y_center2 - block_box_height/2, x_center2 + block_box_width/2, y_center2 + block_box_height/2, "NB");
    } else {
      text2 = new TPaveText(2*(x_center2 - block_box_width/2), 2*(y_center2 - block_box_height/2), 2*(x_center2 + block_box_width/2), 2*(y_center2 + block_box_height/2), "NB");
    }
    text2->SetTextFont(42);
    text2->SetFillColorAlpha(0, 0);
    text2->AddText(Form("%d", i + 1));
    text2->SetTextAlign(32);
    text2->Draw();
  }

  double label_box_width = 1;
  double label_box_height = 0.8;

  // south ABCD
  std::vector<std::string> south_labels = {"A", "B", "C", "D"};
  for (unsigned int i = 0; i < 32; i++) {
    double x_start = 4*i + 0.5;
    for (unsigned int j = 0; j < 4; j++) {
      double x_center = x_start + j;
      double y_center = 48.5;
      TPaveText *text;
      if (!channel_lvl) {
        text = new TPaveText(x_center - label_box_width/2, y_center - label_box_height/2, x_center + label_box_width/2, y_center + label_box_height/2, "NB");
      } else {
        text = new TPaveText(2*(x_center - label_box_width/2), 2*(y_center - label_box_height/2), 2*(x_center + label_box_width/2), 2*(y_center + label_box_height/2), "NB");
      }
      text->SetTextFont(82);
      text->SetFillColorAlpha(0, 0);
      text->AddText(south_labels[j].c_str());
      text->SetTextAlign(22);
      text->Draw();
    }
  }
  // north ABCD
  std::vector<std::string> north_labels = {"D", "C", "B", "A"};
  for (unsigned int i = 0; i < 32; i++) {
    double x_start = 4*i + 0.5;
    for (unsigned int j = 0; j < 4; j++) {
      double x_center = x_start + j;
      double y_center = -0.5;
      TPaveText *text;
      if (!channel_lvl) {
        text = new TPaveText(x_center - label_box_width/2, y_center - label_box_height/2, x_center + label_box_width/2, y_center + label_box_height/2, "NB");
      } else {
        text = new TPaveText(2*(x_center - label_box_width/2), 2*(y_center - label_box_height/2), 2*(x_center + label_box_width/2), 2*(y_center + label_box_height/2), "NB");
      }
      text->SetTextFont(82);
      text->SetFillColorAlpha(0, 0);
      text->AddText(north_labels[j].c_str());
      text->SetTextAlign(22);
      text->Draw();
    }
  }
}

/**
 * @brief Plots EMCal plots for a single value (e.g., density).
 * 
 * @param all_blocks all blocks in the EMCal.
 * @param cfg plot configuration for the value to plot.
 */
void plot_helper(std::vector<Block> all_blocks, PlotConfig cfg) {
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.5);
  
  TH2D* h = new TH2D("", "", 128, 0, 128, 48, 0, 48);
  TH2D* h_psr = new TH2D("", "", 128, -0.1, 0.1, 48, 0, 48);

  double min_value = std::numeric_limits<double>::infinity();
  for (const Block& block : all_blocks) {
    auto offsets = get_plot_indices(block);
    unsigned int x_offset = offsets.first;
    unsigned int y_offset = offsets.second;
    std::pair<bool, double> p = cfg.get_value(block);
    if (p.first) {
      double value = p.second;
      h->SetBinContent(x_offset + 1, y_offset + 1, value);
      h_psr->SetBinContent(x_offset + 1, y_offset + 1, value);
      if (value < min_value) {
        min_value = value;
      }
      // printf("set bin content for %3d, %3d = sector %2d, block_num %2d\n", x_offset, y_offset, block.sector, block.block_number);
    } else {
      // printf("SKIPPED bin content for %3d, %3d\n", x_offset, y_offset);
    }
  }
  printf("min value for %s was %f\n", cfg.file_name.c_str(), min_value);

  std::string title = Form("sPHENIX EMCAL %s;#phi [Blocks];#eta [Blocks];%s%s", cfg.title.c_str(), cfg.title.c_str(), cfg.units.c_str());

  TCanvas* c0 = new TCanvas("c0", "", 700, 500);
  c0->SetRightMargin(0.125);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  h->SetTitle(title.c_str());
  // h->SetMinimum(min_value);
  h->SetMinimum(cfg.plot_min);
  h->SetAxisRange(0, 128 - 1, "X");
  h->SetAxisRange(0, 48 - 1, "Y");
  // h->GetXaxis()->SetRange(0, 128);
  // h->GetYaxis()->SetRange(0, 48);
  h->GetXaxis()->SetNdivisions(32, false);
  h->GetYaxis()->SetNdivisions(2, false);
  h->GetXaxis()->SetLabelOffset(999.0);
  h->GetYaxis()->SetLabelOffset(999.0);
  h->GetXaxis()->SetTickLength(0);
  h->GetYaxis()->SetTickLength(0);
  h->Draw("COLZ0");
  plot_sector_and_block_labels();
  c0->SaveAs(Form("emcal_plots/plot_colz_%s.pdf", cfg.file_name.c_str()));

  TCanvas* c1 = new TCanvas("c1", "", 1200, 500);
  c1->SetRightMargin(0.125);
  c1->SetGrid();
  gStyle->SetOptStat(0);
  h->SetTitle(title.c_str());
  // h->SetMinimum(min_value);
  h->SetMinimum(cfg.plot_min);
  h->SetAxisRange(0, 128 - 1, "X");
  h->SetAxisRange(0, 48 - 1, "Y");
  // h->GetXaxis()->SetRange(0, 128);
  // h->GetYaxis()->SetRange(0, 48);
  h->GetXaxis()->SetNdivisions(32, false);
  h->GetYaxis()->SetNdivisions(2, false);
  h->GetXaxis()->SetLabelOffset(999.0);
  h->GetYaxis()->SetLabelOffset(999.0);
  h->GetXaxis()->SetTickLength(0);
  h->GetYaxis()->SetTickLength(0);
  h->Draw("COLZ0");

  plot_sector_and_block_labels();

  // dbns
  unsigned int max_strlen = 0;
  for (const Block& block : all_blocks) {
    auto offsets = get_plot_indices(block);
    double x_center = offsets.first + 0.5;
    double y_center = offsets.second + 0.5;
    TPaveText *text = new TPaveText(x_center - 2.5, y_center - 2.5, x_center + 2.5, y_center + 2.5, "NB");
    text->SetTextFont(102);
    text->SetTextSize(0.0035);
    text->SetFillColorAlpha(0, 0);
    text->AddText(Form("%s", block.dbn.c_str()));
    text->SetTextAlign(22);
    text->Draw();
  }

  c1->SaveAs(Form("emcal_plots/plot_colz_dbn_%s.pdf", cfg.file_name.c_str()));
  
  gStyle->SetImageScaling(2.0);
  // h->SetMinimum(0);
  
  h->GetZaxis()->SetTitleOffset(1.2);

  TCanvas* c2 = new TCanvas("c2", "", 900, 500);
  h->SetLineColorAlpha(kBlack, 1.0);
  h->SetTitle(title.c_str());
  h->Draw("LEGO2 0");
  c2->SaveAs(Form("emcal_plots/plot_lego_%s.png", cfg.file_name.c_str()));

  TCanvas* c3 = new TCanvas("c3", "", 900, 500);
  h->SetLineColorAlpha(kBlack, 1.0);
  h->SetFillColorAlpha(kBlue, 1.0);
  h->SetTitle(title.c_str());
  h->Draw("LEGO1 CYL 0");
  c3->SaveAs(Form("emcal_plots/plot_cyl_%s.png", cfg.file_name.c_str()));
  
  TCanvas* c4 = new TCanvas("c4", "", 900, 500);
  h_psr->SetLineColorAlpha(kBlack, 1.0);
  h_psr->SetFillColorAlpha(kBlue, 1.0);
  h_psr->SetTitle(title.c_str());
  h_psr->Draw("LEGO1 PSR 0");
  c4->SaveAs(Form("emcal_plots/plot_psr_%s.png", cfg.file_name.c_str()));
  
  delete c0;
  delete c1;
  delete c2;
  delete c3;
  delete c4;

  gStyle->SetImageScaling(1.0);
}

// TODO: check that the channel mapping is correct and matches Tim's (treat north and south separately)
void plot_channel_lvl(std::vector<Block> all_blocks) {
  std::vector<std::vector<double>> chnl_mpvs = get_chnl_mpv();
  TH2D *h_chnl_mpv = new TH2D("", "", 256, 0, 256, 96, 0, 96);
  for (Block &block : all_blocks) {
    auto xy = get_plot_indices(block);
    unsigned int x = 2*xy.first + 1;
    unsigned int y = 2*xy.second + 1;

    // the first block in a sector is in column D and is on the low rapidity edge
    //    and also has only contributions from channel 

    // in Tim's html, we are looking down at the narrow end of the block
    // in Tim's html, blocks are arranged by 1-based number as below:
    // D 01 05 09 ...
    // C 02 06 10 ...
    // B 03 07 11 ...
    // A 04 08 12 ...
    // in Tim's html, channels are arranged within each block as below:
    //  2 3
    //  0 1
    // so, low rapidity edge has 0 and 2 and 2 -> 0 points in the direction of increasing block number
    // so, (along long edge of sector towards low rapidity) cross (along short edge of sector towards increasing block number, or D -> A) 
    //    gives direction of the inside of the detector (wide -> narrow)
    // so, here, we are looking down on the outside of the detector (wide ends of blocks)

    // here, top half (higher y) is South, where LHS of sector is column A and bottom is low rapidity (block 123)
    //    bottom half (lower y) is North, where LHS of sectors is column D and top is low rapidity (block 123)    
    // here, on the top half (South), blocks are arranged by 1-based number as below:
    // .  .  .  .
    // .  .  .  .
    // .  .  .  .
    // 12 11 10 09
    // 08 07 06 05
    // 04 03 02 01
    // A  B  C  D
    // here, on the top half (South), channels are arranged within each block as below:
    //  1 3
    //  0 2
    // here, on the bottom half (North), blocks are arranged by 1-based number as below:
    // D  C  B  A
    // 01 02 03 04
    // 05 06 07 08
    // 09 10 11 12
    // .  .  .  .
    // .  .  .  .
    // .  .  .  .
    // here, on the bottom half (North), channels are arranged within each block as below:
    //  2 0
    //  3 1

    // Tim's    This N    This S
    //  2 3      2 0       1 3
    //  0 1      3 1       0 2
    auto chnls = block_to_channel(block.block_number - 1);
    double ch0 = chnl_mpvs[block.sector - 1][chnls[0]];
    double ch1 = chnl_mpvs[block.sector - 1][chnls[1]];
    double ch2 = chnl_mpvs[block.sector - 1][chnls[2]];
    double ch3 = chnl_mpvs[block.sector - 1][chnls[3]];
    if (block.sector % 2 == 0) { // north sector
      if (ch0 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y + 1, ch0);
      }
      if (ch1 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y, ch1);
      }
      if (ch2 > 0) {
        h_chnl_mpv->SetBinContent(x, y + 1, ch2);
      }
      if (ch3 > 0) {
        h_chnl_mpv->SetBinContent(x, y, ch3);
      }
    } else { // south sector
      if (ch0 > 0) {
        h_chnl_mpv->SetBinContent(x, y, ch0);
      }
      if (ch1 > 0) {
        h_chnl_mpv->SetBinContent(x, y + 1, ch1);
      }
      if (ch2 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y, ch2);
      }
      if (ch3 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y + 1, ch3);
      }
    }
    
    // printf("sector %2d block %2d:\n\tch0: %f\n\tch1: %f\n\tch2: %f\n\tch3: %f\n", block.sector, block.block_number, ch0, ch1, ch2, ch3);
  }
  TCanvas *c_chnl_mpv = new TCanvas();
  gStyle->SetOptStat(0);
  h_chnl_mpv->GetXaxis()->SetNdivisions(32, false);
  h_chnl_mpv->GetYaxis()->SetNdivisions(2, false);
  h_chnl_mpv->GetXaxis()->SetLabelOffset(999.0);
  h_chnl_mpv->GetYaxis()->SetLabelOffset(999.0);
  h_chnl_mpv->GetXaxis()->SetTickLength(0);
  h_chnl_mpv->GetYaxis()->SetTickLength(0);
  h_chnl_mpv->Draw("COLZ0");
  plot_sector_and_block_labels(true);
  c_chnl_mpv->SaveAs("emcal_plots/chnl_mpv.pdf");

  // the normal of the top surface of block points to high rapidity 
  // 
  // *looking at the wide end of the block*
  //     top
  //  ---------
  //  \  3 4  /
  //   \ 1 2 /
  //    -----
  //   bottom
  
  // (wide -> narrow) cross (high -> low rapidity) gives direction of increasing block number
  //    and as before, direction of increasing block number is direction of 2 -> 0 of channels within block 
  // in the figure above, wide -> narrow is into the page and high -> low rapidity is down, so
  //    increasing block number is to the left
  // this means we have the following mapping:
  // channel <-> fiber tower
  //       0  |  1
  //       2  |  2
  //       1  |  3
  //       3  |  4

  // in Tim's html, fiber towers are arranged as below:
  //  2 4
  //  1 3
  // here, on the top half (South), fiber towers are arranged within each block as below:
  //  3 4
  //  1 2
  // here, on the bottom half (North), fiber towers are arranged within each block as below:
  //  2 1
  //  4 3

  TH2D *h_chnl_fiber = new TH2D("", "", 256, 0, 256, 96, 0, 96);
  for (Block &block : all_blocks) {
    auto xy = get_plot_indices(block);
    unsigned int x = 2*xy.first + 1;
    unsigned int y = 2*xy.second + 1;

    double t1 = block.fiber_t1_count;
    double t2 = block.fiber_t2_count;
    double t3 = block.fiber_t3_count;
    double t4 = block.fiber_t4_count;
    if (block.sector % 2 == 0) { // north sector
      if (t1 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y + 1, t1);
      }
      if (t2 > 0) {
        h_chnl_fiber->SetBinContent(x, y + 1, t2);
      }
      if (t3 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y, t3);
      }
      if (t4 > 0) {
        h_chnl_fiber->SetBinContent(x, y, t4);
      }
    } else { // south sector
      if (t1 > 0) {
        h_chnl_fiber->SetBinContent(x, y, t1);
      }
      if (t2 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y, t2);
      }
      if (t3 > 0) {
        h_chnl_fiber->SetBinContent(x, y + 1, t3);
      }
      if (t4 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y + 1, t4);
      }
    }
    
    // printf("sector %2d block %2d:\n\tch0: %f\n\tch1: %f\n\tch2: %f\n\tch3: %f\n", block.sector, block.block_number, ch0, ch1, ch2, ch3);
  }
  TCanvas *c_chnl_fiber = new TCanvas();
  gStyle->SetOptStat(0);
  h_chnl_fiber->GetXaxis()->SetNdivisions(32, false);
  h_chnl_fiber->GetYaxis()->SetNdivisions(2, false);
  h_chnl_fiber->GetXaxis()->SetLabelOffset(999.0);
  h_chnl_fiber->GetYaxis()->SetLabelOffset(999.0);
  h_chnl_fiber->GetXaxis()->SetTickLength(0);
  h_chnl_fiber->GetYaxis()->SetTickLength(0);
  h_chnl_fiber->Draw("COLZ0");
  plot_sector_and_block_labels(true);
  c_chnl_fiber->SaveAs("emcal_plots/chnl_fiber_count.pdf");
}

/**
 * @brief Body of macro (called when macro is executed). 
 */
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
    double mpv_err;
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
    if (tmp == "") {
      mpv = -1;
    } else {
      std::istringstream(tmp) >> mpv;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      mpv_err = -1;
    } else {
      std::istringstream(tmp) >> mpv_err;
    }

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


    // dbn.erase(std::remove(dbn.begin(), dbn.end(), '\n'), dbn.end());
    // dbn.erase(std::remove(dbn.begin(), dbn.end(), '\r'), dbn.end());
    // dbn.erase(std::remove(dbn.begin(), dbn.end(), ' '), dbn.end());

    all_blocks.push_back({sector, block_number, dbn, mpv, mpv_err, density, fiber_count, fiber_t1_count, fiber_t2_count, fiber_t3_count, fiber_t4_count, scint_ratio});
    // printf("sector %2d block %2d mpv %f +/- %f\n", sector, block_number, mpv, mpv_err);
    line_num++;
  }

  std::vector<PlotConfig> cfgs = {
    {
      "mpv", "MPV", "", [](Block block){
        double value = block.mpv;
        return std::make_pair(value >= 0, value);
      },
      0.0
    },
    {
      "scint_ratio", "Scintillation Ratio", "", [](Block block){
        double value = block.scint_ratio;
        return std::make_pair(value >= 0, value);
      },
      0.7
    },
    {
      "fiber_count", "Fiber Count", " [%]", [](Block block){
        double value = block.fiber_count;
        return std::make_pair(value >= 0, value);
      },
      95.0
    },
    {
      "density", "Density", " [g/mL]", [](Block block){
        double value = block.density;
        return std::make_pair(value >= 0, value);
      },
      8.4
    }
  };

  plot_channel_lvl(all_blocks);
  // for (const PlotConfig& cfg : cfgs) {
  //   plot_helper(all_blocks, cfg);
  // }
}
