#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include <map>

#include <TH2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TF1.h>

#include "../includes/mpv_dbn.h"
#include "../includes/utils.h"

/**
 * @brief The default (naive) sector mapping (from Caroline's sector sheet). NOTE: we are looking down on the narrow ends of blocks / inside of detector!
 */
const std::vector<int> pseudo_sector_mapping = {
   1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, // top of plot
   2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64  // bottom of plot
};

/**
 * @brief The true sector mapping as they should appear in the final plots. NOTE: we are looking down on the narrow ends of blocks / inside of detector!
 */
const std::vector<int> true_sector_mapping = {
//><phi=0 is the middle of first sector on list; increasing phi to the right -->
  12, 44, 42, 20, 58, 18, 62, 48,  8, 10, 16, 28, 50, 24, 38, 60,  6, 64, 14, 32, 36, 26, 46, 54,  4,  2, 22, 30, 56, 34, 40, 52, // top of plot = North
   5, 49, 37, 19, 51, 35, 39, 43,  3,  1, 29, 23, 59, 31, 41, 53,  9, 63, 17, 33, 57, 15, 13, 45,  7, 11, 25, 27, 61, 21, 47, 55  // bottom of plot = South
};

/**
 * @brief ...
 */
const std::map<std::string, std::string> fiber_type_compressor = {
  {"", ""},
  {"I-K", "K"},
  {"K", "K"},
  {"P-SG", "SG"},
  {"PSG+IK+K", ""},
  {"SG", "SG"},
  {"SG-B", "SG"},
  {"SG47", "SG"},
};

/**
 * @brief Check if a sector mapping is valid.
 * 
 * @param mapping 
 */
void check_sector_mapping(const std::vector<int> mapping) {
  std::set<int> sectors;
  if (mapping.size() != 64) {
    throw std::runtime_error(Form("sector mapping has %zu != 64 entries", mapping.size()));
  }
  for (const int &sector : mapping) {
    if (sectors.find(sector) != sectors.end()) {
      throw std::runtime_error(Form("found repeated sector (%d) while checking mapping", sector));
    }
    if (sector < 1 || sector > 64) {
      throw std::runtime_error(Form("found invalid sector (%d) while checking mapping", sector));
    }
  }
}

/**
 * @brief Struct that packages all information associated with a block.
 */
typedef struct Block {
  unsigned int sector;
  unsigned int block_number;
  std::string dbn;
  double mpv;
  double mpv_err;
  double ch0_mpv;
  double ch0_mpv_err;
  double ch1_mpv;
  double ch1_mpv_err;
  double ch2_mpv;
  double ch2_mpv_err;
  double ch3_mpv;
  double ch3_mpv_err;
  double density;
  double fiber_count;
  double fiber_t1_count;
  double fiber_t2_count;
  double fiber_t3_count;
  double fiber_t4_count;
  double scint_ratio;
  std::string fiber_type;
  std::string fiber_batch;
  std::string w_powder;
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
 * @param sector_mapping 
 * @return std::pair<unsigned int, unsigned int> (x, y), zero-based.
 */
std::pair<unsigned int, unsigned int> get_block_loc(Block block, std::vector<int> sector_mapping) {
  auto it = std::find(sector_mapping.begin(), sector_mapping.end(), block.sector);
  int pseudo_sector = pseudo_sector_mapping[it - sector_mapping.begin()];
  // printf("true sector %2d -> pseudo sector %2d\n", block.sector, pseudo_sector);
  bool is_top_half_of_plot = pseudo_sector % 2 == 1;
  unsigned int x_offset = (block.block_number - 1) % 4;
  if (!is_top_half_of_plot) {
    x_offset = 3 - x_offset;
  }
  x_offset += 4*((pseudo_sector - 1) / 2);
  unsigned int y_offset;
  unsigned int y_idx = (block.block_number - 1) / 4;
  if (!is_top_half_of_plot) {
    y_offset = 23 - y_idx;
  } else {
    y_offset = 24 + y_idx;
  }
  return std::make_pair(x_offset, y_offset);
}

/**
 * @brief Plot sector numbers (1 - 64) and block type (1 - 24) onto the current canvas. Draws sector numbers above/below sector.
 *    Draws block type left of y axis. NOTE: drawn numbers overlap with x and y axis number, so  one must ...->SetLabelOffset(999.0). 
 * 
 * @param channel_lvl whether plotting on a channel/tower-level canvas (TRUE) or block-level (FALSE).
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
 * @brief Plot true sector numbers (1 - 64) onto the current canvas. Draws sector numbers in the middle of sector.
 * 
 * @param channel_lvl whether plotting on a channel/tower-level canvas (TRUE) or block-level (FALSE).
 */
void plot_sector_labels_debug(bool channel_lvl = false) {
  double sector_box_width = 3;
  double sector_box_height = 1.5;
  // SOUTH (odd pseudo-sectors)
  for (unsigned int i = 0; i < 32; i++) {
    double x_center;
    x_center = 2 + 4*i;
    double y_center = 36;
    TPaveText *text;
    if (!channel_lvl) {
      text = new TPaveText(x_center - sector_box_width/2, y_center - sector_box_height/2, x_center + sector_box_width/2, y_center + sector_box_height/2, "NB");
    } else {
      text = new TPaveText(2*(x_center - sector_box_width/2), 2*(y_center - sector_box_height/2), 2*(x_center + sector_box_width/2), 2*(y_center + sector_box_height/2), "NB");
    }
    text->SetTextFont(42);
    text->SetFillColorAlpha(0, 0);
    text->AddText(Form("%d", true_sector_mapping[i]));
    text->SetTextAlign(22);
    text->Draw();
  }
  // NORTH (even pseudo-sectors)
  for (unsigned int i = 0; i < 32; i++) {
    double x_center;
    x_center = 2 + 4*i;
    double y_center = 12;
    TPaveText *text;
    if (!channel_lvl) {
      text = new TPaveText(x_center - sector_box_width/2, y_center - sector_box_height/2, x_center + sector_box_width/2, y_center + sector_box_height/2, "NB");
    } else {
      text = new TPaveText(2*(x_center - sector_box_width/2), 2*(y_center - sector_box_height/2), 2*(x_center + sector_box_width/2), 2*(y_center + sector_box_height/2), "NB");
    }
    text->SetTextFont(42);
    text->SetFillColorAlpha(0, 0);
    text->AddText(Form("%d", true_sector_mapping[i + 32]));
    text->SetTextAlign(22);
    text->Draw();
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
    auto offsets = get_block_loc(block, pseudo_sector_mapping);
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
    auto offsets = get_block_loc(block, pseudo_sector_mapping);
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

/**
 * @brief Plots EMCal fiber count over towers and MPV over channels.
 * 
 * @param all_blocks all blocks in the EMCal.
 * @param drop_low_rap_edge whether to exclude low rapidity edge like all other edges (TRUE, better for calibration)
 *    or keep it (plot it) (FALSE, default behavior of h_allblocks)
 */
void plot_channel_lvl(std::vector<Block> all_blocks, std::string mode) {
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.5);

  // in Tim's html, we are looking down at the narrow end of the block and blocks are arranged by 1-based number as below:
  // (D) 01 05 09 ...
  // (C) 02 06 10 ...
  // (B) 03 07 11 ...
  // (A) 04 08 12 ...
  // and channels are arranged within each block as below:
  //  2 3
  //  0 1
      
  // here, on the top half (North), blocks are arranged by 1-based number as below:
  // .  .  .  .
  // .  .  .  .
  // .  .  .  .
  // 09 10 11 12
  // 05 06 07 08
  // 01 02 03 04
  // D  C  B  A
  // and channels are arranged within each block as below:
  //  3 1
  //  2 0
  // here, on the bottom half (South), blocks are arranged by 1-based number as below:
  
  // A  B  C  D
  // 04 03 02 01
  // 08 07 06 05
  // 12 11 10 09
  // .  .  .  .
  // .  .  .  .
  // .  .  .  .
  // and channels are arranged within each block as below:
  //  0 2
  //  1 3

  // MPV CHANNEL MAPPING:
  // Tim's  here top/N  here bottom/S
  //  2 3      3 1          0 2
  //  0 1      2 0          1 3

  // the surface normal of the top surface of block points to high rapidity 
  // the fiber towers are mapped within a block as follows:
  //     top
  //  ---------
  //  \  3 4  /     *looking at the WIDE end of the block*
  //   \ 1 2 /
  //    -----
  //   bottom
  
  // channel <-> fiber tower
  //       0  |  1
  //       2  |  2
  //       1  |  3
  //       3  |  4

  // FIBER TOWER MAPPING:
  // Tim's  here top/N  here bottom/S
  //  2 4      4 3          1 2
  //  1 3      2 1          3 4

  gStyle->SetOptStat(1);
  TH1D *h_mpv_dist = new TH1D("h_chnl_mpv", "Distribution of EMCal Channel MPV;MPV;Count [Channels]", 80, 0, 1000);
  TH1D *h_fiber_count_dist = new TH1D("h_chnl_fiber", "Distribution of EMCal Tower Fiber Count;Fiber Count [%];Count [Towers]", 80, 80, 120);
  THStack *hs_mpv_dist = new THStack("hs_chnl_mpv", "Distribution of EMCal Channel MPV;MPV;Count [Channels]");
  TH1D *h_mpv_dist_uiuc = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_dist_china_sg = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_dist_china_k = new TH1D("", "", 80, 0, 800);
  THStack *hs_fiber_count_dist = new THStack("hs_chnl_fiber", "Distribution of EMCal Tower Fiber Count;Fiber Count [%];Count [Towers]");
  TH1D *h_fiber_count_dist_uiuc = new TH1D("", "", 80, 90, 105);
  TH1D *h_fiber_count_dist_china = new TH1D("", "", 80, 90, 105);
  THStack *hs_fiber_count_block_dist = new THStack("hs_block_fiber", "Distribution of EMCal Block Fiber Count;Fiber Count [%];Count [Blocks]");
  TH1D *h_fiber_count_block_dist_uiuc = new TH1D("", "", 80, 95, 101);
  TH1D *h_fiber_count_block_dist_china = new TH1D("", "", 80, 95, 101);
  gStyle->SetOptStat(0);

  double min_fiber_count = std::numeric_limits<double>().infinity();
  double min_mpv = std::numeric_limits<double>().infinity();
  double max_fiber_count = -std::numeric_limits<double>().infinity();
  double max_mpv = -std::numeric_limits<double>().infinity();
  
  std::map<std::string, int> fiber_type_counter;

  for (Block &block : all_blocks) {
    if (fiber_type_counter.find(block.fiber_type) == fiber_type_counter.end()) {
      fiber_type_counter[block.fiber_type] = 1;
    } else {
      fiber_type_counter[block.fiber_type]++;
    }
  }

  printf("FIBER TYPES:\n");
  for (auto const& p : fiber_type_counter) {
    printf("\t%s: %i\n", p.first.c_str(), p.second);
  }

  std::map<std::string, int> fiber_type_compressed;

  for (auto const& p : fiber_type_counter) {
    if (fiber_type_compressor.find(p.first) == fiber_type_compressor.end()) {
      throw std::runtime_error(Form("unknown fiber type: '%s'", p.first.c_str()));
    }
    std::string fiber_type = fiber_type_compressor.at(p.first);
    if (fiber_type_compressed.find(fiber_type) == fiber_type_compressed.end()) {
      fiber_type_compressed[fiber_type] = p.second;
    } else {
      fiber_type_compressed[fiber_type] += p.second;
    }
  }

  printf("COMPRESSED FIBER TYPES:\n");
  for (auto const& p : fiber_type_compressed) {
    printf("\t%s: %i\n", p.first.c_str(), p.second);
  }

  // std::vector<std::vector<double>> chnl_mpvs = get_chnl_mpv_with_err().first;
  TH2D *h_chnl_mpv = new TH2D("", "sPHENIX EMCal Channel MPV;#phi [Channels];#eta [Channels];MPV", 256, 0, 256, 96, 0, 96);
  for (Block &block : all_blocks) {
    auto xy = get_block_loc(block, true_sector_mapping);
    unsigned int x = 2*xy.first + 1;

    unsigned int y = 2*xy.second + 1;


    if (fiber_type_counter.find(block.fiber_type) == fiber_type_counter.end()) {
      fiber_type_counter[block.fiber_type] = 1;
    } else {
      fiber_type_counter[block.fiber_type]++;
    }

    // auto chnls = block_to_channel(block.block_number - 1);
    double ch0 = block.ch0_mpv;
    double ch1 = block.ch1_mpv;
    double ch2 = block.ch2_mpv;
    double ch3 = block.ch3_mpv;
    
    std::vector<double> ch_mpv_added;

    if (block.sector % 2 == 0) { // north sector = top half
      if (ch0 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y, ch0);
        ch_mpv_added.push_back(ch0);
      }
      if (ch1 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y + 1, ch1);
        ch_mpv_added.push_back(ch1);
      }
      if (ch2 > 0) {
        h_chnl_mpv->SetBinContent(x, y, ch2);
        ch_mpv_added.push_back(ch2);
      }
      if (ch3 > 0) {
        h_chnl_mpv->SetBinContent(x, y + 1, ch3);
        ch_mpv_added.push_back(ch3);
      }
    } else { // south sector = bottom half
      if (ch0 > 0) {
        h_chnl_mpv->SetBinContent(x, y + 1, ch0);
        ch_mpv_added.push_back(ch0);
      }
      if (ch1 > 0) {
        h_chnl_mpv->SetBinContent(x, y, ch1);
        ch_mpv_added.push_back(ch1);
      }
      if (ch2 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y + 1, ch2);
        ch_mpv_added.push_back(ch2);
      }
      if (ch3 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y, ch3);
        ch_mpv_added.push_back(ch3);
      }
    }
    for (const double& ch : ch_mpv_added) {
      h_mpv_dist->Fill(ch);
      if (block.dbn[0] != 'F' && block.dbn[0] != 'C') {
        h_mpv_dist_uiuc->Fill(ch);
      } else {
        if (fiber_type_compressor.at(block.fiber_type) == "SG") {
          h_mpv_dist_china_sg->Fill(ch);
        } else if (fiber_type_compressor.at(block.fiber_type) == "K") {
          h_mpv_dist_china_k->Fill(ch);
        }
      }
      if (ch < min_mpv) {
        min_mpv = ch;
      }
      if (ch > max_mpv) {
        max_mpv = ch;
      }
    }
    // printf("sector %2d block %2d:\n\tch0: %f\n\tch1: %f\n\tch2: %f\n\tch3: %f\n", block.sector, block.block_number, ch0, ch1, ch2, ch3);
  }


  TH2D *h_chnl_fiber = new TH2D("", "sPHENIX EMCal Tower Fiber Count;#phi [Towers];#eta [Towers];Fiber Count [%]", 256, 0, 256, 96, 0, 96);
  for (Block &block : all_blocks) {
    auto xy = get_block_loc(block, true_sector_mapping);
    unsigned int x = 2*xy.first + 1;
    unsigned int y = 2*xy.second + 1;

    double t1 = block.fiber_t1_count;
    double t2 = block.fiber_t2_count;
    double t3 = block.fiber_t3_count;
    double t4 = block.fiber_t4_count;

    std::vector<double> tower_counts_added;

    if (block.dbn[0] != 'F' && block.dbn[0] != 'C') {
      h_fiber_count_block_dist_uiuc->Fill(block.fiber_count);
    } else {
      h_fiber_count_block_dist_china->Fill(block.fiber_count);
    }

    if (block.sector % 2 == 0) { // north sector = top half
      if (t1 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y, t1);
        tower_counts_added.push_back(t1);
      }
      if (t2 > 0) {
        h_chnl_fiber->SetBinContent(x, y, t2);
        tower_counts_added.push_back(t2);
      }
      if (t3 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y + 1, t3);
        tower_counts_added.push_back(t3);
      }
      if (t4 > 0) {
        h_chnl_fiber->SetBinContent(x, y + 1, t4);
        tower_counts_added.push_back(t4);
      }
    } else { // south sector = bottom half
      if (t1 > 0) {
        h_chnl_fiber->SetBinContent(x, y + 1, t1);
        tower_counts_added.push_back(t1);
      }
      if (t2 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y + 1, t2);
        tower_counts_added.push_back(t2);
      }
      if (t3 > 0) {
        h_chnl_fiber->SetBinContent(x, y, t3);
        tower_counts_added.push_back(t3);
      }
      if (t4 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y, t4);
        tower_counts_added.push_back(t4);
      }
    }
    
    for (const double& t : tower_counts_added) {
      h_fiber_count_dist->Fill(t);
      if (block.dbn[0] != 'F' && block.dbn[0] != 'C') {
        h_fiber_count_dist_uiuc->Fill(t);
      } else {
        h_fiber_count_dist_china->Fill(t);
      }
      if (t < min_fiber_count) {
        min_fiber_count = t;
      }
      if (t > max_fiber_count) {
        max_fiber_count = t;
      }
    }

    // printf("sector %2d block %2d:\n\tch0: %f\n\tch1: %f\n\tch2: %f\n\tch3: %f\n", block.sector, block.block_number, ch0, ch1, ch2, ch3);
  }

  printf("MIN CHNL MPV = %f\n", min_mpv);
  printf("MIN CHNL FIBER COUNT = %f\n", min_fiber_count);
  printf("MAX CHNL MPV = %f\n", max_mpv);
  printf("MAX CHNL FIBER COUNT = %f\n", max_fiber_count);

  h_mpv_dist_uiuc->Fit("gaus");
  h_mpv_dist_china_sg->Fit("gaus");
  h_mpv_dist_china_k->Fit("gaus");

  h_mpv_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_mpv_dist_uiuc->SetLineWidth(1.0);
  h_mpv_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_mpv_dist_uiuc->SetFillStyle(3354);
  h_mpv_dist_uiuc->GetFunction("gaus")->SetLineColor(kBlue);
  h_mpv_dist_uiuc->GetFunction("gaus")->SetLineWidth(4.0);
  h_mpv_dist_uiuc->GetFunction("gaus")->SetLineStyle(kDashed);
  h_mpv_dist_china_sg->SetLineColorAlpha(kRed, 1);
  h_mpv_dist_china_sg->SetLineWidth(1.0);
  h_mpv_dist_china_sg->SetFillColorAlpha(kRed, 0.3);
  // h_mpv_dist_china_sg->SetFillStyle(3354);
  h_mpv_dist_china_sg->GetFunction("gaus")->SetLineColor(kRed);
  h_mpv_dist_china_sg->GetFunction("gaus")->SetLineWidth(4.0);
  h_mpv_dist_china_sg->GetFunction("gaus")->SetLineStyle(kDashed);
  h_mpv_dist_china_k->SetLineColorAlpha(kGreen, 1);
  h_mpv_dist_china_k->SetLineWidth(1.0);
  h_mpv_dist_china_k->SetFillColorAlpha(kGreen, 0.3);
  // h_mpv_dist_china_k->SetFillStyle(3354);
  h_mpv_dist_china_k->GetFunction("gaus")->SetLineColor(kGreen);
  h_mpv_dist_china_k->GetFunction("gaus")->SetLineWidth(4.0);
  h_mpv_dist_china_k->GetFunction("gaus")->SetLineStyle(kDashed);

  hs_mpv_dist->Add(h_mpv_dist_uiuc);
  hs_mpv_dist->Add(h_mpv_dist_china_sg);
  hs_mpv_dist->Add(h_mpv_dist_china_k);
  
  TLegend *hs_mpv_leg = new TLegend(.7, .7, .85, .85);
  hs_mpv_leg->AddEntry(h_mpv_dist_uiuc, "UIUC", "f");
  hs_mpv_leg->AddEntry(h_mpv_dist_china_sg, "China SG", "f");
  hs_mpv_leg->AddEntry(h_mpv_dist_china_k, "China K", "f");
  
  TCanvas *cs_mpv_dist = new TCanvas();
  hs_mpv_dist->Draw("NOSTACKB");
  hs_mpv_leg->Draw();
  cs_mpv_dist->SaveAs("emcal_plots/chnl_mpv_dist.pdf");

  // h_fiber_count_dist_uiuc->Fit("gaus");
  // h_fiber_count_dist_china->Fit("gaus");

  h_fiber_count_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_fiber_count_dist_uiuc->SetLineWidth(1.0);
  h_fiber_count_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_fiber_count_dist_uiuc->SetFillStyle(3354);
  // h_fiber_count_dist_uiuc->GetFunction("gaus")->SetLineColor(kBlue);
  // h_fiber_count_dist_uiuc->GetFunction("gaus")->SetLineWidth(4.0);
  // h_fiber_count_dist_uiuc->GetFunction("gaus")->SetLineStyle(kDashed);
  h_fiber_count_dist_china->SetLineColorAlpha(kRed, 1);
  h_fiber_count_dist_china->SetLineWidth(1.0);
  h_fiber_count_dist_china->SetFillColorAlpha(kRed, 0.3);
  // h_fiber_count_dist_china->SetFillStyle(3354);
  // h_fiber_count_dist_china->GetFunction("gaus")->SetLineColor(kRed);
  // h_fiber_count_dist_china->GetFunction("gaus")->SetLineWidth(4.0);
  // h_fiber_count_dist_china->GetFunction("gaus")->SetLineStyle(kDashed);

  // h_fiber_count_dist_uiuc->Fit("gaus");
  // h_fiber_count_dist_china->Fit("gaus");

  // h_fiber_count_block_dist_uiuc->Fit("gaus");
  // h_fiber_count_block_dist_china->Fit("gaus");

  // TF1 *skewed_guas = new TF1("sguas","2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)",95,101);

  // h_fiber_count_block_dist_uiuc->Fit(skewed_guas, "R");
  // h_fiber_count_block_dist_china->Fit(skewed_guas, "R");

  h_fiber_count_block_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_fiber_count_block_dist_uiuc->SetLineWidth(1.0);
  h_fiber_count_block_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_fiber_count_block_dist_uiuc->SetFillStyle(3354);
  // h_fiber_count_block_dist_uiuc->GetFunction("gaus")->SetLineColor(kBlue);
  // h_fiber_count_block_dist_uiuc->GetFunction("gaus")->SetLineWidth(4.0);
  // h_fiber_count_block_dist_uiuc->GetFunction("gaus")->SetLineStyle(kDashed);
  h_fiber_count_block_dist_china->SetLineColorAlpha(kRed, 1);
  h_fiber_count_block_dist_china->SetLineWidth(1.0);
  h_fiber_count_block_dist_china->SetFillColorAlpha(kRed, 0.3);
  // h_fiber_count_block_dist_china->SetFillStyle(3354);
  // h_fiber_count_block_dist_china->GetFunction("gaus")->SetLineColor(kRed);
  // h_fiber_count_block_dist_china->GetFunction("gaus")->SetLineWidth(4.0);
  // h_fiber_count_block_dist_china->GetFunction("gaus")->SetLineStyle(kDashed);

  hs_fiber_count_dist->Add(h_fiber_count_dist_uiuc);
  hs_fiber_count_dist->Add(h_fiber_count_dist_china);

  hs_fiber_count_block_dist->Add(h_fiber_count_block_dist_uiuc);
  hs_fiber_count_block_dist->Add(h_fiber_count_block_dist_china);
  
  TLegend *hs_fiber_count_leg = new TLegend(.15, .7, .3, .85);
  hs_fiber_count_leg->AddEntry(h_fiber_count_dist_uiuc, "UIUC", "f");
  hs_fiber_count_leg->AddEntry(h_fiber_count_dist_china, "China", "f");
  
  TLegend *hs_fiber_count_block_leg = new TLegend(.15, .7, .3, .85);
  hs_fiber_count_block_leg->AddEntry(h_fiber_count_block_dist_uiuc, "UIUC", "f");
  hs_fiber_count_block_leg->AddEntry(h_fiber_count_block_dist_china, "China", "f");
  
  TCanvas *cs_fiber_count_dist = new TCanvas();
  hs_fiber_count_dist->Draw("NOSTACKB");
  hs_fiber_count_leg->Draw();
  cs_fiber_count_dist->SaveAs("emcal_plots/chnl_fiber_count_dist.pdf");
  
  TCanvas *cs_fiber_count_block_dist = new TCanvas();
  hs_fiber_count_block_dist->Draw("NOSTACKB");
  hs_fiber_count_block_leg->Draw();
  cs_fiber_count_block_dist->SaveAs("emcal_plots/chnl_fiber_count_block_dist.pdf");

  TFile *chnl_dists_file = new TFile("emcal_plots/histograms.root", "RECREATE");
  chnl_dists_file->WriteObject(h_mpv_dist, "h_chnl_mpv");
  chnl_dists_file->WriteObject(h_fiber_count_dist, "h_chnl_fiber_count");
  chnl_dists_file->WriteObject(hs_mpv_dist, "hs_chnl_mpv");
  chnl_dists_file->WriteObject(hs_fiber_count_dist, "hs_chnl_fiber_count");
  chnl_dists_file->WriteObject(hs_fiber_count_block_dist, "hs_chnl_fiber_count_block");

  if (mode == "caroline") {
    TCanvas *c_chnl_mpv = new TCanvas();
    c_chnl_mpv->SetRightMargin(0.125);
    c_chnl_mpv->SetGrid();
    h_chnl_mpv->SetAxisRange(0, 128*2 - 1, "X");
    h_chnl_mpv->SetAxisRange(0, 48*2 - 1, "Y");
    h_chnl_mpv->GetXaxis()->SetNdivisions(32, false);
    h_chnl_mpv->GetYaxis()->SetNdivisions(2, false);
    h_chnl_mpv->GetXaxis()->SetLabelOffset(999.0);
    h_chnl_mpv->GetYaxis()->SetLabelOffset(999.0);
    h_chnl_mpv->GetXaxis()->SetTickLength(0);
    h_chnl_mpv->GetYaxis()->SetTickLength(0);
    h_chnl_mpv->Draw("COLZ0");
    plot_sector_and_block_labels(true);
    c_chnl_mpv->SaveAs("emcal_plots/chnl_mpv.pdf");
    
    TCanvas *c_chnl_fiber = new TCanvas();
    c_chnl_fiber->SetRightMargin(0.125);
    c_chnl_fiber->SetGrid();
    h_chnl_mpv->SetAxisRange(0, 128*2 - 1, "X");
    h_chnl_mpv->SetAxisRange(0, 48*2 - 1, "Y");
    h_chnl_fiber->GetXaxis()->SetNdivisions(32, false);
    h_chnl_fiber->GetYaxis()->SetNdivisions(2, false);
    h_chnl_fiber->GetXaxis()->SetLabelOffset(999.0);
    h_chnl_fiber->GetYaxis()->SetLabelOffset(999.0);
    h_chnl_fiber->GetXaxis()->SetTickLength(0);
    h_chnl_fiber->GetYaxis()->SetTickLength(0);
    h_chnl_fiber->Draw("COLZ0");
    plot_sector_and_block_labels(true);
    c_chnl_fiber->SaveAs("emcal_plots/chnl_fiber_count.pdf");
  } else if (mode == "tim") {
    TCanvas *c_chnl_mpv = new TCanvas();
    c_chnl_mpv->SetGrid();
    c_chnl_mpv->SetRightMargin(0.125);
    h_chnl_mpv->SetAxisRange(0, 128*2 - 1, "X");
    h_chnl_mpv->SetAxisRange(0, 48*2 - 1, "Y");
    h_chnl_mpv->GetXaxis()->SetLabelSize(0.025);
    h_chnl_mpv->GetYaxis()->SetLabelSize(0.025);
    h_chnl_mpv->Draw("COLZ0");
    h_chnl_mpv->GetXaxis()->SetNdivisions(32, false);
    h_chnl_mpv->GetYaxis()->SetNdivisions(2, false);
    
    h_chnl_mpv->GetXaxis()->SetLabelOffset(999.0);
    h_chnl_mpv->GetXaxis()->SetTickLength(0);

    // draw axes separately (since they don't match the gridlines, which are sector boundaries)
    gPad->Update();
    TGaxis *x_axis_mpv = new TGaxis(gPad->GetUxmin(),
                                  gPad->GetUymin(),
                                  gPad->GetUxmax(),
                                  gPad->GetUymin(),
                                  -4,
                                  256 - 4,
                                  16 + 4*100,"N+");            
    x_axis_mpv->SetLabelOffset(0.01);
    x_axis_mpv->SetLabelFont(42);
    x_axis_mpv->SetLabelSize(0.025);
    x_axis_mpv->Draw();

    h_chnl_mpv->GetYaxis()->SetLabelOffset(999.0);
    h_chnl_mpv->GetYaxis()->SetTickLength(0);
    gPad->Update();
    TGaxis *y_axis_mpv = new TGaxis(gPad->GetUxmin(),
                                    gPad->GetUymin(),
                                    gPad->GetUxmin(),
                                    gPad->GetUymax(),
                                    -48,
                                    48,
                                    12 + 2*100,"N-");            
    y_axis_mpv->SetLabelOffset(0.01);
    y_axis_mpv->SetLabelFont(42);
    y_axis_mpv->SetLabelSize(0.025);
    y_axis_mpv->Draw();

    plot_sector_labels_debug(true);
    c_chnl_mpv->SaveAs("emcal_plots/chnl_mpv.pdf");
    
    TCanvas *c_chnl_fiber = new TCanvas();
    c_chnl_fiber->SetGrid();
    c_chnl_fiber->SetRightMargin(0.125);
    h_chnl_fiber->SetAxisRange(0, 128*2 - 1, "X");
    h_chnl_fiber->SetAxisRange(0, 48*2 - 1, "Y");
    h_chnl_fiber->GetXaxis()->SetLabelSize(0.025);
    h_chnl_fiber->GetYaxis()->SetLabelSize(0.025);
    h_chnl_fiber->SetMinimum(90.0);
    h_chnl_fiber->SetMaximum(105.0);
    h_chnl_fiber->Draw("COLZ0");
    h_chnl_fiber->GetXaxis()->SetNdivisions(32, false);
    h_chnl_fiber->GetYaxis()->SetNdivisions(2, false);

    h_chnl_fiber->GetXaxis()->SetLabelOffset(999.0);
    h_chnl_fiber->GetXaxis()->SetTickLength(0);

    // draw axes separately (since they don't match the gridlines, which are sector boundaries)
    gPad->Update();
    TGaxis *x_axis_fiber = new TGaxis(gPad->GetUxmin(),
                                  gPad->GetUymin(),
                                  gPad->GetUxmax(),
                                  gPad->GetUymin(),
                                  -4,
                                  256 - 4,
                                  16 + 4*100,"N+");            
    x_axis_fiber->SetLabelOffset(0.01);
    x_axis_fiber->SetLabelFont(42);
    x_axis_fiber->SetLabelSize(0.025);
    x_axis_fiber->Draw();

    h_chnl_fiber->GetYaxis()->SetLabelOffset(999.0);
    h_chnl_fiber->GetYaxis()->SetTickLength(0);
    gPad->Update();
    TGaxis *y_axis_fiber = new TGaxis(gPad->GetUxmin(),
                                    gPad->GetUymin(),
                                    gPad->GetUxmin(),
                                    gPad->GetUymax(),
                                    -48,
                                    48,
                                    12 + 2*100,"N-");            
    y_axis_fiber->SetLabelOffset(0.01);
    y_axis_fiber->SetLabelFont(42);
    y_axis_fiber->SetLabelSize(0.025);
    y_axis_fiber->Draw();

    plot_sector_labels_debug(true);
    c_chnl_fiber->SaveAs("emcal_plots/chnl_fiber_count.pdf");
  } else {
    throw std::runtime_error(Form("unknown mode '%s'", mode.c_str()));
  }
}

/**
 * @brief Body of macro (called when macro is executed). 
 */
void plot() {
  check_sector_mapping(pseudo_sector_mapping);
  check_sector_mapping(true_sector_mapping);
  
  std::vector<Block> all_blocks;
  std::fstream database;
  database.open("files/sPHENIX_EMCal_blocks - dbn_mpv.csv", std::ios::in);
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
    double ch0_mpv;
    double ch0_mpv_err;
    double ch1_mpv;
    double ch1_mpv_err;
    double ch2_mpv;
    double ch2_mpv_err;
    double ch3_mpv;
    double ch3_mpv_err;
    double density;
    double fiber_count;
    double fiber_t1_count;
    double fiber_t2_count;
    double fiber_t3_count;
    double fiber_t4_count;
    double scint_ratio;
    std::string fiber_type;
    std::string fiber_batch;
    std::string w_powder;

    std::string tmp;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> sector;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> block_number;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> dbn;
    // printf("DBN WAS '%s'\n", dbn.c_str());

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
      ch0_mpv = -1;
    } else {
      std::istringstream(tmp) >> ch0_mpv;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch0_mpv_err = -1;
    } else {
      std::istringstream(tmp) >> ch0_mpv_err;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch1_mpv = -1;
    } else {
      std::istringstream(tmp) >> ch1_mpv;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch1_mpv_err = -1;
    } else {
      std::istringstream(tmp) >> ch1_mpv_err;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch2_mpv = -1;
    } else {
      std::istringstream(tmp) >> ch2_mpv;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch2_mpv_err = -1;
    } else {
      std::istringstream(tmp) >> ch2_mpv_err;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch3_mpv = -1;
    } else {
      std::istringstream(tmp) >> ch3_mpv;
    }

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      ch3_mpv_err = -1;
    } else {
      std::istringstream(tmp) >> ch3_mpv_err;
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

    std::getline(csvStream, tmp, ',');
    if (tmp == "") {
      scint_ratio = -1;
    } else {
      std::istringstream(tmp) >> scint_ratio;
    }

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> fiber_type;

    std::getline(csvStream, tmp, ',');
    std::istringstream(tmp) >> fiber_batch;

    std::getline(csvStream, tmp);
    if (tmp == "" || tmp == "\r") {
      // if (tmp == "") {
      //   printf("EMPTY STRING\n");
      // }  else {
      //   printf("CARRIAGE RETURN\n");
      // }
      w_powder = "";
    } else {
      std::istringstream(tmp) >> w_powder;
    }

    // printf("DBN %s: type '%s' batch '%s'\n", dbn.c_str(), fiber_type.c_str(), fiber_batch.c_str());

    // dbn.erase(std::remove(dbn.begin(), dbn.end(), '\n'), dbn.end());
    // dbn.erase(std::remove(dbn.begin(), dbn.end(), '\r'), dbn.end());
    // dbn.erase(std::remove(dbn.begin(), dbn.end(), ' '), dbn.end());

    all_blocks.push_back({
      sector, block_number, dbn, mpv, mpv_err,
      ch0_mpv, ch0_mpv_err, ch1_mpv, ch1_mpv_err, ch2_mpv, ch2_mpv_err, ch3_mpv, ch3_mpv_err,
      density, fiber_count, fiber_t1_count, fiber_t2_count, fiber_t3_count, fiber_t4_count, scint_ratio, fiber_type
    });
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

  plot_channel_lvl(all_blocks, "tim");
  for (const PlotConfig& cfg : cfgs) {
    plot_helper(all_blocks, cfg);
  }
}
