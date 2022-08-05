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
#include <TLatex.h>
#include <TPaveStats.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1S.h>

#include "../includes/mpv_dbn.h"
#include "../includes/utils.h"

/**
 * TODO:
 * 
 * maybe correlation plots
 * 
 * maybe do a "string hist" of fiber batch? (the bins in the correct order)
 */

/**
 * @brief list of fiber batches to treat as if they were empty (e.g., we don't have fiber batch information about the block)
 */
const std::vector<std::string> FIBER_BATCH_TREAT_EMPTY = {"0", "none"};

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

class FiberBatch {
  public:
  FiberBatch(std::string fiber_batch) : str(fiber_batch) {
    if (fiber_batch == "" || std::find(FIBER_BATCH_TREAT_EMPTY.begin(), FIBER_BATCH_TREAT_EMPTY.end(), fiber_batch) != FIBER_BATCH_TREAT_EMPTY.end()) {
      valid = false;
      batch_number = -1;
      batch_letter = 'X';
      return;
    } else {
      valid = true;
      // REALLY valid if the constructor exits without throwing a runtime_error...
    }
    size_t dash_pos = fiber_batch.find('-');
    if (dash_pos == std::string::npos) {
      throw std::runtime_error(Form("unable to locate '-' in fiber batch: %s", fiber_batch.c_str()));
    }
    std::string numeric = fiber_batch.substr(0, dash_pos);
    std::string alpha = fiber_batch.substr(dash_pos + 1, fiber_batch.length() - dash_pos);
    // check numeric is numeric
    for (const char &x: numeric) {
      if (!isdigit(x)) {
        throw std::runtime_error(Form("fiber type lhs contains a non-numeric character: %s", fiber_batch.c_str()));
      }
    }
    batch_number = std::stoi(numeric);
    // check alpha is exactly one character
    if (alpha.length() != 1) {
      throw std::runtime_error(Form("fiber type rhs contained more than one character: %s", fiber_batch.c_str()));
    }
    // check alpha is alpha and is uppercase
    for (const char &x: alpha) {
      if (!isalpha(x) || x != toupper(x)) {
        throw std::runtime_error(Form("fiber type rhs contains a non-letter or is not uppercase: %s", fiber_batch.c_str()));
      }
    }
    batch_letter = alpha[0];
  }
  bool operator ==(const FiberBatch &rhs) const {
    return batch_number == rhs.batch_number && batch_letter == rhs.batch_letter;
  }
  bool operator !=(const FiberBatch &rhs) const {
    return !(*this == rhs); 
  }
  bool operator <(const FiberBatch &rhs) const {
    if (batch_number < rhs.batch_number) {
      return true;
    } else if (batch_number == rhs.batch_number && batch_letter < rhs.batch_letter) {
      return true;
    } else {
      return false;
    }
  }
  bool operator >(const FiberBatch &rhs) const {
    return !(*this < rhs);
  }
  bool operator <=(const FiberBatch &rhs) const {
    return (*this < rhs) || (*this == rhs);
  }
  bool operator >=(const FiberBatch &rhs) const {
    return (*this > rhs) || (*this == rhs);
  }

  bool valid;
  std::string str;
  int batch_number;
  char batch_letter;
};

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
  FiberBatch fiber_batch;
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
  int color;
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
  
  TH2D* h_pseudo = new TH2D("", "", 128, 0, 128, 48, 0, 48);
  TH2D* h_true = new TH2D("", "", 128, -2, 126, 48, -24, 24);

  for (const Block& block : all_blocks) {
    auto pseudo_offsets = get_block_loc(block, pseudo_sector_mapping);
    auto true_offsets = get_block_loc(block, true_sector_mapping);
    unsigned int pseudo_x_offset = pseudo_offsets.first;
    unsigned int pseudo_y_offset = pseudo_offsets.second;
    unsigned int true_x_offset = true_offsets.first;
    unsigned int true_y_offset = true_offsets.second;
    std::pair<bool, double> p = cfg.get_value(block);
    if (p.first) {
      double value = p.second;
      h_pseudo->SetBinContent(pseudo_x_offset + 1, pseudo_y_offset + 1, value);
      h_true->SetBinContent(true_x_offset + 1, true_y_offset + 1, value);
      // printf("set bin content for %3d, %3d = sector %2d, block_num %2d\n", x_offset, y_offset, block.sector, block.block_number);
    } else {
      // printf("SKIPPED bin content for %3d, %3d\n", x_offset, y_offset);
    }
  }

  std::string title = Form("sPHENIX EMCAL %s;#phi [Blocks];#eta [Blocks];%s%s", cfg.title.c_str(), cfg.title.c_str(), cfg.units.c_str());

  TCanvas* c0 = new TCanvas("c0", "", 700, 500);
  c0->SetRightMargin(0.125);
  c0->SetGrid();
  gStyle->SetOptStat(0);
  h_pseudo->SetTitle(title.c_str());
  h_pseudo->SetMinimum(cfg.plot_min);
  h_pseudo->SetAxisRange(0, 128 - 1, "X");
  h_pseudo->SetAxisRange(0, 48 - 1, "Y");
  // h_pseudo->GetXaxis()->SetRange(0, 128);
  // h_pseudo->GetYaxis()->SetRange(0, 48);
  h_pseudo->GetXaxis()->SetNdivisions(32, false);
  h_pseudo->GetYaxis()->SetNdivisions(2, false);
  h_pseudo->GetXaxis()->SetLabelOffset(999.0);
  h_pseudo->GetYaxis()->SetLabelOffset(999.0);
  h_pseudo->GetXaxis()->SetTickLength(0);
  h_pseudo->GetYaxis()->SetTickLength(0);
  h_pseudo->Draw("COLZ0");
  plot_sector_and_block_labels();
  c0->SaveAs(Form("emcal_plots/plot_colz_%s.pdf", cfg.file_name.c_str()));

  TCanvas* c1 = new TCanvas("c1", "", 1200, 500);
  c1->SetRightMargin(0.125);
  c1->SetGrid();
  gStyle->SetOptStat(0);
  h_pseudo->SetTitle(title.c_str());
  h_pseudo->SetMinimum(cfg.plot_min);
  h_pseudo->SetAxisRange(0, 128 - 1, "X");
  h_pseudo->SetAxisRange(0, 48 - 1, "Y");
  // h_pseudo->GetXaxis()->SetRange(0, 128);
  // h_pseudo->GetYaxis()->SetRange(0, 48);
  h_pseudo->GetXaxis()->SetNdivisions(32, false);
  h_pseudo->GetYaxis()->SetNdivisions(2, false);
  h_pseudo->GetXaxis()->SetLabelOffset(999.0);
  h_pseudo->GetYaxis()->SetLabelOffset(999.0);
  h_pseudo->GetXaxis()->SetTickLength(0);
  h_pseudo->GetYaxis()->SetTickLength(0);
  h_pseudo->Draw("COLZ0");

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
  
  h_true->GetXaxis()->SetTitleOffset(1.7);
  h_true->GetYaxis()->SetTitleOffset(1.7);
  h_true->GetZaxis()->SetTitleOffset(1.2);

  TCanvas* c2 = new TCanvas("c2", "", 900, 500);
  h_true->SetLineColorAlpha(kBlack, 1.0);
  h_true->SetTitle(title.c_str());
  h_true->Draw("LEGO2 0");
  c2->SaveAs(Form("emcal_plots/plot_lego_%s.png", cfg.file_name.c_str()));

  TCanvas* c3 = new TCanvas("c3", "", 900, 500);
  h_true->SetLineColorAlpha(kBlack, 1.0);
  h_true->SetFillColorAlpha(cfg.color, 1.0);
  h_true->SetTitle(title.c_str());
  h_true->Draw("LEGO1 CYL 0");
  c3->SaveAs(Form("emcal_plots/plot_cyl_%s_front.png", cfg.file_name.c_str()));
  gPad->SetPhi(gPad->GetPhi() + 180);
  gPad->Update();
  c3->SaveAs(Form("emcal_plots/plot_cyl_%s_back.png", cfg.file_name.c_str()));
  
  TCanvas* c4 = new TCanvas("c4", "", 900, 500);
  h_true->SetLineColorAlpha(kBlack, 1.0);
  h_true->SetFillColorAlpha(cfg.color, 1.0);
  h_true->SetTitle(title.c_str());
  h_true->Draw("LEGO1 PSR 0");
  c4->SaveAs(Form("emcal_plots/plot_psr_%s_front.png", cfg.file_name.c_str()));
  gPad->SetPhi(gPad->GetPhi() + 180);
  gPad->Update();
  c4->SaveAs(Form("emcal_plots/plot_psr_%s_back.png", cfg.file_name.c_str()));
  
  delete c0;
  delete c1;
  delete c2;
  delete c3;
  delete c4;

  gStyle->SetImageScaling(1.0);
}

/**
 * @brief 
 * 
 * @param all_blocks 
 */
void make_histograms(std::vector<Block> all_blocks) {
  TH1D *h_mpv_block_dist = new TH1D("h_mpv_block_dist", "Distribution of EMCal Block MPV;MPV;Count [Blocks]", 80, 0, 1000);
  TH1D *h_mpv_chnl_dist = new TH1D("h_mpv_chnl_dist", "Distribution of EMCal Channel MPV;MPV;Count [Channels]", 80, 0, 1000);
  
  TH1D *h_fiber_count_block_dist = new TH1D("h_fiber_count_block_dist", "Distribution of EMCal Block Fiber Count;Fiber Count [%];Count [Block]", 80, 80, 120);
  TH1D *h_fiber_count_tower_dist = new TH1D("h_fiber_count_tower_dist", "Distribution of EMCal Tower Fiber Count;Fiber Count [%];Count [Towers]", 80, 80, 120);
  
  THStack *hs_mpv_block_dist = new THStack("hs_mpv_block_dist", "Distribution of EMCal Block MPV;MPV;Count [Blocks]");
  TH1D *h_mpv_block_dist_uiuc_1 = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_block_dist_uiuc_2 = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_block_dist_china_sg = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_block_dist_china_k = new TH1D("", "", 80, 0, 800);
  
  THStack *hs_mpv_chnl_dist = new THStack("hs_mpv_chnl_dist", "Distribution of EMCal Channel MPV;MPV;Count [Channels]");
  TH1D *h_mpv_chnl_dist_uiuc_1 = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_chnl_dist_uiuc_2 = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_chnl_dist_china_sg = new TH1D("", "", 80, 0, 800);
  TH1D *h_mpv_chnl_dist_china_k = new TH1D("", "", 80, 0, 800);

  THStack *hs_fiber_count_block_dist = new THStack("hs_fiber_count_block_dist", "Distribution of EMCal Block Fiber Count;Fiber Count [%];Count [Blocks]");
  TH1D *h_fiber_count_block_dist_uiuc = new TH1D("", "", 80, 95, 101);
  TH1D *h_fiber_count_block_dist_fudan = new TH1D("", "", 80, 95, 101);
  TH1D *h_fiber_count_block_dist_ciae = new TH1D("", "", 80, 95, 101);
  
  THStack *hs_fiber_count_tower_dist = new THStack("hs_fiber_count_tower_dist", "Distribution of EMCal Tower Fiber Count;Fiber Count [%];Count [Towers]");
  TH1D *h_fiber_count_tower_dist_uiuc = new TH1D("", "", 80, 90, 105);
  TH1D *h_fiber_count_tower_dist_fudan = new TH1D("", "", 80, 90, 105);
  TH1D *h_fiber_count_tower_dist_ciae = new TH1D("", "", 80, 90, 105);
  
  THStack *hs_density_dist = new THStack("hs_density_dist", "Distribution of EMCal Block Density;Density [g/mL];Count [Blocks]");
  TH1D *h_density_dist_uiuc = new TH1D("", "", 80, 8, 11);
  TH1D *h_density_dist_fudan = new TH1D("", "", 80, 8, 11);
  TH1D *h_density_dist_ciae = new TH1D("", "", 80, 8, 11);
  
  THStack *hs_scint_ratio_dist = new THStack("hs_scint_ratio_dist", "Distribution of EMCal Block Scintillation Ratio;Scintillation Ratio;Count [Blocks]");
  TH1D *h_scint_ratio_dist_uiuc = new TH1D("", "", 80, 0, 4);
  TH1D *h_scint_ratio_dist_fudan = new TH1D("", "", 80, 0, 4);
  TH1D *h_scint_ratio_dist_ciae = new TH1D("", "", 80, 0, 4);


  std::map<std::string, int> fiber_type_counter;
  std::map<std::string, int> fiber_batch_counter;
  std::set<std::string> fiber_batch_set;
  std::map<int, int> fiber_batch_number_counter;
  std::set<int> fiber_batch_number_set;

  for (Block &block : all_blocks) {
    if (fiber_type_counter.find(block.fiber_type) == fiber_type_counter.end()) {
      fiber_type_counter[block.fiber_type] = 1;
    } else {
      fiber_type_counter[block.fiber_type]++;
    }
    if (block.fiber_batch.valid) {
      fiber_batch_set.insert(block.fiber_batch.str);
      if (fiber_batch_counter.find(block.fiber_batch.str) == fiber_batch_counter.end()) {
        fiber_batch_counter[block.fiber_batch.str] = 1;
      } else {
        fiber_batch_counter[block.fiber_batch.str]++;
      }
      fiber_batch_number_set.insert(block.fiber_batch.batch_number);
      if (fiber_batch_number_counter.find(block.fiber_batch.batch_number) == fiber_batch_number_counter.end()) {
        fiber_batch_number_counter[block.fiber_batch.batch_number] = 1;
      } else {
        fiber_batch_number_counter[block.fiber_batch.batch_number]++;
      }
    }
  }

  std::vector<std::string> fiber_batch_list;
  for (const std::string &x : fiber_batch_set) {
    fiber_batch_list.push_back(x);
  }
  std::sort(fiber_batch_list.begin(), fiber_batch_list.end());
  TH1S *h_fiber_batches = new TH1S("h_fiber_batches", "Distribution of EMCal Fiber Batch; Fiber Batch; Count [Blocks]", fiber_batch_set.size(), 0, fiber_batch_set.size());
  for (size_t i = 0; i < fiber_batch_set.size(); i++) {
    h_fiber_batches->SetBinContent(i + 1, fiber_batch_counter.at(fiber_batch_list.at(i)));
  }
  std::cout << "THERE ARE " << fiber_batch_set.size() << " UNIQUE FIBER BATCHES" << std::endl;
  
  std::vector<int> fiber_batch_number_list;
  for (const int &x : fiber_batch_number_set) {
    fiber_batch_number_list.push_back(x);
  }
  std::sort(fiber_batch_number_list.begin(), fiber_batch_number_list.end());
  TH1S *h_fiber_batch_numbers = new TH1S("h_fiber_batch_numbers", "Distribution of EMCal Fiber Batch Number; Fiber Batch Number; Count [Blocks]", fiber_batch_number_set.size(), 0, fiber_batch_number_set.size());
  for (size_t i = 0; i < fiber_batch_number_set.size(); i++) {
    h_fiber_batch_numbers->SetBinContent(i + 1, fiber_batch_number_counter.at(fiber_batch_number_list.at(i)));
  }
  std::cout << "THERE ARE " << fiber_batch_number_set.size() << " UNIQUE FIBER BATCHES" << std::endl;

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

  for (Block &block : all_blocks) {
    auto xy = get_block_loc(block, true_sector_mapping);
    unsigned int x = 2*xy.first + 1;

    unsigned int y = 2*xy.second + 1;


    if (fiber_type_counter.find(block.fiber_type) == fiber_type_counter.end()) {
      fiber_type_counter[block.fiber_type] = 1;
    } else {
      fiber_type_counter[block.fiber_type]++;
    }
  }

  for (const Block &block : all_blocks) {
    std::vector<double> chnl_mpvs = {block.ch0_mpv, block.ch1_mpv, block.ch2_mpv, block.ch3_mpv};
    std::vector<double> tower_counts = {block.fiber_t1_count, block.fiber_t2_count, block.fiber_t3_count, block.fiber_t4_count};
    
    if (block.mpv > 0) {
      h_mpv_block_dist->Fill(block.mpv);
    }
    for (const double &chnl : chnl_mpvs) {
      if (chnl > 0) {
        h_mpv_chnl_dist->Fill(chnl);
      }
    }
    
    if (block.fiber_count > 0) {
      h_fiber_count_block_dist->Fill(block.fiber_count);
    }
    for (const double &tower : tower_counts) {
      if (tower > 0) {
        h_fiber_count_tower_dist->Fill(tower);
      }
    }

    if (block.dbn[0] != 'F' && block.dbn[0] != 'C') {
      // UIUC
      if (block.fiber_batch.valid) {
        if (block.fiber_batch <= FiberBatch("16-B")) {
          // "S1-12" 
          if (block.mpv > 0) {
            h_mpv_block_dist_uiuc_1->Fill(block.mpv);
          }
          for (const double &chnl : chnl_mpvs) {
            if (chnl > 0) {
              h_mpv_chnl_dist_uiuc_1->Fill(chnl);
            }
          }
        } else if (block.fiber_batch >= FiberBatch("21-A")) {
          // "S13-64"
          if (block.mpv > 0) {
            h_mpv_block_dist_uiuc_2->Fill(block.mpv);
          }
          for (const double &chnl : chnl_mpvs) {
            if (chnl > 0) {
              h_mpv_chnl_dist_uiuc_2->Fill(chnl);
            }
          }
        } else {
          // panic?
          std::cout << "PANIC: Fiber Batch " << block.fiber_batch.str << std::endl;
        }
      }
      // NOT UIUC
    } else if (fiber_type_compressor.at(block.fiber_type) == "SG") {
      if (block.mpv > 0) {
        h_mpv_block_dist_china_sg->Fill(block.mpv);
      }
      for (const double &chnl : chnl_mpvs) {
        if (chnl > 0) {
          h_mpv_chnl_dist_china_sg->Fill(chnl);
        }
      }
    } else if (fiber_type_compressor.at(block.fiber_type) == "K") {
      if (block.mpv > 0) {
        h_mpv_block_dist_china_k->Fill(block.mpv);
      }
      for (const double &chnl : chnl_mpvs) {
        if (chnl > 0) {
          h_mpv_chnl_dist_china_k->Fill(chnl);
        }
      }
    } else {
      // panic?
      std::cout << "PANIC: Fiber Type " << block.fiber_type << std::endl;
    }
    
    if (block.dbn[0] != 'F' && block.dbn[0] != 'C') {
      // UIUC
      if (block.fiber_count > 0) {
        h_fiber_count_block_dist_uiuc->Fill(block.fiber_count);
      }
      for (const double &tower : tower_counts) {
        if (tower > 0) {
          h_fiber_count_tower_dist_uiuc->Fill(tower);
        }
      }
      if (block.density > 0) {
        h_density_dist_uiuc->Fill(block.density);
      }
      if (block.scint_ratio > 0) {
        h_scint_ratio_dist_uiuc->Fill(block.scint_ratio);
      }
    } else if (block.dbn[0] == 'F') {
      // FUDAN
      if (block.fiber_count > 0) {
        h_fiber_count_block_dist_fudan->Fill(block.fiber_count);
      }
      for (const double &tower : tower_counts) {
        if (tower > 0) {
          h_fiber_count_tower_dist_fudan->Fill(tower);
        }
      }
      if (block.density > 0) {
        h_density_dist_fudan->Fill(block.density);
      }
      if (block.scint_ratio > 0) {
        h_scint_ratio_dist_fudan->Fill(block.scint_ratio);
      }
    } else if (block.dbn[0] == 'C') {
      // CIAE
      if (block.fiber_count > 0) {
        h_fiber_count_block_dist_ciae->Fill(block.fiber_count);
      }
      for (const double &tower : tower_counts) {
        if (tower > 0) {
          h_fiber_count_tower_dist_ciae->Fill(tower);
        }
      }
      if (block.density > 0) {
        h_density_dist_ciae->Fill(block.density);
      }
      if (block.scint_ratio > 0) {
        h_scint_ratio_dist_ciae->Fill(block.scint_ratio);
      }
    } else {
      // panic?
      std::cout << "PANIC: DBN " << block.dbn << std::endl;
    }
  }

  h_mpv_block_dist_uiuc_1->Fit("gaus");
  h_mpv_block_dist_uiuc_2->Fit("gaus");
  h_mpv_block_dist_china_sg->Fit("gaus");
  h_mpv_block_dist_china_k->Fit("gaus");
  
  h_mpv_chnl_dist_uiuc_1->Fit("gaus");
  h_mpv_chnl_dist_uiuc_2->Fit("gaus");
  h_mpv_chnl_dist_china_sg->Fit("gaus");
  h_mpv_chnl_dist_china_k->Fit("gaus");

  h_mpv_block_dist_uiuc_1->SetLineColorAlpha(kBlue, 1);
  h_mpv_block_dist_uiuc_1->SetLineWidth(1.0);
  h_mpv_block_dist_uiuc_1->SetFillColorAlpha(kBlue, 0.3);
  // h_mpv_block_dist_uiuc->SetFillStyle(3354);
  h_mpv_block_dist_uiuc_1->GetFunction("gaus")->SetLineColor(kBlue);
  h_mpv_block_dist_uiuc_1->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_block_dist_uiuc_1->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_block_dist_uiuc_2->SetLineColorAlpha(kOrange + 7, 1);
  h_mpv_block_dist_uiuc_2->SetLineWidth(1.0);
  h_mpv_block_dist_uiuc_2->SetFillColorAlpha(kOrange + 7, 0.3);
  // h_mpv_block_dist_uiuc->SetFillStyle(3354);
  h_mpv_block_dist_uiuc_2->GetFunction("gaus")->SetLineColor(kOrange + 7);
  h_mpv_block_dist_uiuc_2->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_block_dist_uiuc_2->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_block_dist_china_sg->SetLineColorAlpha(kRed, 1);
  h_mpv_block_dist_china_sg->SetLineWidth(1.0);
  h_mpv_block_dist_china_sg->SetFillColorAlpha(kRed, 0.3);
  // h_mpv_block_dist_china_sg->SetFillStyle(3354);
  h_mpv_block_dist_china_sg->GetFunction("gaus")->SetLineColor(kRed);
  h_mpv_block_dist_china_sg->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_block_dist_china_sg->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_block_dist_china_k->SetLineColorAlpha(kGreen, 1);
  h_mpv_block_dist_china_k->SetLineWidth(1.0);
  h_mpv_block_dist_china_k->SetFillColorAlpha(kGreen, 0.3);
  // h_mpv_block_dist_china_k->SetFillStyle(3354);
  h_mpv_block_dist_china_k->GetFunction("gaus")->SetLineColor(kGreen);
  h_mpv_block_dist_china_k->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_block_dist_china_k->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_chnl_dist_uiuc_1->SetLineColorAlpha(kBlue, 1);
  h_mpv_chnl_dist_uiuc_1->SetLineWidth(1.0);
  h_mpv_chnl_dist_uiuc_1->SetFillColorAlpha(kBlue, 0.3);
  // h_mpv_chnl_dist_uiuc->SetFillStyle(3354);
  h_mpv_chnl_dist_uiuc_1->GetFunction("gaus")->SetLineColor(kBlue);
  h_mpv_chnl_dist_uiuc_1->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_chnl_dist_uiuc_1->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_chnl_dist_uiuc_2->SetLineColorAlpha(kOrange + 7, 1);
  h_mpv_chnl_dist_uiuc_2->SetLineWidth(1.0);
  h_mpv_chnl_dist_uiuc_2->SetFillColorAlpha(kOrange + 7, 0.3);
  // h_mpv_chnl_dist_uiuc->SetFillStyle(3354);
  h_mpv_chnl_dist_uiuc_2->GetFunction("gaus")->SetLineColor(kOrange + 7);
  h_mpv_chnl_dist_uiuc_2->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_chnl_dist_uiuc_2->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_chnl_dist_china_sg->SetLineColorAlpha(kRed, 1);
  h_mpv_chnl_dist_china_sg->SetLineWidth(1.0);
  h_mpv_chnl_dist_china_sg->SetFillColorAlpha(kRed, 0.3);
  // h_mpv_chnl_dist_china_sg->SetFillStyle(3354);
  h_mpv_chnl_dist_china_sg->GetFunction("gaus")->SetLineColor(kRed);
  h_mpv_chnl_dist_china_sg->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_chnl_dist_china_sg->GetFunction("gaus")->SetLineStyle(kDashed);
  
  h_mpv_chnl_dist_china_k->SetLineColorAlpha(kGreen, 1);
  h_mpv_chnl_dist_china_k->SetLineWidth(1.0);
  h_mpv_chnl_dist_china_k->SetFillColorAlpha(kGreen, 0.3);
  // h_mpv_chnl_dist_china_k->SetFillStyle(3354);
  h_mpv_chnl_dist_china_k->GetFunction("gaus")->SetLineColor(kGreen);
  h_mpv_chnl_dist_china_k->GetFunction("gaus")->SetLineWidth(2.0);
  h_mpv_chnl_dist_china_k->GetFunction("gaus")->SetLineStyle(kDashed);

  hs_mpv_block_dist->Add(h_mpv_block_dist_uiuc_1);
  hs_mpv_block_dist->Add(h_mpv_block_dist_uiuc_2);
  hs_mpv_block_dist->Add(h_mpv_block_dist_china_sg);
  hs_mpv_block_dist->Add(h_mpv_block_dist_china_k);
  
  hs_mpv_chnl_dist->Add(h_mpv_chnl_dist_uiuc_1);
  hs_mpv_chnl_dist->Add(h_mpv_chnl_dist_uiuc_2);
  hs_mpv_chnl_dist->Add(h_mpv_chnl_dist_china_sg);
  hs_mpv_chnl_dist->Add(h_mpv_chnl_dist_china_k);
  
  TLegend *hs_mpv_leg = new TLegend(.7, .7, .85, .85);
  hs_mpv_leg->AddEntry(h_mpv_block_dist_uiuc_1, "UIUC S1-12", "f");
  hs_mpv_leg->AddEntry(h_mpv_block_dist_uiuc_2, "UIUC S13-64", "f");
  hs_mpv_leg->AddEntry(h_mpv_block_dist_china_sg, "China SG", "f");
  hs_mpv_leg->AddEntry(h_mpv_block_dist_china_k, "China K", "f");
  
  TCanvas *cs_mpv_block_dist = new TCanvas();
  hs_mpv_block_dist->Draw("NOSTACKB");
  hs_mpv_leg->Draw();
  TPaveStats *pt_mpv_block_dist = new TPaveStats(0.65, 0.2, 0.875, 0.65, "NDC");
  pt_mpv_block_dist->SetTextFont(42);
  pt_mpv_block_dist->AddText("EMCal Block MPV")->SetTextFont(62);
  pt_mpv_block_dist->AddText("UIUC S1-12")->SetTextFont(62);
  pt_mpv_block_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_block_dist_uiuc_1->GetFunction("gaus")->GetParameter(1), h_mpv_block_dist_uiuc_1->GetFunction("gaus")->GetParError(1)));
  pt_mpv_block_dist->AddText(Form("Mean = %7.3f", h_mpv_block_dist_uiuc_1->GetMean()));
  pt_mpv_block_dist->AddText("UIUC S13-64")->SetTextFont(62);
  pt_mpv_block_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_block_dist_uiuc_2->GetFunction("gaus")->GetParameter(1), h_mpv_block_dist_uiuc_2->GetFunction("gaus")->GetParError(1)));
  pt_mpv_block_dist->AddText(Form("Mean = %7.3f", h_mpv_block_dist_uiuc_2->GetMean()));
  pt_mpv_block_dist->AddText("China SG")->SetTextFont(62);
  pt_mpv_block_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_block_dist_china_sg->GetFunction("gaus")->GetParameter(1), h_mpv_block_dist_china_sg->GetFunction("gaus")->GetParError(1)));
  pt_mpv_block_dist->AddText(Form("Mean = %7.3f", h_mpv_block_dist_china_sg->GetMean()));
  pt_mpv_block_dist->AddText("China K")->SetTextFont(62);
  pt_mpv_block_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_block_dist_china_k->GetFunction("gaus")->GetParameter(1), h_mpv_block_dist_china_k->GetFunction("gaus")->GetParError(1)));
  pt_mpv_block_dist->AddText(Form("Mean = %7.3f", h_mpv_block_dist_china_k->GetMean()));
  pt_mpv_block_dist->SetFillColorAlpha(0, 0);
  pt_mpv_block_dist->SetTextAlign(12);
  pt_mpv_block_dist->SetTextSize(0.02);
  pt_mpv_block_dist->SetShadowColor(kWhite);
  pt_mpv_block_dist->Draw();
  cs_mpv_block_dist->SaveAs("emcal_plots/mpv_block_dist.pdf");

  TCanvas *cs_mpv_chnl_dist = new TCanvas();
  hs_mpv_chnl_dist->Draw("NOSTACKB");
  hs_mpv_leg->Draw();
  TPaveStats *pt_mpv_chnl_dist = new TPaveStats(0.65, 0.2, 0.875, 0.65, "NDC");
  pt_mpv_chnl_dist->SetTextFont(42);
  pt_mpv_chnl_dist->AddText("EMCal Channel MPV")->SetTextFont(62);
  pt_mpv_chnl_dist->AddText("UIUC S1-12")->SetTextFont(62);
  pt_mpv_chnl_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_chnl_dist_uiuc_1->GetFunction("gaus")->GetParameter(1), h_mpv_chnl_dist_uiuc_1->GetFunction("gaus")->GetParError(1)));
  pt_mpv_chnl_dist->AddText(Form("Mean = %7.3f", h_mpv_chnl_dist_uiuc_1->GetMean()));
  pt_mpv_chnl_dist->AddText("UIUC S13-64")->SetTextFont(62);
  pt_mpv_chnl_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_chnl_dist_uiuc_2->GetFunction("gaus")->GetParameter(1), h_mpv_chnl_dist_uiuc_2->GetFunction("gaus")->GetParError(1)));
  pt_mpv_chnl_dist->AddText(Form("Mean = %7.3f", h_mpv_chnl_dist_uiuc_2->GetMean()));
  pt_mpv_chnl_dist->AddText("China SG")->SetTextFont(62);
  pt_mpv_chnl_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_chnl_dist_china_sg->GetFunction("gaus")->GetParameter(1), h_mpv_chnl_dist_china_sg->GetFunction("gaus")->GetParError(1)));
  pt_mpv_chnl_dist->AddText(Form("Mean = %7.3f", h_mpv_chnl_dist_china_sg->GetMean()));
  pt_mpv_chnl_dist->AddText("China K")->SetTextFont(62);
  pt_mpv_chnl_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_mpv_chnl_dist_china_k->GetFunction("gaus")->GetParameter(1), h_mpv_chnl_dist_china_k->GetFunction("gaus")->GetParError(1)));
  pt_mpv_chnl_dist->AddText(Form("Mean = %7.3f", h_mpv_chnl_dist_china_k->GetMean()));
  pt_mpv_chnl_dist->SetFillColorAlpha(0, 0);
  pt_mpv_chnl_dist->SetTextAlign(12);
  pt_mpv_chnl_dist->SetTextSize(0.02);
  pt_mpv_chnl_dist->SetShadowColor(kWhite);
  pt_mpv_chnl_dist->Draw();
  cs_mpv_block_dist->SaveAs("emcal_plots/mpv_block_dist.pdf");
  cs_mpv_chnl_dist->SaveAs("emcal_plots/mpv_chnl_dist.pdf");

  h_fiber_count_block_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_fiber_count_block_dist_uiuc->SetLineWidth(1.0);
  h_fiber_count_block_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_fiber_count_block_dist_uiuc->SetFillStyle(3354);

  h_fiber_count_block_dist_fudan->SetLineColorAlpha(kRed, 1);
  h_fiber_count_block_dist_fudan->SetLineWidth(1.0);
  h_fiber_count_block_dist_fudan->SetFillColorAlpha(kRed, 0.3);
  // h_fiber_count_block_dist_fudan->SetFillStyle(3354);

  h_fiber_count_block_dist_ciae->SetLineColorAlpha(kGreen, 1);
  h_fiber_count_block_dist_ciae->SetLineWidth(1.0);
  h_fiber_count_block_dist_ciae->SetFillColorAlpha(kGreen, 0.3);
  // h_fiber_count_block_dist_ciae->SetFillStyle(3354);

  h_fiber_count_tower_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_fiber_count_tower_dist_uiuc->SetLineWidth(1.0);
  h_fiber_count_tower_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_fiber_count_tower_dist_uiuc->SetFillStyle(3354);

  h_fiber_count_tower_dist_fudan->SetLineColorAlpha(kRed, 1);
  h_fiber_count_tower_dist_fudan->SetLineWidth(1.0);
  h_fiber_count_tower_dist_fudan->SetFillColorAlpha(kRed, 0.3);
  // h_fiber_count_tower_dist_fudan->SetFillStyle(3354);

  h_fiber_count_tower_dist_ciae->SetLineColorAlpha(kGreen, 1);
  h_fiber_count_tower_dist_ciae->SetLineWidth(1.0);
  h_fiber_count_tower_dist_ciae->SetFillColorAlpha(kGreen, 0.3);
  // h_fiber_count_tower_dist_ciae->SetFillStyle(3354);

  h_density_dist_uiuc->Fit("gaus");
  h_density_dist_fudan->Fit("gaus");
  h_density_dist_ciae->Fit("gaus");
  h_scint_ratio_dist_uiuc->Fit("gaus");
  h_scint_ratio_dist_fudan->Fit("gaus");
  h_scint_ratio_dist_ciae->Fit("gaus");

  h_density_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_density_dist_uiuc->SetLineWidth(1.0);
  h_density_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_density_dist_uiuc->SetFillStyle(3354);
  h_density_dist_uiuc->GetFunction("gaus")->SetLineColor(kBlue);
  h_density_dist_uiuc->GetFunction("gaus")->SetLineWidth(2.0);
  h_density_dist_uiuc->GetFunction("gaus")->SetLineStyle(kDashed);

  h_density_dist_fudan->SetLineColorAlpha(kRed, 1);
  h_density_dist_fudan->SetLineWidth(1.0);
  h_density_dist_fudan->SetFillColorAlpha(kRed, 0.3);
  // h_density_dist_fudan->SetFillStyle(3354);
  h_density_dist_fudan->GetFunction("gaus")->SetLineColor(kRed);
  h_density_dist_fudan->GetFunction("gaus")->SetLineWidth(2.0);
  h_density_dist_fudan->GetFunction("gaus")->SetLineStyle(kDashed);

  h_density_dist_ciae->SetLineColorAlpha(kGreen, 1);
  h_density_dist_ciae->SetLineWidth(1.0);
  h_density_dist_ciae->SetFillColorAlpha(kGreen, 0.3);
  // h_density_dist_ciae->SetFillStyle(3354);
  h_density_dist_ciae->GetFunction("gaus")->SetLineColor(kGreen);
  h_density_dist_ciae->GetFunction("gaus")->SetLineWidth(2.0);
  h_density_dist_ciae->GetFunction("gaus")->SetLineStyle(kDashed);

  h_scint_ratio_dist_uiuc->SetLineColorAlpha(kBlue, 1);
  h_scint_ratio_dist_uiuc->SetLineWidth(1.0);
  h_scint_ratio_dist_uiuc->SetFillColorAlpha(kBlue, 0.3);
  // h_scint_ratio_dist_uiuc->SetFillStyle(3354);
  h_scint_ratio_dist_uiuc->GetFunction("gaus")->SetLineColor(kBlue);
  h_scint_ratio_dist_uiuc->GetFunction("gaus")->SetLineWidth(2.0);
  h_scint_ratio_dist_uiuc->GetFunction("gaus")->SetLineStyle(kDashed);

  h_scint_ratio_dist_fudan->SetLineColorAlpha(kRed, 1);
  h_scint_ratio_dist_fudan->SetLineWidth(1.0);
  h_scint_ratio_dist_fudan->SetFillColorAlpha(kRed, 0.3);
  // h_scint_ratio_dist_fudan->SetFillStyle(3354);
  h_scint_ratio_dist_fudan->GetFunction("gaus")->SetLineColor(kRed);
  h_scint_ratio_dist_fudan->GetFunction("gaus")->SetLineWidth(2.0);
  h_scint_ratio_dist_fudan->GetFunction("gaus")->SetLineStyle(kDashed);

  h_scint_ratio_dist_ciae->SetLineColorAlpha(kGreen, 1);
  h_scint_ratio_dist_ciae->SetLineWidth(1.0);
  h_scint_ratio_dist_ciae->SetFillColorAlpha(kGreen, 0.3);
  // h_scint_ratio_dist_ciae->SetFillStyle(3354);
  h_scint_ratio_dist_ciae->GetFunction("gaus")->SetLineColor(kGreen);
  h_scint_ratio_dist_ciae->GetFunction("gaus")->SetLineWidth(2.0);
  h_scint_ratio_dist_ciae->GetFunction("gaus")->SetLineStyle(kDashed);

  hs_fiber_count_block_dist->Add(h_fiber_count_block_dist_uiuc);
  hs_fiber_count_block_dist->Add(h_fiber_count_block_dist_fudan);
  hs_fiber_count_block_dist->Add(h_fiber_count_block_dist_ciae);

  hs_fiber_count_tower_dist->Add(h_fiber_count_tower_dist_uiuc);
  hs_fiber_count_tower_dist->Add(h_fiber_count_tower_dist_fudan);
  hs_fiber_count_tower_dist->Add(h_fiber_count_tower_dist_ciae);
  
  hs_density_dist->Add(h_density_dist_uiuc);
  hs_density_dist->Add(h_density_dist_fudan);
  hs_density_dist->Add(h_density_dist_ciae);
  
  hs_scint_ratio_dist->Add(h_scint_ratio_dist_uiuc);
  hs_scint_ratio_dist->Add(h_scint_ratio_dist_fudan);
  hs_scint_ratio_dist->Add(h_scint_ratio_dist_ciae);
  
  TLegend *hs_fiber_count_leg = new TLegend(.15, .7, .3, .85);
  hs_fiber_count_leg->AddEntry(h_fiber_count_block_dist_uiuc, "UIUC", "f");
  hs_fiber_count_leg->AddEntry(h_fiber_count_block_dist_fudan, "Fudan", "f");
  hs_fiber_count_leg->AddEntry(h_fiber_count_block_dist_ciae, "CIAE", "f");
  
  TCanvas *cs_fiber_count_block_dist = new TCanvas();
  hs_fiber_count_block_dist->Draw("NOSTACKB");
  hs_fiber_count_leg->Draw();
  TPaveStats *pt_fiber_count_block_dist = new TPaveStats(0.15, 0.45, 0.35, 0.65, "NDC");
  pt_fiber_count_block_dist->SetTextFont(42);
  pt_fiber_count_block_dist->AddText("EMCal Block Fiber Count")->SetTextFont(62);
  pt_fiber_count_block_dist->AddText("UIUC")->SetTextFont(62);
  pt_fiber_count_block_dist->AddText(Form("Mean = %7.3f", h_fiber_count_block_dist_uiuc->GetMean()));
  pt_fiber_count_block_dist->AddText("Fudan")->SetTextFont(62);
  pt_fiber_count_block_dist->AddText(Form("Mean = %7.3f", h_fiber_count_block_dist_fudan->GetMean()));
  pt_fiber_count_block_dist->AddText("CIAE")->SetTextFont(62);
  pt_fiber_count_block_dist->AddText(Form("Mean = %7.3f", h_fiber_count_block_dist_ciae->GetMean()));
  pt_fiber_count_block_dist->SetFillColorAlpha(0, 0);
  pt_fiber_count_block_dist->SetTextAlign(12);
  pt_fiber_count_block_dist->SetTextSize(0.02);
  pt_fiber_count_block_dist->SetShadowColor(kWhite);
  pt_fiber_count_block_dist->Draw();
  cs_fiber_count_block_dist->SaveAs("emcal_plots/fiber_count_block_dist.pdf");
  
  TCanvas *cs_fiber_count_tower_dist = new TCanvas();
  hs_fiber_count_tower_dist->Draw("NOSTACKB");
  hs_fiber_count_leg->Draw();
  TPaveStats *pt_fiber_count_tower_dist = new TPaveStats(0.15, 0.45, 0.35, 0.65, "NDC");
  pt_fiber_count_tower_dist->SetTextFont(42);
  pt_fiber_count_tower_dist->AddText("EMCal Tower Fiber Count")->SetTextFont(62);
  pt_fiber_count_tower_dist->AddText("UIUC")->SetTextFont(62);
  pt_fiber_count_tower_dist->AddText(Form("Mean = %7.3f", h_fiber_count_tower_dist_uiuc->GetMean()));
  pt_fiber_count_tower_dist->AddText("Fudan")->SetTextFont(62);
  pt_fiber_count_tower_dist->AddText(Form("Mean = %7.3f", h_fiber_count_tower_dist_fudan->GetMean()));
  pt_fiber_count_tower_dist->AddText("CIAE")->SetTextFont(62);
  pt_fiber_count_tower_dist->AddText(Form("Mean = %7.3f", h_fiber_count_tower_dist_ciae->GetMean()));
  pt_fiber_count_tower_dist->SetFillColorAlpha(0, 0);
  pt_fiber_count_tower_dist->SetTextAlign(12);
  pt_fiber_count_tower_dist->SetTextSize(0.02);
  pt_fiber_count_tower_dist->SetShadowColor(kWhite);
  pt_fiber_count_tower_dist->Draw();
  cs_fiber_count_tower_dist->SaveAs("emcal_plots/fiber_count_tower_dist.pdf");

  TCanvas *cs_density_dist = new TCanvas();
  hs_density_dist->Draw("NOSTACKB");
  hs_fiber_count_leg->Draw();
  TPaveStats *pt_density_dist = new TPaveStats(0.65, 0.55, 0.875, 0.875, "NDC");
  pt_density_dist->SetTextFont(42);
  pt_density_dist->AddText("EMCal Block Density")->SetTextFont(62);
  pt_density_dist->AddText("UIUC")->SetTextFont(62);
  pt_density_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_density_dist_uiuc->GetFunction("gaus")->GetParameter(1), h_density_dist_uiuc->GetFunction("gaus")->GetParError(1)));
  pt_density_dist->AddText(Form("Mean = %7.3f", h_density_dist_uiuc->GetMean()));
  pt_density_dist->AddText("Fudan")->SetTextFont(62);
  pt_density_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_density_dist_fudan->GetFunction("gaus")->GetParameter(1), h_density_dist_fudan->GetFunction("gaus")->GetParError(1)));
  pt_density_dist->AddText(Form("Mean = %7.3f", h_density_dist_fudan->GetMean()));
  pt_density_dist->AddText("CIAE")->SetTextFont(62);
  pt_density_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_density_dist_ciae->GetFunction("gaus")->GetParameter(1), h_density_dist_ciae->GetFunction("gaus")->GetParError(1)));
  pt_density_dist->AddText(Form("Mean = %7.3f", h_density_dist_ciae->GetMean()));
  pt_density_dist->SetFillColorAlpha(0, 0);
  pt_density_dist->SetTextAlign(12);
  pt_density_dist->SetTextSize(0.02);
  pt_density_dist->SetShadowColor(kWhite);
  pt_density_dist->Draw();
  cs_mpv_block_dist->SaveAs("emcal_plots/mpv_block_dist.pdf");
  cs_density_dist->SaveAs("emcal_plots/density_dist.pdf");

  TCanvas *cs_scint_ratio_dist = new TCanvas();
  hs_scint_ratio_dist->Draw("NOSTACKB");
  hs_fiber_count_leg->Draw();
  TPaveStats *pt_scint_ratio_dist = new TPaveStats(0.65, 0.55, 0.875, 0.875, "NDC");
  pt_scint_ratio_dist->SetTextFont(42);
  pt_scint_ratio_dist->AddText("EMCal Block Scintillation Ratio")->SetTextFont(62);
  pt_scint_ratio_dist->AddText("UIUC")->SetTextFont(62);
  pt_scint_ratio_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_scint_ratio_dist_uiuc->GetFunction("gaus")->GetParameter(1), h_scint_ratio_dist_uiuc->GetFunction("gaus")->GetParError(1)));
  pt_scint_ratio_dist->AddText(Form("Mean = %7.3f", h_scint_ratio_dist_uiuc->GetMean()));
  pt_scint_ratio_dist->AddText("Fudan")->SetTextFont(62);
  pt_scint_ratio_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_scint_ratio_dist_fudan->GetFunction("gaus")->GetParameter(1), h_scint_ratio_dist_fudan->GetFunction("gaus")->GetParError(1)));
  pt_scint_ratio_dist->AddText(Form("Mean = %7.3f", h_scint_ratio_dist_fudan->GetMean()));
  pt_scint_ratio_dist->AddText("CIAE")->SetTextFont(62);
  pt_scint_ratio_dist->AddText(Form("Gaus Mean = %7.3f #pm %5.3f", h_scint_ratio_dist_ciae->GetFunction("gaus")->GetParameter(1), h_scint_ratio_dist_ciae->GetFunction("gaus")->GetParError(1)));
  pt_scint_ratio_dist->AddText(Form("Mean = %7.3f", h_scint_ratio_dist_ciae->GetMean()));
  pt_scint_ratio_dist->SetFillColorAlpha(0, 0);
  pt_scint_ratio_dist->SetTextAlign(12);
  pt_scint_ratio_dist->SetTextSize(0.02);
  pt_scint_ratio_dist->SetShadowColor(kWhite);
  pt_scint_ratio_dist->Draw();
  cs_scint_ratio_dist->SaveAs("emcal_plots/scint_ratio_dist.pdf");


  TFile *chnl_dists_file = new TFile("emcal_plots/histograms.root", "RECREATE");
  chnl_dists_file->WriteObject(h_mpv_block_dist, "h_mpv_block_dist");
  chnl_dists_file->WriteObject(h_mpv_chnl_dist, "h_mpv_chnl_dist");
  chnl_dists_file->WriteObject(h_fiber_count_block_dist, "h_fiber_count_block_dist");
  chnl_dists_file->WriteObject(h_fiber_count_tower_dist, "h_fiber_count_tower_dist");
  chnl_dists_file->WriteObject(hs_mpv_block_dist, "hs_mpv_block_dist");
  chnl_dists_file->WriteObject(hs_mpv_chnl_dist, "hs_mpv_chnl_dist");
  chnl_dists_file->WriteObject(hs_fiber_count_block_dist, "hs_fiber_count_block_dist");
  chnl_dists_file->WriteObject(hs_fiber_count_tower_dist, "hs_fiber_count_tower_dist");
  chnl_dists_file->WriteObject(hs_density_dist, "hs_density_dist");
  chnl_dists_file->WriteObject(hs_scint_ratio_dist, "hs_scint_ratio_dist");
  chnl_dists_file->WriteObject(h_fiber_batches, "h_fiber_batches");
  chnl_dists_file->WriteObject(h_fiber_batch_numbers, "h_fiber_batch_numbers");
} 

/**
 * @brief 
 * 
 * @param channel_lvl 
 * @param mode 
 */
void draw_axes(bool channel_lvl, std::string mode, TCanvas* c, TH2D *h) {
    if (mode == "caroline") {
      if (channel_lvl) {
        c->SetRightMargin(0.125);
        c->SetGrid();
        h->SetAxisRange(0, 128*2 - 1, "X");
        h->SetAxisRange(0, 48*2 - 1, "Y");
        h->GetXaxis()->SetNdivisions(32, false);
        h->GetYaxis()->SetNdivisions(2, false);
        h->GetXaxis()->SetLabelOffset(999.0);
        h->GetYaxis()->SetLabelOffset(999.0);
        h->GetXaxis()->SetTickLength(0);
        h->GetYaxis()->SetTickLength(0);
        h->Draw("COLZ0");
        plot_sector_and_block_labels(true);
      } else {
        c->SetRightMargin(0.125);
        c->SetGrid();
        h->SetAxisRange(0, 128 - 1, "X");
        h->SetAxisRange(0, 48 - 1, "Y");
        h->GetXaxis()->SetNdivisions(32, false);
        h->GetYaxis()->SetNdivisions(2, false);
        h->GetXaxis()->SetLabelOffset(999.0);
        h->GetYaxis()->SetLabelOffset(999.0);
        h->GetXaxis()->SetTickLength(0);
        h->GetYaxis()->SetTickLength(0);
        h->Draw("COLZ0");
        plot_sector_and_block_labels(false);
      }
    } else if (mode == "tim") {
      if (channel_lvl) {
        c->SetGrid();
        c->SetRightMargin(0.125);
        h->SetAxisRange(0, 128*2 - 1, "X");
        h->SetAxisRange(0, 48*2 - 1, "Y");
        h->GetXaxis()->SetLabelSize(0.025);
        h->GetYaxis()->SetLabelSize(0.025);
        h->Draw("COLZ0");
        h->GetXaxis()->SetNdivisions(32, false);
        h->GetYaxis()->SetNdivisions(2, false);
        
        h->GetXaxis()->SetLabelOffset(999.0);
        h->GetXaxis()->SetTickLength(0);

        // draw axes separately (since they don't match the gridlines, which are sector boundaries)
        gPad->Update();
        TGaxis *x_axis = new TGaxis(gPad->GetUxmin(),
                                      gPad->GetUymin(),
                                      gPad->GetUxmax(),
                                      gPad->GetUymin(),
                                      -4,
                                      256 - 4,
                                      16 + 4*100,"N+");            
        x_axis->SetLabelOffset(0.01);
        x_axis->SetLabelFont(42);
        x_axis->SetLabelSize(0.025);
        x_axis->Draw();

        h->GetYaxis()->SetLabelOffset(999.0);
        h->GetYaxis()->SetTickLength(0);
        gPad->Update();
        TGaxis *y_axis = new TGaxis(gPad->GetUxmin(),
                                        gPad->GetUymin(),
                                        gPad->GetUxmin(),
                                        gPad->GetUymax(),
                                        -48,
                                        48,
                                        12 + 2*100,"N-");            
        y_axis->SetLabelOffset(0.01);
        y_axis->SetLabelFont(42);
        y_axis->SetLabelSize(0.025);
        y_axis->Draw();

        plot_sector_labels_debug(true);
      } else {
        c->SetGrid();
        c->SetRightMargin(0.125);
        h->SetAxisRange(0, 128 - 1, "X");
        h->SetAxisRange(0, 48 - 1, "Y");
        h->GetXaxis()->SetLabelSize(0.025);
        h->GetYaxis()->SetLabelSize(0.025);
        h->Draw("COLZ0");
        h->GetXaxis()->SetNdivisions(32, false);
        h->GetYaxis()->SetNdivisions(2, false);
        
        h->GetXaxis()->SetLabelOffset(999.0);
        h->GetXaxis()->SetTickLength(0);

        // draw axes separately (since they don't match the gridlines, which are sector boundaries)
        gPad->Update();
        TGaxis *x_axis = new TGaxis(gPad->GetUxmin(),
                                      gPad->GetUymin(),
                                      gPad->GetUxmax(),
                                      gPad->GetUymin(),
                                      -2,
                                      128 - 2,
                                      16 + 4*100,"N+");            
        x_axis->SetLabelOffset(0.01);
        x_axis->SetLabelFont(42);
        x_axis->SetLabelSize(0.025);
        x_axis->Draw();

        h->GetYaxis()->SetLabelOffset(999.0);
        h->GetYaxis()->SetTickLength(0);
        gPad->Update();
        TGaxis *y_axis = new TGaxis(gPad->GetUxmin(),
                                        gPad->GetUymin(),
                                        gPad->GetUxmin(),
                                        gPad->GetUymax(),
                                        -24,
                                        24,
                                        12 + 2*100,"N-");            
        y_axis->SetLabelOffset(0.01);
        y_axis->SetLabelFont(42);
        y_axis->SetLabelSize(0.025);
        y_axis->Draw();

        plot_sector_labels_debug(false);
      }
    } else {
      throw std::runtime_error(Form("unknown mode '%s'", mode.c_str()));
    }
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
  
  TH2D *h_chnl_mpv = new TH2D("", "sPHENIX EMCal Channel MPV;#phi [Channels];#eta [Channels];MPV", 256, 0, 256, 96, 0, 96);
  TH2D *h_chnl_fiber = new TH2D("", "sPHENIX EMCal Tower Fiber Count;#phi [Towers];#eta [Towers];Fiber Count [%]", 256, 0, 256, 96, 0, 96);

  for (Block &block : all_blocks) {
    auto xy = get_block_loc(block, true_sector_mapping);
    unsigned int x = 2*xy.first + 1;
    unsigned int y = 2*xy.second + 1;

    // auto chnls = block_to_channel(block.block_number - 1);
    double ch0 = block.ch0_mpv;
    double ch1 = block.ch1_mpv;
    double ch2 = block.ch2_mpv;
    double ch3 = block.ch3_mpv;

    if (block.sector % 2 == 0) { // north sector = top half
      if (ch0 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y, ch0);
      }
      if (ch1 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y + 1, ch1);
      }
      if (ch2 > 0) {
        h_chnl_mpv->SetBinContent(x, y, ch2);
      }
      if (ch3 > 0) {
        h_chnl_mpv->SetBinContent(x, y + 1, ch3);
      }
    } else { // south sector = bottom half
      if (ch0 > 0) {
        h_chnl_mpv->SetBinContent(x, y + 1, ch0);
      }
      if (ch1 > 0) {
        h_chnl_mpv->SetBinContent(x, y, ch1);
      }
      if (ch2 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y + 1, ch2);
      }
      if (ch3 > 0) {
        h_chnl_mpv->SetBinContent(x + 1, y, ch3);
      }
    }

    double t1 = block.fiber_t1_count;
    double t2 = block.fiber_t2_count;
    double t3 = block.fiber_t3_count;
    double t4 = block.fiber_t4_count;

    if (block.sector % 2 == 0) { // north sector = top half
      if (t1 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y, t1);
      }
      if (t2 > 0) {
        h_chnl_fiber->SetBinContent(x, y, t2);
      }
      if (t3 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y + 1, t3);
      }
      if (t4 > 0) {
        h_chnl_fiber->SetBinContent(x, y + 1, t4);
      }
    } else { // south sector = bottom half
      if (t1 > 0) {
        h_chnl_fiber->SetBinContent(x, y + 1, t1);
      }
      if (t2 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y + 1, t2);
      }
      if (t3 > 0) {
        h_chnl_fiber->SetBinContent(x, y, t3);
      }
      if (t4 > 0) {
        h_chnl_fiber->SetBinContent(x + 1, y, t4);
      }
    }

    // printf("sector %2d block %2d:\n\tch0: %f\n\tch1: %f\n\tch2: %f\n\tch3: %f\n", block.sector, block.block_number, ch0, ch1, ch2, ch3);
  }

  TCanvas *c_chnl_mpv = new TCanvas();
  h_chnl_mpv->SetMaximum(800);
  draw_axes(true, mode, c_chnl_mpv, h_chnl_mpv);
  c_chnl_mpv->SaveAs("emcal_plots/mpv_chnl_map.pdf");
  
  TCanvas *c_chnl_fiber = new TCanvas();
  h_chnl_fiber->SetMinimum(94);
  h_chnl_fiber->SetMaximum(104);
  draw_axes(true, mode, c_chnl_fiber, h_chnl_fiber);
  c_chnl_fiber->SaveAs("emcal_plots/fiber_count_tower_map.pdf");
}

void plot_block_lvl(std::vector<Block> all_blocks, std::string mode) {
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.5);
  
  TH2D *h_block_mpv = new TH2D("", "sPHENIX EMCal Block MPV;#phi [Blocks];#eta [Blocks];MPV", 128, 0, 128, 48, 0, 48);
  TH2D *h_block_fiber = new TH2D("", "sPHENIX EMCal Block Fiber Count;#phi [Blocks];#eta [Blocks];Fiber Count [%]", 128, 0, 128, 48, 0, 48);
  TH2D *h_block_density = new TH2D("", "sPHENIX EMCal Block Density;#phi [Blocks];#eta [Blocks];Density [g/mL]", 128, 0, 128, 48, 0, 48);
  TH2D *h_block_scint_ratio = new TH2D("", "sPHENIX EMCal Block Scintillation Ratio;#phi [Blocks];#eta [Blocks];Scintillation Ratio", 128, 0, 128, 48, 0, 48);
  
  for (Block &block : all_blocks) {
    auto xy = get_block_loc(block, true_sector_mapping);
    unsigned int x = xy.first + 1;
    unsigned int y = xy.second + 1;

    if (block.mpv > 0) {
      h_block_mpv->SetBinContent(x, y, block.mpv);
    }

    if (block.fiber_count > 0) {
      h_block_fiber->SetBinContent(x, y, block.fiber_count);
    }

    if (block.density > 0) {
      h_block_density->SetBinContent(x, y, block.density);
    }

    if (block.scint_ratio > 0) {
      h_block_scint_ratio->SetBinContent(x, y, block.scint_ratio);
    }

    // printf("sector %2d block %2d:\n\tch0: %f\n\tch1: %f\n\tch2: %f\n\tch3: %f\n", block.sector, block.block_number, ch0, ch1, ch2, ch3);
  }

  TCanvas *c_block_mpv = new TCanvas();
  draw_axes(false, mode, c_block_mpv, h_block_mpv);
  c_block_mpv->SaveAs("emcal_plots/mpv_block_map.pdf");
  
  TCanvas *c_block_fiber = new TCanvas();
  draw_axes(false, mode, c_block_fiber, h_block_fiber);
  c_block_fiber->SaveAs("emcal_plots/fiber_count_block_map.pdf");

  TCanvas *c_block_density = new TCanvas();
  draw_axes(false, mode, c_block_density, h_block_density);
  c_block_density->SaveAs("emcal_plots/density_map.pdf");

  TCanvas *c_block_scint_ratio = new TCanvas();
  h_block_scint_ratio->SetMaximum(3);
  draw_axes(false, mode, c_block_scint_ratio, h_block_scint_ratio);
  c_block_scint_ratio->SaveAs("emcal_plots/scint_ratio_map.pdf");
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

    all_blocks.push_back({
      sector, block_number, dbn, mpv, mpv_err,
      ch0_mpv, ch0_mpv_err, ch1_mpv, ch1_mpv_err, ch2_mpv, ch2_mpv_err, ch3_mpv, ch3_mpv_err,
      density, fiber_count, fiber_t1_count, fiber_t2_count, fiber_t3_count, fiber_t4_count, scint_ratio, fiber_type, fiber_batch, w_powder
    });
    // printf("sector %2d block %2d mpv %f +/- %f\n", sector, block_number, mpv, mpv_err);
    line_num++;
  }

  // TEST SOME THINGS 
  // std::vector<std::string> BASIC_BATCHES;
  // std::vector<FiberBatch> BATCHES;
  // for (const std::string &x : TEST_BATCHES) {
  //   BASIC_BATCHES.push_back(x);
  //   BATCHES.push_back(FiberBatch(x));
  // } 
  
  // std::cout << "unsorted strings:" << std::endl;
  // for (const std::string &x : BASIC_BATCHES) {
  //   std::cout << x << ", ";
  // }
  // std::cout << std::endl;
  // std::sort(BASIC_BATCHES.begin(), BASIC_BATCHES.end());
  // std::cout << "sorted strings:" << std::endl;
  // for (const std::string &x : BASIC_BATCHES) {
  //   std::cout << x << ", ";
  // }
  // std::cout << std::endl;
  // std::cout << "unsorted batches:" << std::endl;
  // for (const FiberBatch &x : BATCHES) {
  //   std::cout << x.str << ", ";
  // }
  // std::sort(BATCHES.begin(), BATCHES.end());
  // std::cout << std::endl;
  // std::cout << "sorted batches:" << std::endl;
  // for (const FiberBatch &x : BATCHES) {
  //   std::cout << x.str << ", ";
  // }
  // std::cout << std::endl;
  
  std::vector<PlotConfig> cfgs = {
    {
      "mpv", "MPV", "", [](Block block){
        double value = block.mpv;
        return std::make_pair(value > 0, value);
      },
      0.0,
      kAzure + 6
    },
    {
      "scint_ratio", "Scintillation Ratio", "", [](Block block){
        double value = block.scint_ratio;
        return std::make_pair(value > 0, value);
      },
      0.7,
      kBlue - 9
    },
    {
      "fiber_count", "Fiber Count", " [%]", [](Block block){
        double value = block.fiber_count;
        return std::make_pair(value > 0, value);
      },
      95.0,
      kGreen - 9
    },
    {
      "density", "Density", " [g/mL]", [](Block block){
        double value = block.density;
        return std::make_pair(value > 0, value);
      },
      8.4,
      kRed - 9
    }
  };

  make_histograms(all_blocks);
  plot_channel_lvl(all_blocks, "tim");
  plot_block_lvl(all_blocks, "tim");
  for (const PlotConfig& cfg : cfgs) {
    plot_helper(all_blocks, cfg);
  }
}
