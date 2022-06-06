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

typedef std::function<std::pair<bool, double>(Block)> value_getter;

typedef struct PlotConfig {
  std::string file_name;
  std::string title;
  std::string units;
  value_getter get_value;
  double plot_min;
} PlotConfig;

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

void plot_sector_and_block_labels() {
  double sector_box_width = 3;
  double sector_box_height = 1.5;
  // odd sectors
  for (unsigned int i = 0; i < 32; i++) {
    double x_center = 2 + 4*i;
    double y_center = 48 + 1.5;
    TPaveText *text = new TPaveText(x_center - sector_box_width/2, y_center - sector_box_height/2, x_center + sector_box_width/2, y_center + sector_box_height/2, "NB");
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
    TPaveText *text = new TPaveText(x_center - sector_box_width/2, y_center - sector_box_height/2, x_center + sector_box_width/2, y_center + sector_box_height/2, "NB");
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
    TPaveText *text1 = new TPaveText(x_center1 - block_box_width/2, y_center1 - block_box_height/2, x_center1 + block_box_width/2, y_center1 + block_box_height/2, "NB");
    text1->SetTextFont(42);
    text1->SetFillColorAlpha(0, 0);
    text1->AddText(Form("%d", i + 1));
    text1->SetTextAlign(32);
    text1->Draw();

    double x_center2 = -2;
    double y_center2 = 24.5 + i;
    TPaveText *text2 = new TPaveText(x_center2 - block_box_width/2, y_center2 - block_box_height/2, x_center2 + block_box_width/2, y_center2 + block_box_height/2, "NB");
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
      TPaveText *text = new TPaveText(x_center - label_box_width/2, y_center - label_box_height/2, x_center + label_box_width/2, y_center + label_box_height/2, "NB");
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
      TPaveText *text = new TPaveText(x_center - label_box_width/2, y_center - label_box_height/2, x_center + label_box_width/2, y_center + label_box_height/2, "NB");
      text->SetTextFont(82);
      text->SetFillColorAlpha(0, 0);
      text->AddText(north_labels[j].c_str());
      text->SetTextAlign(22);
      text->Draw();
    }
  }
}

void plot_helper(std::vector<Block> all_blocks, PlotConfig cfg) {
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
    unsigned int len = strlen(block.dbn.c_str());
    if (len > max_strlen) {
      max_strlen = len;
      // printf("'%s' IS NOW THE LONGEST AT %u\n", block.dbn.c_str(), len);
    }
    text->AddText(Form("%s", block.dbn.c_str()));
    text->SetTextAlign(22);
    text->Draw();
  }
  // printf("MAX STRLEN IS %u\n", max_strlen);

  c1->SaveAs(Form("emcal_plots/plot_colz_dbn_%s.pdf", cfg.file_name.c_str()));
  
  gStyle->SetImageScaling(2.0);
  // h->SetMinimum(0);
  
  h->GetZaxis()->SetTitleOffset(1.2);

  TCanvas* c2 = new TCanvas("c2", "", 900, 500);
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.5);
  h->SetLineColorAlpha(kBlack, 1.0);
  h->SetTitle(title.c_str());
  h->Draw("LEGO2 0");
  c2->SaveAs(Form("emcal_plots/plot_lego_%s.png", cfg.file_name.c_str()));

  TCanvas* c3 = new TCanvas("c3", "", 900, 500);
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.5);
  h->SetLineColorAlpha(kBlack, 1.0);
  h->SetFillColorAlpha(kBlue, 1.0);
  h->SetTitle(title.c_str());
  h->Draw("LEGO1 CYL 0");
  c3->SaveAs(Form("emcal_plots/plot_cyl_%s.png", cfg.file_name.c_str()));
  
  TCanvas* c4 = new TCanvas("c4", "", 900, 500);
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.5);
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
}

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


    // dbn.erase(std::remove(dbn.begin(), dbn.end(), '\n'), dbn.end());
    // dbn.erase(std::remove(dbn.begin(), dbn.end(), '\r'), dbn.end());
    // dbn.erase(std::remove(dbn.begin(), dbn.end(), ' '), dbn.end());

    all_blocks.push_back({sector, block_number, dbn, mpv, density, fiber_count, fiber_t1_count, fiber_t2_count, fiber_t3_count, fiber_t4_count, scint_ratio});
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
      8.5
    }
  };

  for (const PlotConfig& cfg : cfgs) {
    plot_helper(all_blocks, cfg);
  }
}