#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <sstream>
#include <TFile.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <Fit/Fitter.h>
#include <HFitInterface.h>
#include <Fit/BinData.h>
//#include <Fit/UnBinData.h>

#include <cmath>
#include <stdexcept>
#include <queue>
#include <sys/stat.h>
#include <pthread.h>

#define SAVE_PLOTS false

// TODO: multithreading

pthread_mutex_t queue_lock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
  int channel_num;
  TH1D *hist;
} task_t;

std::queue<task_t> tasks;

char *file_prefix = NULL;
Double_t gap_spacing[384][2];
Double_t first_ampl[384][2];
Double_t first_mean[384][2];
Double_t first_sigma[384][2];
Double_t other_ampl[384][5][2];
Double_t other_sigma[384][5][2];

bool fit_success[384] = {false};

/**
 * The function used for fitting
 */
Double_t fitf (Double_t *x, const Double_t *par) {
  Double_t gaussians = 0;
  gaussians += par[1] * TMath::Exp(-TMath::Power(x[0]-par[2],2)/(2*par[3]*par[3]));
  for (int i = 0; i < 5; i++) {
    // par[2*i+6] = par[6];
    gaussians += par[2*i+4] * TMath::Exp(-TMath::Power(x[0]-((i+1)*par[0]+par[2]),2)/(2*par[2*i+5]*par[2*i+5]));
  }
  //par[0] = gap spacing
  //par[1] = first peak amplitude
  //par[2] = first peak mean
  //par[3] = first peak sigma
  //par[4,6,8,10,12] = other peak amplitudes
  //par[5,7,9,11,13] = other peak sigmas

  return gaussians;
}

class six_gauss : public ROOT::Math::IParamFunction { 
  public:
    void SetParameters(const Double_t* p) { std::copy(p,p+NPar(),fp);}
    const Double_t* Parameters() const { return fp; }
    ROOT::Math::IGenFunction* Clone() const { 
      six_gauss* g = new six_gauss(); 
      g->SetParameters(fp);
      return g;
    };
    unsigned int NPar() const { return 14; }
  private:
    Double_t DoEvalPar(Double_t x, const Double_t* par) const { 
      return fitf(&x, fp);
    }
    Double_t fp[14];  
};

void *parallel_fit(void *ptr) {
  while (1) {
    pthread_mutex_lock(&queue_lock);
    if (tasks.empty()) {
      pthread_mutex_unlock(&queue_lock);
      break;
    }
    // else there are tasks :)
    task_t my_task = tasks.front();
    tasks.pop();
    pthread_mutex_unlock(&queue_lock);

    TH1D *my_hist = my_task.hist;
    int i = my_task.channel_num;
    pthread_t my_id = pthread_self();

    printf("(thread %lu) channel %i: started...\n", my_id, i);

    if (my_hist->GetEntries() < 1e2) {
      printf("(thread %lu) channel %i: too few entries!\n", my_id, i);
      continue;
    }

    /*
    TF1 *my_fit = new TF1(Form("f_%i", i), fitf, 1200.0, 2400.0, 14);


    my_fit->SetParameter(0, 28.0); // par[0] = gap spacing
    my_fit->SetParLimits(0, 18.0, 36.0);
    my_fit->SetParameter(1, 1e4); // par[1] = first peak amplitude
    my_fit->SetParLimits(1, 10.0, 1e6);
    my_fit->SetParameter(2, 1500.0); // par[2] = first peak mean
    my_fit->SetParLimits(2, 1200.0, 1800.0);
    my_fit->SetParameter(3, 8.0); // par[3] = first peak sigma
    my_fit->SetParLimits(3, 5.0, 10.0);
    for (int j = 0; j < 5; j++) {
      my_fit->SetParameter(2*j+4, 1e4); // par[4,6,8,10,12] = other peak amplitudes
      my_fit->SetParLimits(2*j+4, 10.0, 1e6);
      my_fit->SetParameter(2*j+5, 8.0); // par[5,7,9,11,13] = other peak sigmas
      my_fit->SetParLimits(2*j+5, 5.0, 10.0);
    }

    
    my_fit->SetLineWidth(2);

    TFitResultPtr sp_fit = my_hist->Fit(my_fit, "Q S", "", 1200.0, 2000.0);

    printf("thread %lu successfully fit\n", my_id);
    
    if (sp_fit->IsEmpty()) {
      printf("rejected channel %i for null or empty sp_fit\n", i);
      continue;
    }
    */

    ROOT::Fit::Fitter fitter;
    six_gauss g;
    fitter.Config().SetMinimizer("GSLMultiMin");
    fitter.SetFunction(g);
    ROOT::Fit::BinData d; 
    ROOT::Fit::FillData(d, my_hist);

    fitter.Config().ParSettings(0).SetValue(28.0); // par[0] = gap spacing
    fitter.Config().ParSettings(0).SetLimits(18.0, 36.0);
    fitter.Config().ParSettings(1).SetValue(1e4); // par[1] = first peak amplitude
    fitter.Config().ParSettings(1).SetLimits(10.0, 1e6);
    fitter.Config().ParSettings(2).SetValue(1500.0); // par[2] = first peak mean
    fitter.Config().ParSettings(2).SetLimits(1200.0, 1800.0);
    fitter.Config().ParSettings(3).SetValue(8.0); // par[3] = first peak sigma
    fitter.Config().ParSettings(3).SetLimits(5.0, 10.0);
    for (int j = 0; j < 5; j++) {
      fitter.Config().ParSettings(2*j+4).SetValue(1e4); // par[4,6,8,10,12] = other peak amplitudes
      fitter.Config().ParSettings(2*j+4).SetLimits(10.0, 1e6);
      fitter.Config().ParSettings(2*j+5).SetValue(8.0); // par[5,7,9,11,13] = other peak sigmas
      fitter.Config().ParSettings(2*j+5).SetLimits(5.0, 10.0);
    }

    bool sp_fit = fitter.Fit(d);

    bool success = fitter.Result().IsValid();

    if (!success) {
      printf("(thread %lu) channel %i: fit did not converge!\n", my_id, i);
      continue;
    } else {
      fit_success[i] = true;
    }

    //printf("thread %lu fitted channel %i \n", my_id, i);

    gap_spacing[i][0] = fitter.Result().Parameter(0); //par[0] = gap spacing
    gap_spacing[i][1] = fitter.Result().ParError(0); //par[0] = gap spacing
    first_ampl[i][0] = fitter.Result().Parameter(1); //par[1] = first peak amplitude
    first_ampl[i][1] = fitter.Result().ParError(1); //par[1] = first peak amplitude
    first_mean[i][0] = fitter.Result().Parameter(2); //par[2] = first peak mean
    first_mean[i][1] = fitter.Result().ParError(2); //par[2] = first peak mean
    first_sigma[i][0] = fitter.Result().Parameter(3); //par[3] = first peak sigma
    first_sigma[i][1] = fitter.Result().ParError(3); //par[3] = first peak sigma
    for (int j = 0; j < 5; j++) {
      other_ampl[i][j][0] = fitter.Result().Parameter(2*j+4);
      other_ampl[i][j][1] = fitter.Result().ParError(2*j+4);
      other_sigma[i][j][0] = fitter.Result().Parameter(2*j+5);
      other_sigma[i][j][1] = fitter.Result().ParError(2*j+5);
    }

    //printf("thread %lu saved fit params for channel %i \n", my_id, i);

    

    /*
    TF1 *landau = new TF1("landau", "landau", 0.0, 200.0);
    landau->SetParameter(0, landau_amplitudes[i]);
    landau->SetParameter(1, landau_mpvs[i]);
    landau->SetParameter(2, landau_sigmas[i]);
    landau->SetLineStyle(1); landau->SetLineColor(kGreen); landau->SetLineWidth(1);
    
    //TLine* landau_ampl = new TLine(0.0, landau_amplitudes[i], 200.0, landau_amplitudes[i]);
    //landau_ampl->SetLineColor(kGreen);
    //landau_ampl->SetLineStyle(2);
    //landau_ampl->SetLineWidth(1);
    

    TF1 *gaus1 = new TF1("gaus1", "gaus", 0.0, 200.0);
    gaus1->SetParameter(0, first_gauss_amplitude[i]);
    gaus1->SetParameter(1, first_gauss_mean[i]);
    gaus1->SetParameter(2, first_gauss_sigma[i]);
    gaus1->SetLineStyle(1); gaus1->SetLineColor(kBlue); gaus1->SetLineWidth(1);

    //TLine* gauss1_ampl = new TLine(0.0, first_gauss_amplitude[i], 200.0, first_gauss_amplitude[i]);
    //gauss1_ampl->SetLineColor(kBlue);
    //gauss1_ampl->SetLineStyle(2);
    //gauss1_ampl->SetLineWidth(1);

    TF1 *gaus2 = new TF1("gaus2", "gaus", 0.0, 200.0);
    gaus2->SetParameter(0, other_gauss_amplitudes[i][0]);
    gaus2->SetParameter(1, first_gauss_mean[i] + sp_gap);
    gaus2->SetParameter(2, other_gauss_sigmas[i][0]);
    gaus2->SetLineStyle(1); gaus2->SetLineColor(kBlue); gaus2->SetLineWidth(1);

    TF1 *gaus3 = new TF1("gaus3", "gaus", 0.0, 200.0);
    gaus3->SetParameter(0, other_gauss_amplitudes[i][1]);
    gaus3->SetParameter(1, first_gauss_mean[i] + 2*sp_gap);
    gaus3->SetParameter(2, other_gauss_sigmas[i][1]);
    gaus3->SetLineStyle(1); gaus3->SetLineColor(kBlue); gaus3->SetLineWidth(1);

    TF1 *gaus4 = new TF1("gaus4", "gaus", 0.0, 200.0);
    gaus4->SetParameter(0, other_gauss_amplitudes[i][2]);
    gaus4->SetParameter(1, first_gauss_mean[i] + 3*sp_gap);
    gaus4->SetParameter(2, other_gauss_sigmas[i][2]);
    gaus4->SetLineStyle(1); gaus4->SetLineColor(kBlue); gaus4->SetLineWidth(1);
    
    TF1 *gaus5 = new TF1("gaus5", "gaus", 0.0, 200.0);
    gaus5->SetParameter(0, other_gauss_amplitudes[i][3]);
    gaus5->SetParameter(1, first_gauss_mean[i] + 4*sp_gap);
    gaus5->SetParameter(2, other_gauss_sigmas[i][3]);
    gaus5->SetLineStyle(1); gaus5->SetLineColor(kBlue); gaus5->SetLineWidth(1);

    */
    
    if (SAVE_PLOTS) {
      double max_bin = my_hist->GetMaximumBin();

      my_hist->GetXaxis()->SetRangeUser(max_bin - 50.0, max_bin + 200.0);    
      //my_hist->GetXaxis()->SetRangeUser(1200.0, 2000.0);    

      TCanvas* c1 = new TCanvas(Form("c%i_%lu", i, my_id), "", 700, 500);

      Double_t means[6];
      means[0] = first_mean[i][0];
      for (int j = 1; j < 6; j++) {
        means[j] = means[0] + j * gap_spacing[i][0];
      }

      my_hist->Draw();

      gPad->SetLogy();
      gPad->Update();

      //printf("uymin %f, uymax %f\n", c1->GetUymin(), c1->GetUymax());

      for (int j = 0; j < 6; j++) {
        TLine* line = new TLine(means[j], pow(10.0, c1->GetUymin()), means[j], pow(10.0, c1->GetUymax()));
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("SAME");
      }

      //landau->Draw("SAME");
      //landau_ampl->Draw("SAME");
      //gaus1->Draw("SAME");
      //gauss1_ampl->Draw("SAME");
      //gaus2->Draw("SAME");
      //gaus3->Draw("SAME");
      //gaus4->Draw("SAME");
      //gaus5->Draw("SAME");

      /*
      Double_t yndc = 1e5;
      TLine* p1_line = new TLine(p1, 0.0, p1, yndc);
      p1_line->SetLineColor(kRed);
      p1_line->SetLineStyle(2);
      p1_line->SetLineWidth(2);
      p1_line->Draw("SAME");
      TLine* p2_line = new TLine(p2, 0.0, p2, yndc);
      p2_line->SetLineColor(kRed);
      p2_line->SetLineStyle(2);
      p2_line->SetLineWidth(2);
      p2_line->Draw("SAME");
      */

      /*
      TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
      legend->AddEntry(f_landaugaus, Form("SP Gap: %.3f", sp_gap), "l");
      legend->AddEntry("", Form("ChiSqr/NDF: %.3f", chisqr_ndfs[i]));
      legend->Draw();
      */

      c1->SaveAs(Form("%s/channel_%i.png", file_prefix, i));
    }

    //printf("thread %lu finished channel %i and gap is %f ± %f\n", my_id, i, gap_spacing[i][0], gap_spacing[i][1]);
  }
  printf("thread %lu EXITED!!!\n", pthread_self());
  return NULL;
}

void fit_no_pedestal_multithread(int run_num, int n_threads) {
  if (n_threads < 1) {
    throw std::runtime_error("n_threads should be >= 1");
  }
  ROOT::EnableThreadSafety();
  //ROOT::EnableImplicitMT();
  int n_channels = 0;
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);   //make the plot list all the fit information
  
  TFile *hist_file = TFile::Open(Form("./all_runs/qa_output_nopedestal_000%i/histograms.root", run_num));
  if (!hist_file) {
    throw std::runtime_error("unable to open histogram file");
  }

  file_prefix = strdup(Form("./no_pedestal_hists_%i", run_num));

  // retrieve all histograms from the root file
  TH1D *h_alladc[384];
  for (int i = 0; i < 384;i++) {
    h_alladc[i] = (TH1D*) hist_file->Get(Form("h_alladc_%i", i));
  }
  
  TF1 *f_singlepixels[384];


  for (int i = 0; i < 384; i++) {
    task_t t = {i, h_alladc[i]};
    tasks.push(t);
  }
  
  // execute multithreaded code:

  pthread_t pids[n_threads];
  for (int i = 0; i < n_threads; i++) {
    if (pthread_create(pids + i, NULL, parallel_fit, NULL)) {
      throw std::runtime_error("pthread_create failed!");
    }
  }

  // wait for theads!
  for (int i = 0; i < n_threads; i++) {
    if (pthread_join(pids[i], NULL) != 0) {
      // pthread_join failed!
      throw std::runtime_error("pthread_join failed!");
    }
  }



  // threads complete

  // now we can single-threadedly save all histograms

  if (!SAVE_PLOTS) {
    for (int i = 0; i < 384; i++) {
      char *hist_file_name = Form("%s/channel_%i.png", file_prefix, i);
      if (!fit_success[i]) {
        remove(hist_file_name);
        continue;
      }
      TH1D *my_hist = h_alladc[i];
      double max_bin = my_hist->GetMaximumBin();

      my_hist->GetXaxis()->SetRangeUser(max_bin - 50.0, max_bin + 200.0);    
      //my_hist->GetXaxis()->SetRangeUser(1200.0, 2000.0);    

      TCanvas* c1 = new TCanvas(Form("c%i", i), "", 700, 500);

      if (fit_success[i]) {
        Double_t means[6];
        means[0] = first_mean[i][0];
        for (int j = 1; j < 6; j++) {
          means[j] = means[0] + j * gap_spacing[i][0];
        }

        my_hist->Draw();

        gPad->SetLogy();
        gPad->Update();

        //printf("uymin %f, uymax %f\n", c1->GetUymin(), c1->GetUymax());

        for (int j = 0; j < 6; j++) {
          TLine* line = new TLine(means[j], pow(10.0, c1->GetUymin()), means[j], pow(10.0, c1->GetUymax()));
          line->SetLineColor(kBlack);
          line->SetLineStyle(2);
          line->SetLineWidth(1);
          line->Draw("SAME");
        }



        //landau->Draw("SAME");
        //landau_ampl->Draw("SAME");
        //gaus1->Draw("SAME");
        //gauss1_ampl->Draw("SAME");
        //gaus2->Draw("SAME");
        //gaus3->Draw("SAME");
        //gaus4->Draw("SAME");
        //gaus5->Draw("SAME");

        /*
        Double_t yndc = 1e5;
        TLine* p1_line = new TLine(p1, 0.0, p1, yndc);
        p1_line->SetLineColor(kRed);
        p1_line->SetLineStyle(2);
        p1_line->SetLineWidth(2);
        p1_line->Draw("SAME");
        TLine* p2_line = new TLine(p2, 0.0, p2, yndc);
        p2_line->SetLineColor(kRed);
        p2_line->SetLineStyle(2);
        p2_line->SetLineWidth(2);
        p2_line->Draw("SAME");
        */

        /*
        TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
        legend->AddEntry(f_landaugaus, Form("SP Gap: %.3f", sp_gap), "l");
        legend->AddEntry("", Form("ChiSqr/NDF: %.3f", chisqr_ndfs[i]));
        legend->Draw();
        */
      }

      c1->SaveAs(hist_file_name);
    }
  }

  FILE *outfile = fopen(Form("%s/no_pedestal_params.csv", file_prefix), "w+");
  fprintf(outfile, "channel, gap spacing, gap spacing err, 1st peak amplitude, 1st peak amplitude err, 1st peak mean, 1st peak mean err, 1st peak sigma, 1st peak sigma err, 2nd peak amplitudes, 2nd peak amplitudes err, 2nd peak sigmas, 2nd peak sigmas err, 3rd peak amplitudes, 3rd peak amplitudes err, 3rd peak sigmas, 3rd peak sigmas err, 4th peak amplitudes, 4th peak amplitudes err, 4th peak sigmas, 4th peak sigmas err, 5th peak amplitudes, 5th peak amplitudes err, 5th peak sigmas, 5th peak sigmas err, 6th peak amplitudes, 6th peak amplitudes err, 6th peak sigmas, 6th peak sigmas err");
  for (short i = 0; i < 384; i++) {
    fprintf(outfile, "\n%i,%f,%f,%f,%f,%f,%f,%f,%f", i, gap_spacing[i][0], gap_spacing[i][1], first_ampl[i][0], first_ampl[i][1], first_mean[i][0], first_mean[i][1], first_sigma[i][0], first_sigma[i][1]);
    for (short j = 0; j < 5; j++) {
      fprintf(outfile, ",%f,%f,%f,%f", other_ampl[i][j][0], other_ampl[i][j][1], other_sigma[i][j][0], other_sigma[i][j][1]);
    }
  }
  fclose(outfile);


  // compute block gaps
  Double_t blocks[96];
  Double_t blocks_err[96];
  Double_t block_gaps[96];
  Double_t block_gaps_err[96];
  for (int i = 0; i < 96; i++) {
    blocks[i] = i + 0.5;
    blocks_err[i] = 0;
    block_gaps[i] = 0.0;
    block_gaps_err[i] = 0.0;
    for (int j = 0; j < 4; j++) {
      block_gaps[i] += gap_spacing[4*i + j][0];
      block_gaps_err[i] += gap_spacing[4*i + j][1];
    }
    block_gaps[i] /= 4;
    block_gaps_err[i] /= 4;
  }

  TGraphErrors *gap_graph = new TGraphErrors(96, blocks, block_gaps, blocks_err, block_gaps_err);
  gap_graph->SetTitle(Form("Run %i Single Pixel Gaps (No Pedestal)", run_num));
  gap_graph->GetXaxis()->SetTitle("Block");
  gap_graph->GetYaxis()->SetTitle("Single Pixel Gap");
  gap_graph->GetXaxis()->SetLimits(0.0, 96.0);
  gap_graph->SetMarkerStyle(21);
  gap_graph->SetMarkerColor(kBlue);
  TCanvas *c2 = new TCanvas("c2", "", 700, 500);
  gap_graph->Draw("AP");
  gPad->Update();
  for (int j = 1; j < 6; j++) {
    TLine* line = new TLine(j*16, c2->GetUymin(), j*16, c2->GetUymax());
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");
  }
  c2->SetGrid();
  c2->SaveAs(Form("%s/gaps.png", file_prefix));

  //printf("max chisqr_ndf: %f (channel %i); min chisqr_ndf: %f (channel %i)\n", max_chisqr_ndf, max_chisqr_idx, min_chisqr_ndf, min_chisqr_idx);

  // compute averages for IBs 0-2 and 3-5
  Double_t avg_gap = 0;
  for (int i = 0; i < 192; i++) {
    avg_gap += gap_spacing[i][0];
  }
  Double_t ib0_2_avg_gap = avg_gap / 192;
  avg_gap = 0;
  for (int i = 192; i < 384; i++) {
    avg_gap += gap_spacing[i][0];
  }
  Double_t ib3_5_avg_gap = avg_gap / 192;
  printf("IBs 0-2: avg sp gap = %f; IBs 3-5 avg sp gap = %f\n", ib0_2_avg_gap, ib3_5_avg_gap);
  free(file_prefix);
}