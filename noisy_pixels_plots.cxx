// MFT study of the number of noisy pixels
// with respect to the last SB stop and last GO_READY
// David Grund, 2024

// cpp headers
#include <sstream>
// root headers
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLegend.h"
// custom headers
#include "utilities.h"

template<typename TH>
void format_histo (TH* h)
{
  h->SetMarkerColor(kBlack);
  h->SetMarkerSize(1.2);
  h->SetMarkerStyle(kFullCircle);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  // x-axis
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetLabelSize(0.04);
  // y-axis
  h->GetYaxis()->SetTitle("#Noisy pixels");
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetYaxis()->SetLabelSize(0.03);
  return;
}

void noisy_pix_correlation_with_delays (std::vector<noise>* noise_runs, std::string folder = "")
{
  // axes and ranges
  const int N_2d = 5;
  float ax_5bins[N_2d+1] = {0, 1, 2, 3, 4, 5};
  float ax_sb_2d[N_2d+1] = {0, 5, 10, 30, 60, 1e4};
  float ax_re_2d[N_2d+1] = {0, 5, 10, 60, 600, 1e4};

  const int N_sb = 8;
  double ax_8bins[N_sb+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  double ax_sb_1d[N_sb+1] = {0, 5, 10, 15, 30, 60, 120, 480, 1e4};

  const int N_re = 7;
  double ax_7bins[N_re+1] = {0, 1, 2, 3, 4, 5, 6, 7};
  double ax_re_1d[N_re+1] = {0, 5, 10, 60, 120, 360, 720, 1e4};

  // histograms
  TH2F* h2_true_total = new TH2F("", "", N_2d, &ax_sb_2d[0], N_2d, &ax_re_2d[0]);
  TH2I* h2_true_count = new TH2I("", "", N_2d, &ax_sb_2d[0], N_2d, &ax_re_2d[0]);

  TH2F* h2_bins_total = new TH2F("", Form("#Noisy pixels (#color[4]{#runs, total: %lu})", noise_runs->size()), N_2d, &ax_5bins[0], N_2d, &ax_5bins[0]);
  TH2I* h2_bins_count = new TH2I("", "#Runs", N_2d, &ax_5bins[0], N_2d, &ax_5bins[0]);

  TH1F* h1_true_total_sb = new TH1F("", "", N_sb, &ax_sb_1d[0]);
  TH1I* h1_true_count_sb = new TH1I("", "", N_sb, &ax_sb_1d[0]);
  TH1F* h1_bins_total_sb = new TH1F("", "#Noisy pixels", N_sb, &ax_8bins[0]);

  TH1F* h1_true_total_re = new TH1F("", "", N_re, &ax_re_1d[0]);
  TH1I* h1_true_count_re = new TH1I("", "", N_re, &ax_re_1d[0]);
  TH1F* h1_bins_total_re = new TH1F("", "#Noisy pixels", N_re, &ax_7bins[0]);

  // fill the true histograms
  for (auto n : *noise_runs) {
    float delay_sb = (std::get<0>(n.timestamps) - std::get<1>(n.timestamps)) / 60; // in minutes
    float delay_re = (std::get<0>(n.timestamps) - std::get<2>(n.timestamps)) / 60; // in minutes
    h2_true_total->Fill(delay_sb, delay_re, std::get<0>(n.noisy_pixels));
    h2_true_count->Fill(delay_sb, delay_re);
    h1_true_total_sb->Fill(delay_sb, std::get<0>(n.noisy_pixels));
    h1_true_count_sb->Fill(delay_sb);
    h1_true_total_re->Fill(delay_re, std::get<0>(n.noisy_pixels));
    h1_true_count_re->Fill(delay_re);
  }

  // fill the bin histograms
  for (int x = 1; x <= N_2d; x++) {
    for (int y = 1; y <= N_2d; y++) {
      int entries = h2_true_count->GetBinContent(x, y);
      if (entries > 0) {
        h2_bins_total->SetBinContent(x, y, h2_true_total->GetBinContent(x, y) / entries);
        h2_bins_count->SetBinContent(x, y, entries);
      }
    }
  }
  for (int x = 1; x <= N_sb; x++) {
    int entries = h1_true_count_sb->GetBinContent(x);
    if (entries > 0) h1_bins_total_sb->SetBinContent(x, h1_true_total_sb->GetBinContent(x) / entries);
  }
  for (int x = 1; x <= N_re; x++) {
    int entries = h1_true_count_re->GetBinContent(x);
    if (entries > 0) h1_bins_total_re->SetBinContent(x, h1_true_total_re->GetBinContent(x) / entries);
  }

  // 2d correlation plot
  TCanvas* c1 = new TCanvas("", "", 900, 700);
  set_margins(c1, 0.07, 0.12, 0.10, 0.145);

  // h2_total
  // axis labels
  for(int i = 1; i < N_2d; i++) {
    h2_bins_total->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_sb_2d[i-1], ax_sb_2d[i]));
    h2_bins_total->GetYaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_re_2d[i-1], ax_re_2d[i]));
  }
  h2_bins_total->GetXaxis()->SetBinLabel(N_2d, Form(">= %.0f", ax_sb_2d[N_2d-1]));
  h2_bins_total->GetYaxis()->SetBinLabel(N_2d, Form(">= %.0f", ax_re_2d[N_2d-1]));
  // x-axis
  h2_bins_total->GetXaxis()->SetTitle("Time since last SB stop [min]");
  h2_bins_total->GetXaxis()->SetTitleSize(0.035);
  h2_bins_total->GetXaxis()->SetTitleOffset(1.2);
  h2_bins_total->GetXaxis()->SetLabelSize(0.04);
  // y-axis
  h2_bins_total->GetYaxis()->SetTitle("Time since last GO_READY [min]");
  h2_bins_total->GetYaxis()->SetTitleSize(0.035);
  h2_bins_total->GetYaxis()->SetTitleOffset(2.2);
  h2_bins_total->GetYaxis()->SetLabelSize(0.04);
  // z-axis
  h2_bins_total->GetZaxis()->SetLabelSize(0.03);
  h2_bins_total->GetZaxis()->SetRangeUser(7500, 9000);
  h2_bins_total->SetBarOffset(+0.08);

  // h2_count
  h2_bins_count->SetBarOffset(-0.08);
  h2_bins_count->SetMarkerColor(kBlue);

  // draw the histograms and print the canvas
  h2_bins_total->Draw("colz text");
  h2_bins_count->Draw("same text");
  c1->Print(Form("%scorr_2d.pdf", folder.data()));

  TCanvas* c2 = new TCanvas("", "", 900, 700);
  set_margins(c2, 0.07, 0.02, 0.10, 0.12);
  format_histo(h1_bins_total_sb);
  h1_bins_total_sb->GetXaxis()->SetTitle("Time since last SB stop [min]");
  for(int i = 1; i < N_sb; i++) h1_bins_total_sb->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_sb_1d[i-1], ax_sb_1d[i]));
  h1_bins_total_sb->GetXaxis()->SetBinLabel(N_sb, Form(">= %.0f", ax_sb_1d[N_sb-1]));
  h1_bins_total_sb->Draw("e0");
  c2->Print(Form("%scorr_sb_stop.pdf", folder.data()));

  TCanvas* c3 = new TCanvas("", "", 900, 700);
  set_margins(c3, 0.07, 0.02, 0.10, 0.12);
  format_histo(h1_bins_total_re);
  h1_bins_total_re->GetXaxis()->SetTitle("Time since last GO_READY [min]");
  for(int i = 1; i < N_re; i++) h1_bins_total_re->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_re_1d[i-1], ax_re_1d[i]));
  h1_bins_total_re->GetXaxis()->SetBinLabel(N_re, Form(">= %.0f", ax_re_1d[N_re-1]));
  h1_bins_total_re->Draw("e0");
  c3->Print(Form("%scorr_go_ready.pdf", folder.data()));

  return;
}

void noisy_pix_trend (std::vector<noise>* noise_runs, std::string plot_opt,
  std::tuple<float, float, float, float> delays = {0, 1e4, 0, 1e4}, // sb_min, sb_max, re_min, re_max
  std::tuple<float, float> force_ranges = {-1, -1},
  std::string folder = "")
{
  TGraph* gr_trend_total = new TGraph();
  TGraph* gr_trend_new = new TGraph();
  TGraph* gr_trend_disapp = new TGraph();

  TF1 *f_total = new TF1("f", "[1]*x + [0]");

  std::vector<int> run_numbers;
  for (auto n : *noise_runs) {
    float delay_sb = (std::get<0>(n.timestamps) - std::get<1>(n.timestamps) ) / 60; // in minutes
    float delay_re = (std::get<0>(n.timestamps) - std::get<2>(n.timestamps) ) / 60; // in minutes
    if(  (delay_sb >= std::get<0>(delays)) && (delay_sb < std::get<1>(delays)) 
      && (delay_re >= std::get<2>(delays)) && (delay_re < std::get<3>(delays)) )
    {
      run_numbers.push_back(n.run);
      gr_trend_total->AddPoint(std::get<0>(n.timestamps), std::get<0>(n.noisy_pixels));
      gr_trend_new->AddPoint(std::get<0>(n.timestamps), std::get<1>(n.noisy_pixels));
      gr_trend_disapp->AddPoint(std::get<0>(n.timestamps), std::get<2>(n.noisy_pixels));
    }
  }

  // fits
  gr_trend_total->Fit(f_total);

  gr_trend_total->SetLineColor(kBlack);
  gr_trend_total->SetLineStyle(1);
  gr_trend_total->SetLineWidth(1);
  gr_trend_total->SetMarkerColor(kBlack);
  gr_trend_total->SetMarkerStyle(kOpenCircle);
  gr_trend_total->SetMarkerSize(0.8);

  gr_trend_new->SetLineColor(kBlue);
  gr_trend_new->SetLineStyle(1);
  gr_trend_new->SetLineWidth(1);
  gr_trend_new->SetMarkerColor(kBlue);
  gr_trend_new->SetMarkerStyle(kOpenSquare);
  gr_trend_new->SetMarkerSize(0.8);

  gr_trend_disapp->SetLineColor(kBlue);
  gr_trend_disapp->SetLineStyle(2);
  gr_trend_disapp->SetLineWidth(1);
  gr_trend_disapp->SetMarkerColor(kBlue);
  gr_trend_disapp->SetMarkerStyle(kFullTriangleUp);
  gr_trend_disapp->SetMarkerSize(0.8);

  float margin_t = 0.02;
  float margin_r = 0.10;
  float margin_b = 0.15;
  float margin_l = 0.12;

  TCanvas* c = new TCanvas("", "", 900, 700);
  c->cd();
  // first pad
  TPad* p1 = new TPad("", "", 0., 0., 1., 1.);
  set_margins(p1, margin_t, margin_r, margin_b, margin_l);
  p1->Draw();
  p1->cd();
  gr_trend_total->GetYaxis()->SetTitle("#Noisy pixels");
  gr_trend_total->GetYaxis()->SetTitleOffset(1.65);
  gr_trend_total->GetHistogram()->SetNdivisions(510, "X");
  gr_trend_total->GetXaxis()->SetMaxDigits(6);
  gr_trend_total->GetXaxis()->SetLabelSize(0);
  gr_trend_total->GetXaxis()->SetTitleSize(0);
  float set_ymin = gr_trend_total->GetHistogram()->GetMinimum() * 0.75;
  float set_ymax = gr_trend_total->GetHistogram()->GetMaximum();
  if (std::get<0>(force_ranges) > 0) set_ymin = std::get<0>(force_ranges);
  if (std::get<1>(force_ranges) > 0) set_ymax = std::get<1>(force_ranges);
  float set_dy = (set_ymax - set_ymin) / (1. - margin_t - margin_b);
  gr_trend_total->GetYaxis()->SetRangeUser(set_ymin, set_ymax);
  gr_trend_total->Draw(Form("A%s", plot_opt.data()));
  // font sizes and styles
  Style_t tfont = gr_trend_total->GetHistogram()->GetYaxis()->GetTitleFont();
  Style_t lfont = gr_trend_total->GetHistogram()->GetYaxis()->GetLabelFont();
  float tsize = gr_trend_total->GetHistogram()->GetYaxis()->GetTitleSize();
  float lsize = gr_trend_total->GetHistogram()->GetYaxis()->GetLabelSize();
  // ranges
  float xmin = gr_trend_total->GetXaxis()->GetXmin();
  float xmax = gr_trend_total->GetXaxis()->GetXmax();
  float dx = (xmax - xmin) / (1. - margin_l - margin_r);
  float ymin = 0; 
  float ymax_new = gr_trend_new->GetHistogram()->GetMaximum() * 1.8;
  float ymax_disapp = gr_trend_disapp->GetHistogram()->GetMaximum() * 1.8;
  float ymax = ymax_new > ymax_disapp ? ymax_new : ymax_disapp;
  float dy = (ymax - ymin) / (1. - margin_t - margin_b);
  // custom x-axis
  double x_incr = (xmax - xmin) / 50;
  double x, y;
  TLatex *t;
  double x_curr = 0;
  for (int i = 0; i < gr_trend_total->GetN(); i++) {
    gr_trend_total->GetPoint(i, x, y);
    if (x > x_curr + x_incr) {
      t = new TLatex(x, set_ymin - set_dy * margin_b / 3, Form("%i", run_numbers[i]));
      t->SetTextSize(0.025);
      t->SetTextFont(42);
      t->SetTextAlign(22);
      t->SetTextAngle(90);
      t->Draw();
      x_curr = x; 
    }
  }
  // create the title
  t = new TLatex(xmax + dx * margin_r / 2, set_ymin - set_dy * margin_b * 4/5, " Time (run numbers are shown)");
  t->SetTextSize(tsize);
  t->SetTextFont(42);
  t->SetTextAlign(32);
  t->Draw();
  // second pad
  TPad* p2 = new TPad("", "", 0., 0., 1., 1.);
  set_margins(p2, margin_t, margin_r, margin_b, margin_t);
  p2->Range(xmin-margin_l*dx, ymin-margin_b*dy, xmax+margin_r*dx, ymax+margin_t*dy);
  p2->SetFillStyle(4000); // transparent
  p2->Draw();
  p2->cd();
  gr_trend_new->GetYaxis()->SetTitle("#New/disappeared pixels");
  gr_trend_new->Draw(Form("%s", plot_opt.data()));
  gr_trend_disapp->Draw(Form("%s SAME", plot_opt.data()));
  // right axis
  TGaxis *ax = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
  ax->SetTitle("#New/disappeared noisy pixels");
  ax->SetTitleOffset(1.5);
  ax->SetTitleFont(tfont);
  ax->SetTitleSize(tsize);
  ax->SetTitleColor(kBlue);
  ax->SetLabelFont(lfont);
  ax->SetLabelSize(lsize);
  ax->SetLabelColor(kBlue);
  ax->SetLineColor(kBlue);
  ax->Draw();
  // legend
  TLegend* l = new TLegend(0.15, 0.83, 0.45, 0.97);
  l->AddEntry(gr_trend_total, Form("total + fit: #it{f}(#it{x}) = %.1e + %.1e #it{x} ", f_total->GetParameter(0), f_total->GetParameter(1)), "LP");
  l->AddEntry(gr_trend_new, "new", "P");
  l->AddEntry(gr_trend_disapp, "disappeared", "P");
  l->AddEntry((TObject*)0, Form("#runs: %i", gr_trend_total->GetN()), "");
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.03);
  l->Draw();
  c->Print(Form("%strend_ready(%.0f-%.0f)_sb(%.0f-%.0f).pdf", 
    folder.data(), std::get<2>(delays), std::get<3>(delays), std::get<0>(delays), std::get<1>(delays)));

  return;
}

void noisy_pixels_plots (long start_min, long start_max)
{
  gStyle->SetOptStat(0);
  std::vector<noise>* noise_runs = read_csv("noise_runs_input.csv", start_min, start_max);

  if(!noise_runs) {
    std::cout << "Noise run vector could not be loaded\n";
    return;
  }

  std::string folder = Form("%s_%s/",
    timestamp_to_str(std::get<0>(noise_runs->front().timestamps), "%d.%m.%y").data(),
    timestamp_to_str(std::get<0>(noise_runs->back().timestamps), "%d.%m.%y").data()
  );

  gSystem->Exec(Form("mkdir -p %s/", folder.data()));

  if (true) noisy_pix_correlation_with_delays(noise_runs, folder);
  if (true) {
    // time since last SB stop: min, max [minutes]
    // time since last GO_READY: min, max [minutes]
    noisy_pix_trend(noise_runs, "LP", {0, 1e4, 0, 1e4}, {5400, 10400}, folder);
    noisy_pix_trend(noise_runs, "LP", {0, 50, 0, 1e4}, {5400, 10400}, folder);
    noisy_pix_trend(noise_runs, "LP", {0, 50, 0, 10}, {5400, 10400}, folder);
    noisy_pix_trend(noise_runs, "LP", {50, 1e4, 0, 10}, {5400, 10400}, folder);
  }
  
  return;
}