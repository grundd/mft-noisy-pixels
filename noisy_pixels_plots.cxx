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

bool first_execution(true);
float xmin(-1), xmax(-1);

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

void format_graph (TGraph* g, Color_t clr, Style_t lst, Style_t mrk)
{
  g->SetLineColor(clr);
  g->SetLineStyle(lst);
  g->SetLineWidth(1);
  g->SetMarkerColor(clr);
  g->SetMarkerStyle(mrk);
  g->SetMarkerSize(0.8);
  return;
}

void noisy_pix_correlation_with_delays (std::vector<noise>* noise_runs, std::string folder = "")
{
  // axes and ranges
  const int N_2d = 5;
  float ax_5bins[N_2d+1] = {0, 1, 2, 3, 4, 5};
  float ax_sb_2d[N_2d+1] = {0, 10, 30, 60, 600, 1e5};
  float ax_re_2d[N_2d+1] = {0, 10, 30, 60, 600, 1e5};

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
  float max = h2_bins_total->GetMaximum();
  float min = h2_bins_total->GetMinimum(0);
  float margin = (max - min) * 0.07;
  h2_bins_total->GetZaxis()->SetRangeUser(min - margin, max + margin);
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
  std::tuple<float, float, float, float> force_ranges = {-1, -1, -1, -1},
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

  format_graph(gr_trend_total, kBlack, 1, kOpenCircle);
  format_graph(gr_trend_new, kBlue, 1, kOpenSquare);
  format_graph(gr_trend_disapp, kBlue, 2, kFullTriangleUp);

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

  // x-axis
  if(first_execution) {
    xmin = gr_trend_total->GetXaxis()->GetXmin();
    xmax = gr_trend_total->GetXaxis()->GetXmax();
  }
  float dx = (xmax - xmin) / (1. - margin_l - margin_r);
  gr_trend_total->GetXaxis()->SetLabelSize(0);
  gr_trend_total->GetXaxis()->SetTitleSize(0);
  gr_trend_total->GetXaxis()->SetTickSize(0);
  TAxis *x_ax = gr_trend_total->GetXaxis();
  x_ax->SetLimits(xmin, xmax);

  // y-axis
  float ymin = gr_trend_total->GetHistogram()->GetMinimum() * 0.75;
  float ymax = gr_trend_total->GetHistogram()->GetMaximum() * 1.05;
  if (std::get<0>(force_ranges) > 0) ymin = std::get<0>(force_ranges);
  if (std::get<1>(force_ranges) > 0) ymax = std::get<1>(force_ranges);
  float dy = (ymax - ymin) / (1. - margin_t - margin_b);
  gr_trend_total->GetYaxis()->SetTitle("#Noisy pixels");
  gr_trend_total->GetYaxis()->SetTitleOffset(1.65);
  gr_trend_total->GetYaxis()->SetRangeUser(ymin, ymax);
  
  gr_trend_total->Draw(Form("A%s", plot_opt.data()));

  // get font sizes and styles from gr_trend_total
  Style_t tfont = gr_trend_total->GetHistogram()->GetYaxis()->GetTitleFont();
  Style_t lfont = gr_trend_total->GetHistogram()->GetYaxis()->GetLabelFont();
  float tsize = gr_trend_total->GetHistogram()->GetYaxis()->GetTitleSize();
  float lsize = gr_trend_total->GetHistogram()->GetYaxis()->GetLabelSize();

  // custom x-axis
  double x_incr = (xmax - xmin) / 50;
  double x, y;
  TLatex *t;
  double x_curr = 0;
  for (int i = 0; i < gr_trend_total->GetN(); i++) {
    gr_trend_total->GetPoint(i, x, y);
    if (x > x_curr + x_incr) {
      // label
      t = new TLatex(x, ymin - dy * margin_b / 3, Form("%i", run_numbers[i]));
      t->SetTextSize(0.025);
      t->SetTextFont(42);
      t->SetTextAlign(22);
      t->SetTextAngle(90);
      t->Draw();
      // tick
      TLine* l = new TLine(x, ymin, x, ymin + (ymax - ymin) * 0.015);
      l->Draw();

      x_curr = x;  
    }
  }

  // create the title
  t = new TLatex(xmax + dx * margin_r / 2, ymin - dy * margin_b * 4/5, " Time (run numbers are shown)");
  t->SetTextSize(tsize);
  t->SetTextFont(42);
  t->SetTextAlign(32);
  t->Draw();

  // set ranges for gr_trend_new & gr_trend_disapp
  float ymin2 = 0; 
  float ymax_new = gr_trend_new->GetHistogram()->GetMaximum() * 1.8;
  float ymax_disapp = gr_trend_disapp->GetHistogram()->GetMaximum() * 1.8;
  float ymax2 = ymax_new > ymax_disapp ? ymax_new : ymax_disapp;
  if (std::get<2>(force_ranges) > 0) ymin2 = std::get<2>(force_ranges);
  if (std::get<3>(force_ranges) > 0) ymax2 = std::get<3>(force_ranges);  
  float dy2 = (ymax2 - ymin2) / (1. - margin_t - margin_b);
  
  // second pad
  TPad* p2 = new TPad("", "", 0., 0., 1., 1.);
  set_margins(p2, margin_t, margin_r, margin_b, margin_t);
  p2->Range(xmin-margin_l*dx, ymin2-margin_b*dy2, xmax+margin_r*dx, ymax2+margin_t*dy2);
  p2->SetFillStyle(4000); // transparent
  p2->Draw();
  p2->cd();
  gr_trend_new->GetYaxis()->SetTitle("#New/disappeared pixels");
  gr_trend_new->Draw(Form("%s", plot_opt.data()));
  gr_trend_disapp->Draw(Form("%s SAME", plot_opt.data()));
  // right axis
  TGaxis *ax = new TGaxis(xmax, ymin2, xmax, ymax2, ymin2, ymax2, 510, "+L");
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
  l->AddEntry(gr_trend_total, Form("total + fit: #it{f}(#it{x}) = %.2e + %.2e #it{x} ", 
    f_total->GetParameter(0), f_total->GetParameter(1)), plot_opt.data());
  l->AddEntry(gr_trend_new, "new", plot_opt.data());
  l->AddEntry(gr_trend_disapp, "disappeared", plot_opt.data());
  l->AddEntry((TObject*)0, Form("#runs: %i", gr_trend_total->GetN()), "");
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextSize(0.03);
  l->Draw();
  c->Print(Form("%strend_ready(%.0f-%.0f)_sb(%.0f-%.0f).pdf", 
    folder.data(), std::get<2>(delays), std::get<3>(delays), std::get<0>(delays), std::get<1>(delays)));

  first_execution = false;
  return;
}

void noisy_pixels_plots (long start_min, long start_max)
{
  gStyle->SetOptStat(0);
  std::vector<noise>* noise_runs = read_csv("input.csv", start_min, start_max);

  if(!noise_runs) {
    std::cout << "Noise run vector could not be loaded\n";
    return;
  }

  std::string folder = Form("%s_%s/",
    timestamp_to_str(std::get<0>(noise_runs->front().timestamps), "%y.%m.%d").data(),
    timestamp_to_str(std::get<0>(noise_runs->back().timestamps), "%y.%m.%d").data()
  );

  gSystem->Exec(Form("mkdir -p %s/", folder.data()));

  if (true) noisy_pix_correlation_with_delays(noise_runs, folder);
  if (true) {
    // time since last SB stop: min, max [minutes]
    // time since last GO_READY: min, max [minutes]
    noisy_pix_trend(noise_runs, "LP", {0, 1e5, 0, 1e4}, {5400, 10600, 0, 4200}, folder); // all runs
    noisy_pix_trend(noise_runs, "LP", {0, 30, 30, 1e4}, {5400, 10600, 0, 4200}, folder); // "standard" noise runs only
    noisy_pix_trend(noise_runs, "LP", {60, 1e5, 0, 30}, {5400, 10600, 0, 4200}, folder); // group 2
    noisy_pix_trend(noise_runs, "LP", {60, 1e5, 0, 5}, {5400, 10600, 0, 4200}, folder); // group 2 extreme
  }
  
  return;
}