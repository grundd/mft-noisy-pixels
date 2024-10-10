// MFT study of the number of noisy pixels
// with respect to the last SB stop and last GO_READY
// David Grund, 2024

// cpp headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
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

// is the string positive integer?
bool is_pos_int (const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}

// convert datestring with the specified format to long (unix timestamp)
long datestring_to_long (const std::string& s, bool verbose = false)
{
  long ts_long;
  if(s.length() < 2) { // empty string or end-of-line characters
    ts_long = 0;
  } else {
    std::tm t;
    std::istringstream ss(s);
    ss >> std::get_time(&t, "%d/%m/%Y, %H:%M:%S");
    if (ss.fail()) {
      throw std::runtime_error(Form("Failed to convert %s to time", s.data()));
    }
    ts_long = (long)std::mktime(&t);
    if(verbose) std::cout << s << " -> " << ts_long << "\n";
  }
  return ts_long;
}

std::string timestamp_to_mmddyyyy (long ts)
{
  std::time_t ts_as_time_t = ts; // convert from long to time_t
  auto ts_as_tm = std::localtime(&ts_as_time_t);
  char buff[80];
  std::strftime(buff, sizeof(buff), "%d-%m-%Y", ts_as_tm);
  std::string date(buff);
  return date;
}

// structure to store noise run information
struct noise
{
  int run;
  int noisy_total;
  int noisy_new;
  int noisy_disapp;
  long ts_trg_start;
  long ts_sb_stop;
  long ts_go_ready;
  noise(int r, int t, int n, int d, long start, long sbs, long ready) {
    run = r; 
    noisy_total = t; 
    noisy_new = n; 
    noisy_disapp = d; 
    ts_trg_start = start;
    ts_sb_stop = sbs; 
    ts_go_ready = ready;
  }
};

std::vector<noise> runs;
std::string folder_name;

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

bool read_csv (string fname, long start_min = 0, long start_max = 2e9, bool verbose = false)
{
  runs.clear();
  std::vector<std::vector<std::string>> table;
  std::vector<std::string> row;
  std::string line, item;
  // open and read the csv
  fstream f(fname,ios::in);
  if(f.is_open()) {
		while(getline(f,line)) {
		  row.clear();
      // https://stackoverflow.com/questions/17859992/reading-from-a-csv-text-file-with-quotes-in-c
		  stringstream str(line);
      // go until the first double quote
      while(getline(str, item, '"')) {
        // collect all items separated by commas
        stringstream up_to_quote(item);
        while(getline(up_to_quote, item, ',')) row.push_back(item);
        // get the full surrounded by double quotes
        if(getline(str, item, '"')) row.push_back(item);
        // move behind the comma after the closing double quote
        getline(str, item, ',');
      }
      table.push_back(row);
		}
	} else {
    std::cout << "Cannot open " << fname << "\n";
    return false;
  }
  // print the loaded values?
  if(verbose) {
    for(int i = 0; i < (int)table.size(); i++) {
      for(int j = 0; j < (int)table[i].size(); j++) std::cout << table[i][j] << "\t";
      std::cout << "\n";
    }
  }
  // fill the vector of structures
  for(int i = 1; i < (int)table.size(); i++) {
    bool pos_integers = true;
    for (int j = 0; j < 4; j++) if (!is_pos_int(table[i][j])) pos_integers = false;
    if (pos_integers) {
      if (std::stoi(table[i][1]) > 0) {
        long ts_trg_start = datestring_to_long(table[i][4]);
        long ts_sb_stop = datestring_to_long(table[i][5]);
        long ts_go_ready = datestring_to_long(table[i][6]);
        if (ts_trg_start > start_min && ts_trg_start < start_max) {
          noise run(std::stoi(table[i][0]), std::stoi(table[i][1]), std::stoi(table[i][2]), std::stoi(table[i][3]), 
            ts_trg_start, ts_sb_stop, ts_go_ready);
          runs.push_back(run);
        }
      }
    }
  }
  std::cout << "CSV read successfully\n"
    << " #noise runs: " << runs.size() << "\n";
  return true; 
}

float margin_top = 0.07;
float margin_right = 0.10;
float margin_bottom = 0.10;
float margin_left = 0.12;

template<typename T>
void set_canvas_for_h1 (T* c)
{
  c->SetTopMargin(margin_top);
  c->SetBottomMargin(margin_bottom);
  c->SetRightMargin(margin_right);
  c->SetLeftMargin(margin_left);
  return;
}

template<typename TH>
void format_h2 (TH* h)
{
  for(int i = 1; i < N_2d; i++) {
    h->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_sb_2d[i-1], ax_sb_2d[i]));
    h->GetYaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_re_2d[i-1], ax_re_2d[i]));
  }
  h->GetXaxis()->SetBinLabel(N_2d, Form(">= %.0f", ax_sb_2d[N_2d-1]));
  h->GetYaxis()->SetBinLabel(N_2d, Form(">= %.0f", ax_re_2d[N_2d-1]));
  // x-axis
  h->GetXaxis()->SetTitle("Time since SB stop [min]");
  h->GetXaxis()->SetTitleSize(0.035);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetLabelSize(0.04);
  // y-axis
  h->GetYaxis()->SetTitle("Time since last GO_READY [min]");
  h->GetYaxis()->SetTitleSize(0.035);
  h->GetYaxis()->SetTitleOffset(2.2);
  h->GetYaxis()->SetLabelSize(0.04);
  // z-axis
  h->GetZaxis()->SetLabelSize(0.03);
  return;
}

template<typename TH>
void format_h1 (TH* h)
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
  h->GetYaxis()->SetRangeUser(7800, 9000);
  return;
}

void noisy_pix_correlation_with_delays ()
{
  TH2F* h2_true_total = new TH2F("", "", N_2d, &ax_sb_2d[0], N_2d, &ax_re_2d[0]);
  TH2I* h2_true_count = new TH2I("", "", N_2d, &ax_sb_2d[0], N_2d, &ax_re_2d[0]);

  TH2F* h2_bins_total = new TH2F("", "#Noisy pixels (#color[4]{#runs})", N_2d, &ax_5bins[0], N_2d, &ax_5bins[0]);
  TH2I* h2_bins_count = new TH2I("", "#Runs", N_2d, &ax_5bins[0], N_2d, &ax_5bins[0]);

  TH1F* h1_true_total_sb = new TH1F("", "", N_sb, &ax_sb_1d[0]);
  TH1I* h1_true_count_sb = new TH1I("", "", N_sb, &ax_sb_1d[0]);
  TH1F* h1_bins_total_sb = new TH1F("", "#Noisy pixels", N_sb, &ax_8bins[0]);

  TH1F* h1_true_total_re = new TH1F("", "", N_re, &ax_re_1d[0]);
  TH1I* h1_true_count_re = new TH1I("", "", N_re, &ax_re_1d[0]);
  TH1F* h1_bins_total_re = new TH1F("", "#Noisy pixels", N_re, &ax_7bins[0]);

  for (auto r : runs) {
    float delay_sb = (r.ts_trg_start - r.ts_sb_stop) / 60; // in minutes
    float delay_re = (r.ts_trg_start - r.ts_go_ready) / 60; // in minutes
    h2_true_total->Fill(delay_sb, delay_re, r.noisy_total);
    h2_true_count->Fill(delay_sb, delay_re);
    h1_true_total_sb->Fill(delay_sb, r.noisy_total);
    h1_true_count_sb->Fill(delay_sb);
    h1_true_total_re->Fill(delay_re, r.noisy_total);
    h1_true_count_re->Fill(delay_re);
  }
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

  TCanvas* c1 = new TCanvas("", "", 900, 700);
  c1->SetTopMargin(0.07);
  c1->SetBottomMargin(0.10);
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(0.145);
  format_h2(h2_bins_total);
  h2_bins_total->GetZaxis()->SetRangeUser(7500, 9000);
  h2_bins_total->SetBarOffset(+0.08);
  h2_bins_total->Draw("colz text");
  h2_bins_count->SetBarOffset(-0.08);
  h2_bins_count->SetMarkerColor(kBlue);
  h2_bins_count->Draw("same text");
  c1->Print(Form("%s/corr_2d.pdf", folder_name.data()));

  TCanvas* c2 = new TCanvas("", "", 900, 700);
  set_canvas_for_h1(c2);
  format_h1(h1_bins_total_sb);
  h1_bins_total_sb->GetXaxis()->SetTitle("Time since SB stop [min]");
  for(int i = 1; i < N_sb; i++) h1_bins_total_sb->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_sb_1d[i-1], ax_sb_1d[i]));
  h1_bins_total_sb->GetXaxis()->SetBinLabel(N_sb, Form(">= %.0f", ax_sb_1d[N_sb-1]));
  h1_bins_total_sb->Draw("e0");
  c2->Print(Form("%s/corr_SB_stop.pdf", folder_name.data()));

  TCanvas* c3 = new TCanvas("", "", 900, 700);
  set_canvas_for_h1(c3);
  format_h1(h1_bins_total_re);
  h1_bins_total_re->GetXaxis()->SetTitle("Time since last GO_READY [min]");
  for(int i = 1; i < N_re; i++) h1_bins_total_re->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_re_1d[i-1], ax_re_1d[i]));
  h1_bins_total_re->GetXaxis()->SetBinLabel(N_re, Form(">= %.0f", ax_re_1d[N_re-1]));
  h1_bins_total_re->Draw("e0");
  c3->Print(Form("%s/corr_last_GO_READY.pdf", folder_name.data()));

  return;
}

void noisy_pix_trend (float delay_sb_min = 0, float delay_sb_max = 1e4, float delay_re_min = 0, float delay_re_max = 1e4,
  float force_ymin = -1, float force_ymax = -1)
{
  TGraph* gr_trend_total = new TGraph();
  TGraph* gr_trend_new = new TGraph();
  TGraph* gr_trend_disapp = new TGraph();

  for (auto r : runs) {
    float delay_sb = (r.ts_trg_start - r.ts_sb_stop) / 60; // in minutes
    float delay_re = (r.ts_trg_start - r.ts_go_ready) / 60; // in minutes
    if(delay_sb >= delay_sb_min && delay_sb < delay_sb_max && delay_re >= delay_re_min && delay_re < delay_re_max) {
      gr_trend_total->AddPoint(r.run, r.noisy_total);
      gr_trend_new->AddPoint(r.run, r.noisy_new);
      gr_trend_disapp->AddPoint(r.run, r.noisy_disapp);
    }
  }

  gr_trend_total->SetLineColor(kBlack);
  gr_trend_total->SetLineStyle(1);
  gr_trend_total->SetLineWidth(2);
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

  TCanvas* c = new TCanvas("", "", 900, 700);
  c->cd();
  // first pad
  TPad* p1 = new TPad("", "", 0., 0., 1., 1.);
  set_canvas_for_h1(p1);
  p1->Draw();
  p1->cd();
  gr_trend_total->GetXaxis()->SetTitle("Run number");
  gr_trend_total->GetXaxis()->SetTitleOffset(1.1);
  gr_trend_total->GetYaxis()->SetTitle("#Noisy pixels");
  gr_trend_total->GetYaxis()->SetTitleOffset(1.65);
  gr_trend_total->GetHistogram()->SetNdivisions(510, "X");
  gr_trend_total->GetXaxis()->SetMaxDigits(6);
  float set_ymin = gr_trend_total->GetHistogram()->GetMinimum();
  float set_ymax = gr_trend_total->GetHistogram()->GetMaximum();
  if (force_ymin > 0) set_ymin = force_ymin;
  if (force_ymax > 0) set_ymax = force_ymax;
  gr_trend_total->GetYaxis()->SetRangeUser(set_ymin * 0.75, set_ymax);
  gr_trend_total->Draw("ALP");
  // font sizes and styles
  Style_t tfont = gr_trend_total->GetHistogram()->GetYaxis()->GetTitleFont();
  Style_t lfont = gr_trend_total->GetHistogram()->GetYaxis()->GetLabelFont();
  float tsize = gr_trend_total->GetHistogram()->GetYaxis()->GetTitleSize();
  float lsize = gr_trend_total->GetHistogram()->GetYaxis()->GetLabelSize();
  // ranges
  float xmin = gr_trend_total->GetXaxis()->GetXmin();
  float xmax = gr_trend_total->GetXaxis()->GetXmax();
  float dx = (xmax - xmin) / (1. - margin_left - margin_right);
  float ymin = 0; 
  float ymax_new = gr_trend_new->GetHistogram()->GetMaximum()*1.8;
  float ymax_disapp = gr_trend_disapp->GetHistogram()->GetMaximum()*1.8;
  float ymax = ymax_new > ymax_disapp ? ymax_new : ymax_disapp;
  float dy = (ymax - ymin) / (1. - margin_top - margin_bottom);
  // second pad
  TPad* p2 = new TPad("", "", 0., 0., 1., 1.);
  set_canvas_for_h1(p2);
  p2->Range(xmin-margin_left*dx, ymin-margin_bottom*dy, xmax+margin_right*dx, ymax+margin_top*dy);
  p2->SetFillStyle(4000); // transparent
  p2->Draw();
  p2->cd();
  gr_trend_new->GetYaxis()->SetTitle("#New/disappeared pixels");
  gr_trend_new->Draw("LP");
  gr_trend_disapp->Draw("LP same");
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
  TLegend* l = new TLegend(0.15, 0.82, 0.45, 0.92);
  l->AddEntry(gr_trend_total, "total", "LP");
  l->AddEntry(gr_trend_new, "new", "LP");
  l->AddEntry(gr_trend_disapp, "disappeared", "LP");
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->Draw();
  c->Print(Form("%s/trend_SB(%.0f-%.0f)_READY(%.0f-%.0f).pdf", 
    folder_name.data(), delay_sb_min, delay_sb_max, delay_re_min, delay_re_max));

  return;
}

void noisy_pixels (long start_min, long start_max)
{
  gStyle->SetOptStat(0);
  read_csv("noise_runs_input.csv", start_min, start_max);
  folder_name = timestamp_to_mmddyyyy(runs.front().ts_trg_start) + "_" + timestamp_to_mmddyyyy(runs.back().ts_trg_start);
  gSystem->Exec(Form("mkdir -p %s/", folder_name.data()));

  if (true) noisy_pix_correlation_with_delays();
  if (true) {
    // time since last SB stop: min, max [minutes]
    // time since last GO_READY: min, max [minutes]
    noisy_pix_trend(0, 1e4, 0, 1e4, 5400, 10800);
    noisy_pix_trend(0, 10, 0, 1e4, 5400, 10800);
    noisy_pix_trend(60, 1e4, 0, 10, 5400, 10800);
  }
  
  return;
}