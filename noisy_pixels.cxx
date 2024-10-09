// MFT study of the number of noisy pixels
// with respect to the last SB stop and last GO_READY
// David Grund, 2024

#include <string>
#include <vector>
#include <fstream>
#include "TH2.h"
#include "TProfile2D.h"
#include "TStyle.h"

bool is_number(const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}

struct noise
{
  int run;
  int noisy_total;
  int noisy_new;
  int noisy_disapp;
  int delay_sb;
  int delay_ready;
  noise(int r, int t, int n, int d, int sb, int re) {
    run = r; noisy_total = t; noisy_new = n; noisy_disapp = d; delay_sb = sb; delay_ready = re;
  }
};

std::vector<noise> runs;

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

bool read_csv (string fname, bool verbose = false)
{
  runs.clear();
  std::vector<std::vector<std::string>> table;
  std::vector<std::string> row;
  std::string line, item;
  // open and read the csv
  fstream f(fname,ios::in);
  if(f.is_open()) {
    while(getline(f, line)) {
      row.clear();
      stringstream str(line);
      while(getline(str, item, ',')) row.push_back(item);
      table.push_back(row);
    }
  } else {
    cout << "Cannot open " << fname << "\n";
    return false;
  }
  // print the loaded values?
  if(verbose) {
    for(int i = 0; i < (int)table.size(); i++) {
      for(int j = 0; j < (int)table[i].size(); j++) cout << table[i][j] << "\t";
      cout << "\n";
    }
  }
  // fill the list of histograms
  for(int i = 1; i < (int)table.size(); i++) {
    bool isNumbers = true;
    for (int j = 0; j < 4; j++) if (!is_number(table[i][j])) isNumbers = false;
    if (isNumbers) {
      if (std::stoi(table[i][1]) > 0) {
        noise run(std::stoi(table[i][0]), std::stoi(table[i][1]), std::stoi(table[i][2]),
          std::stoi(table[i][3]), std::stoi(table[i][4]), std::stoi(table[i][5]));
        runs.push_back(run);
      }
    }
  }
  cout << "CSV read successfully\n"
    << " #noise runs: " << runs.size() << "\n";
  return true;   
}

void set_canvas (TCanvas* c)
{
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.10);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.12);
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

void noisy_pixels ()
{
  gStyle->SetOptStat(0);
  read_csv("noise_runs_input.csv");

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
    h2_true_total->Fill(r.delay_sb, r.delay_ready, r.noisy_total);
    h2_true_count->Fill(r.delay_sb, r.delay_ready);
    h1_true_total_sb->Fill(r.delay_sb, r.noisy_total);
    h1_true_count_sb->Fill(r.delay_sb);
    h1_true_total_re->Fill(r.delay_sb, r.noisy_total);
    h1_true_count_re->Fill(r.delay_sb);
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
  c1->Print("h2_total.pdf");

  TCanvas* c2 = new TCanvas("", "", 900, 700);
  set_canvas(c2);
  format_h1(h1_bins_total_sb);
  h1_bins_total_sb->GetXaxis()->SetTitle("Time since SB stop [min]");
  for(int i = 1; i < N_sb; i++) h1_bins_total_sb->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_sb_1d[i-1], ax_sb_1d[i]));
  h1_bins_total_sb->GetXaxis()->SetBinLabel(N_sb, Form(">= %.0f", ax_sb_1d[N_sb-1]));
  h1_bins_total_sb->Draw("e0");
  c2->Print("h1_delay_sbs.pdf");

  TCanvas* c3 = new TCanvas("", "", 900, 700);
  set_canvas(c3);
  format_h1(h1_bins_total_re);
  h1_bins_total_re->GetXaxis()->SetTitle("Time since last GO_READY [min]");
  for(int i = 1; i < N_re; i++) h1_bins_total_re->GetXaxis()->SetBinLabel(i, Form("[%.0f, %.0f)", ax_re_1d[i-1], ax_re_1d[i]));
  h1_bins_total_re->GetXaxis()->SetBinLabel(N_re, Form(">= %.0f", ax_re_1d[N_re-1]));
  h1_bins_total_re->Draw("e0");
  c3->Print("h1_delay_ready.pdf");

  return;
}