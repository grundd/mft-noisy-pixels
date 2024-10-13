// MFT study of the number of noisy pixels
// with respect to the last SB stop and last GO_READY
// David Grund, 2024

// cpp headers
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include <fstream>
// root headers
#include "TSystem.h"
#include "TFile.h"
// o2 headers
#include "CCDB/CcdbApi.h"
#include "DataFormatsITSMFT/NoiseMap.h"

o2::ccdb::CcdbApi api;

std::tuple<int, long, long> read_noise_maps (long ts, bool rewrite, bool verbose = false)
{
  gSystem->Exec("mkdir -p noise_maps/");

  std::map<std::string, std::string> filter;
  auto calib = api.retrieveFromTFileAny<o2::itsmft::NoiseMap>("MFT/Calib/NoiseMap/", filter, ts);

  std::map<std::string, std::string> headers = api.retrieveHeaders("MFT/Calib/NoiseMap/", filter, ts);
  if (verbose)
  {
    std::map<std::string, std::string>::iterator it;
    for (it = headers.begin(); it != headers.end(); it++) std::cout << it->first << "\t" << it->second << "\n";
  }

  int run_no = std::stoi(headers["runNumber"]);
  long val_from = std::stol(headers["Valid-From"]);
  long val_until = std::stol(headers["Valid-Until"]);

  std::string fname = Form("noise_maps/%i.root", run_no);
  bool exists = !gSystem->AccessPathName(fname.data());
  if(exists && !rewrite) 
  {
    std::cout << "Noise map for " << run_no << " already read\n";
  } 
  else 
  {
    std::cout << "Reading noise map for " << run_no << "\n";
    std::vector<std::vector<int>> noisy_pixs; // chipID, row, col, noise
    int noisy_pix_count = 0;

    const int N_pix = 936; // 936
    for (int chipID = 0; chipID < N_pix; chipID++) {
      if ((chipID + 1) % 50 == 0) std::cout << " " << chipID+1 << " chips read\n";
      for (int row = 0; row < 512; row++) {
        for (int col = 0; col < 1024; col++) {
          int noise = calib->getNoiseLevel(chipID, row, col);
          if (!noise) {
            noise = -1;
          } else {
            noisy_pixs.push_back({chipID, row, col, noise});
            noisy_pix_count++;
          }
        }
      }
    }

    std::sort(noisy_pixs.begin(), noisy_pixs.end(), 
      [](auto const& lhs, auto const& rhs) { return lhs[0] < rhs[0]; }
    );

    std::cout << " Done: " << noisy_pix_count << " noisy pixels found\n";

    TFile* fout = TFile::Open(fname.data(), "recreate");
    fout->WriteObject(&noisy_pixs, "noisy_pixs");
    fout->Close();
  }
  return {run_no, val_from, val_until};
}

void compare_noise_maps (long ts_first, long ts_last)
{
  long ts_curr = ts_first;

  std::vector<std::vector<int>> noise_run_stats;

  while (ts_curr < ts_last)
  {
    // current noise map
    std::tuple<int, long, long> info_curr = read_noise_maps(ts_curr, false);
    int run_curr = std::get<0>(info_curr);
    long val_from_curr = std::get<1>(info_curr);
    long val_until_curr = std::get<2>(info_curr);

    // previous noise run
    std::tuple<int, long, long> info_prev = read_noise_maps(val_from_curr-1, false);
    int run_prev = std::get<0>(info_prev);
    long val_from_prev = std::get<1>(info_prev);
    long val_until_prev = std::get<2>(info_prev);

    std::string fname_curr = Form("noise_maps/%i.root", run_curr);
    std::string fname_prev = Form("noise_maps/%i.root", run_prev);
    bool exists_curr = !gSystem->AccessPathName(fname_curr.data());
    bool exists_prev = !gSystem->AccessPathName(fname_prev.data());
    if(exists_curr && exists_prev) 
    {
      std::cout << "Run " << run_curr << "\n";

      TFile *f_curr = TFile::Open(fname_curr.data(), "read");
      std::vector<std::vector<int>>* pixs_curr;
      f_curr->GetObject("noisy_pixs", pixs_curr); 

      TFile *f_prev = TFile::Open(fname_prev.data(), "read");
      std::vector<std::vector<int>>* pixs_prev;
      f_prev->GetObject("noisy_pixs", pixs_prev); 

      std::vector<std::vector<int>>::iterator it_curr;
      std::vector<std::vector<int>>::iterator it_prev;

      int common_pixs = 0;
      for (it_curr = pixs_curr->begin(); it_curr != pixs_curr->end(); it_curr++) 
      {
        int chip_curr = it_curr->at(0);
        for(it_prev = pixs_prev->begin(); it_prev != pixs_prev->end(); it_prev++) 
        {
          int chip_prev = it_prev->at(0);
          if (chip_prev > chip_curr) continue;
          
          int row_curr = it_curr->at(1);
          int col_curr = it_curr->at(2);
          int row_prev = it_prev->at(1);
          int col_prev = it_prev->at(2);

          if ((chip_curr == chip_prev) && (row_curr == row_prev) && (col_curr == col_prev)) common_pixs++;          
        }
      }
      int noisy_total = pixs_curr->size();
      int noisy_new = pixs_curr->size() - common_pixs;
      int noisy_disapp = pixs_prev->size() - common_pixs;
      std::cout << " Total: " << noisy_total << "\n"
                << " New: " << noisy_new << "\n"
                << " Disappeared: " << noisy_disapp << "\n";
      noise_run_stats.push_back({run_curr, noisy_total, noisy_new, noisy_disapp});
    }
    else 
    {
      std::cout << "Noise map for " << run_curr << " or " << run_prev << " not found\n";
      break;
    }

    ts_curr = val_until_curr+1;
    std::cout << " New timestamp: " << ts_curr << "\n";
  }

  // create and save the csv
  ofstream csv("noisy_pixels.csv");
  csv << "run,total,new,disapp\n";
  for (auto noise : noise_run_stats)
    csv << noise[0] << "," << noise[1] << "," << noise[2] << "," << noise[3] << "\n";
  csv.close();

  return;
}

void noisy_pixels_count ()
{
  api.init("http://alice-ccdb.cern.ch");

  compare_noise_maps(1714531157487, 1760266579564);
  return;
}