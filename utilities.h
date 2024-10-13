// MFT study of the number of noisy pixels
// with respect to the last SB stop and last GO_READY
// David Grund, 2024

// cpp headers
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>

// structure to store noise run information
struct noise
{
  int run;
  std::tuple<int, int, int> noisy_pixels; // total, new, disapp
  std::tuple<long, long, long> timestamps; // trg_start, delay_sb, delay_ready
  noise(int r, 
    std::tuple<int, int, int> px, 
    std::tuple<long, long, long> ts
  ) {
    run = r;
    noisy_pixels = px;
    timestamps = ts;
  }
};

// is the string a positive integer?
bool is_positive_int (const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}

// convert date string with specified format to long (unix timestamp)
long datestring_to_long (const std::string& s, std::string format, bool verbose = false)
{
  long ts_long;
  if(s.length() < 2) { // empty string or end-of-line character
    ts_long = 0;
  } else {
    std::tm t;
    std::istringstream ss(s);
    ss >> std::get_time(&t, format.data());
    if (ss.fail()) {
      throw std::runtime_error(Form("Failed to convert %s to time", s.data()));
    }
    ts_long = (long)std::mktime(&t);
    if(verbose) std::cout << s << " -> " << ts_long << "\n";
  }
  return ts_long;
}

std::string timestamp_to_str (long ts, string format)
{
  std::time_t ts_as_time_t = ts; // convert from long to time_t
  auto ts_as_tm = std::localtime(&ts_as_time_t);
  char buff[80];
  std::strftime(buff, sizeof(buff), format.data(), ts_as_tm);
  std::string date(buff);
  return date;
}

std::vector<noise>* read_csv (string fname, long start_min = 0, long start_max = 2e9, bool verbose = false)
{
  std::vector<noise>* runs = new std::vector<noise>;
  std::vector<std::vector<std::string>> csv;
  std::vector<std::string> row;
  std::string line, item;

  // open and read the csv
  fstream f(fname,ios::in);
  if(f.is_open()) 
  {
		while(getline(f,line)) 
    {
		  row.clear();
		  stringstream str(line);
      // go until the first double quote
      while(getline(str, item, '"')) 
      {
        // collect all items separated by commas
        stringstream up_to_quote(item);
        while(getline(up_to_quote, item, ',')) row.push_back(item);
        // get the full surrounded by double quotes
        if(getline(str, item, '"')) row.push_back(item);
        // move behind the comma after the closing double quote
        getline(str, item, ',');
      }
      csv.push_back(row);
		}
	} 
  else 
  {
    std::cout << "Cannot open " << fname << "\n";
    return NULL;
  }
  // print the loaded values?
  if(verbose) 
  {
    for(int i = 0; i < (int)csv.size(); i++) {
      for(int j = 0; j < (int)csv[i].size(); j++) std::cout << csv[i][j] << "\t";
      std::cout << "\n";
    }
  }
  // fill the vector of structures
  for(int i = 1; i < (int)csv.size(); i++) 
  {
    // run number and noisy_pixels (tuple) need to be positive integers
    bool positive_int = true;
    for (int j = 0; j < 4; j++) if (!is_positive_int(csv[i][j])) positive_int = false;
    if (positive_int) 
    {
      // #total noisy pixels > 0
      if (std::stoi(csv[i][1]) > 0) 
      {
        std::tuple<int, int, int> noisy_pixels = {
          std::stoi(csv[i][1]), std::stoi(csv[i][2]), std::stoi(csv[i][3])
        };
        std::tuple<long, long, long> timestamps = {
          datestring_to_long(csv[i][4], "%d/%m/%Y, %H:%M:%S"), 
          datestring_to_long(csv[i][5], "%d/%m/%Y, %H:%M:%S"),
          datestring_to_long(csv[i][6], "%d/%m/%Y, %H:%M:%S")
        };
        noise run(std::stoi(csv[i][0]), noisy_pixels, timestamps);
        if ((std::get<0>(timestamps) > start_min) && (std::get<0>(timestamps) < start_max)) {
          runs->push_back(run);
        }
      }
    }
  }
  std::cout << "CSV read successfully\n"
    << " #noise runs: " << runs->size() << "\n";
  return runs; 
}