#ifndef __PPREF_GETDATA__
#define __PPREF_GETDATA__

#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>

namespace ppref
{
  class getdata
  {
  public:
    getdata(std::string inputfile, bool nonprompt = false, float scale = 1.);
    getdata(std::map<std::string, std::vector<float>>& v, int n, std::string name);
    std::vector<float>& operator[](std::string var) { return val_[var]; }
    int n() { return n_; }
    float binwidth(int j) { if(j < n_) { return val_["ptmax"][j]-val_["ptmin"][j]; } 
      else { return 0; }}
    void print();
  private:
    int n_;
    std::string name_;
    std::map<std::string, std::vector<float>> val_;
    std::ifstream g_;
    float scale_;
    void calcrel();
  };
}

ppref::getdata::getdata(std::string inputfile, bool nonprompt, float scale) : name_(inputfile), scale_(scale)
{
  std::cout << inputfile << std::endl;
  g_.open(inputfile.c_str());
  n_ = 0;
  while(true)
    {
      float ptmin, ptmax, center, stat, syst;
      g_ >> ptmin;
      if(g_.eof()) break;
      g_ >> ptmax >> center >> stat >> syst;
      if(nonprompt) center = 1 - center;

      val_["ptmin"].push_back(ptmin);
      val_["ptmax"].push_back(ptmax);
      val_["center"].push_back(center * scale);
      val_["stat"].push_back(stat * scale);
      val_["syst"].push_back(syst * scale);

      n_++;
    }
  calcrel();
}

ppref::getdata::getdata(std::map<std::string, std::vector<float>>& v, int n, std::string name) : n_(n), name_(name), val_(v), scale_(1.)
{
  calcrel();
}

void ppref::getdata::calcrel()
{
  for(int i=0; i<n_; i++)
    {
      val_["stat_rel"].push_back(val_["stat"][i]/val_["center"][i]);
      val_["syst_rel"].push_back(val_["syst"][i]/val_["center"][i]);
    }
}

void ppref::getdata::print()
{
  int wid = 10;
  int totalwid = wid*5 + 1;
  std::cout<<name_<<std::endl;
  std::cout << std::string(totalwid, '-') << std::endl;
  std::cout << std::left
            << std::setw(wid) << "| ptmin"
            << std::setw(wid) << "| ptmax"
            << std::setw(wid) << "| center"
            << std::setw(wid) << "| stat"
            << std::setw(wid) << "| syst"
            << "|" <<std::endl;
  std::cout << std::string(totalwid, '-') << std::endl;
  for(int i=0; i<n_; i++)
    {
      std::cout << std::left
                << std::setw(wid) << Form("| %.1f", val_["ptmin"][i])
                << std::setw(wid) << Form("| %.1f", val_["ptmax"][i])
                << std::setw(wid) << Form("| %.4f", val_["center"][i])
                << std::setw(wid) << Form("| %.4f", val_["stat"][i])
                << std::setw(wid) << Form("| %.4f", val_["syst"][i])
                << "|" <<std::endl;
      std::cout << std::string(totalwid, '-') << std::endl;
    }
}

#endif
