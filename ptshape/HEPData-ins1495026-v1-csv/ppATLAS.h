#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TMath.h>

namespace ppRef
{
  class ppATLAS
  {
  public:
    ppATLAS(std::string dirname) : fdir(dirname) { init(); }
    ~ppATLAS() { ; }
    const std::vector<std::string> types = {"promptPsi", "nonpromptPsi", "promptX", "nonpromptX"};
    const std::vector<std::string> err = {"stat", "syst"};
    std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> gg;
    std::map<std::string, std::map<std::string, TH1F*>> hh;
    TH2F* hempty;
    float getycut() { return ycut; }

  private:
    const float ycut = 0.75;
    std::string fdir;
    void init();
  };
}

void ppRef::ppATLAS::init()
{
  for(auto& tt : types)
    {
      std::ifstream inf(fdir+"/"+tt+".dat");
      std::string line;
      std::vector<float> xx, yy, xxel, xxeh, yyel_stat, yyeh_stat, yyel_syst, yyeh_syst;
      std::vector<float> xbin, yhist, yehist_stat, yehist_syst;
      while(std::getline(inf, line))
        {
          if(line.empty() || line.length() <= 1) break;
          if(line.front() == '#') continue;

          float x, y, xel, xeh, yel_stat, yeh_stat, yel_syst, yeh_syst;
          std::istringstream ss(line);
          ss >> x >> xeh >> xel >> y >> yeh_stat >> yel_stat >> yeh_syst >> yel_syst;
          y *= 1.e-3; yel_stat*=1.e-3; yeh_stat*=1.e-3; yel_syst*=1.e-3; yeh_syst*=1.e-3;

          if(xbin.empty()) xbin.push_back(xel);
          xbin.push_back(xeh);
          yhist.push_back(y);
          yehist_stat.push_back(TMath::Abs(yeh_stat));
          yehist_syst.push_back(TMath::Abs(yeh_syst));

          xx.push_back(x);
          yy.push_back(y);
          xxel.push_back(x-xel);
          xxeh.push_back(xeh-x);
          yyel_stat.push_back(TMath::Abs(yel_stat));
          yyeh_stat.push_back(TMath::Abs(yeh_stat));
          yyel_syst.push_back(TMath::Abs(yel_syst));
          yyeh_syst.push_back(TMath::Abs(yeh_syst));
        }
      gg[tt]["stat"] = new TGraphAsymmErrors(xx.size(), xx.data(), yy.data(), xxel.data(), xxeh.data(), yyel_stat.data(), yyeh_stat.data());
      gg[tt]["stat"]->SetName(Form("grATLAS_%s_%s", tt.c_str(), "stat"));
      gg[tt]["syst"] = new TGraphAsymmErrors(xx.size(), xx.data(), yy.data(), xxel.data(), xxeh.data(), yyel_syst.data(), yyeh_syst.data());
      gg[tt]["syst"]->SetName(Form("grATLAS_%s_%s", tt.c_str(), "syst"));
      
      hh[tt]["stat"] = new TH1F(Form("hATLAS_%s_%s", tt.c_str(), "stat"), "", xbin.size()-1, xbin.data());
      hh[tt]["syst"] = new TH1F(Form("hATLAS_%s_%s", tt.c_str(), "syst"), "", xbin.size()-1, xbin.data());
      for(int i=0; i<xbin.size()-1; i++)
        {
          hh[tt]["stat"]->SetBinContent(i+1, yhist[i]);
          hh[tt]["stat"]->SetBinError(i+1, yehist_stat[i]);
          hh[tt]["syst"]->SetBinContent(i+1, yhist[i]);
          hh[tt]["syst"]->SetBinError(i+1, yehist_syst[i]);
        }
      inf.close();
    }
  hempty = new TH2F("hempty_ppATLAS", ";p_{T} (GeV/c);Br #times d#sigma/dp_{T}", 10, 10, 70, 10, 1.e-7, 1.);
}
