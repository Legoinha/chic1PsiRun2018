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
#include <TF1.h>

#include "xjjrootuti.h"

namespace ppref
{
  class ppATLAS
  {
  public:
    ppATLAS(std::string dirname) : fdir(dirname) { init(); }
    ~ppATLAS() { ; }
    const std::vector<std::string> types = {"promptPsi", "nonpromptPsi", "promptX", "nonpromptX", "promptRatio", "nonpromptRatio"};
    const std::vector<std::string> err = {"stat", "syst"};
    std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> gg;
    std::map<std::string, std::map<std::string, TH1F*>> hh;
    TH2F* hempty;
    float getycut() { return ycut; }
    float getweight(float pt, std::string type);
    std::map<std::string, std::string> tpar = {
      std::pair<std::string, std::string>("promptPsi", "Prompt #psi(2S)"),
      std::pair<std::string, std::string>("nonpromptPsi", "Nonprompt #psi(2S)"),
      std::pair<std::string, std::string>("promptX", "Prompt X(3872)"),
      std::pair<std::string, std::string>("nonpromptX", "Nonprompt X(3872)"),
    };
    // std::map<std::string, std::string> formulanum = {
    //   std::pair<std::string, std::string>("promptPsi", ""),
    //   std::pair<std::string, std::string>("nonpromptPsi", ""),
    //   std::pair<std::string, std::string>("promptX", ""),
    //   std::pair<std::string, std::string>("nonpromptX", ""),
    // };
    std::map<std::string, std::string> formula = {
      std::pair<std::string, std::string>("promptPsi", "(0.000017+TMath::Exp(-3.426132-0.202144*x+34.850532/x))/(0.000001+TMath::Exp(-8.128359-0.088624*x+99.751422/x-322.373215/(x*x)))"),
      std::pair<std::string, std::string>("nonpromptPsi", "(0.000033+TMath::Exp(-3.275975-0.174123*x+25.609737/x))/(0.000002+TMath::Exp(-7.853838-0.086783*x+99.999999/x-369.451690/(x*x)))"),
      std::pair<std::string, std::string>("promptX", "(0.000003+TMath::Exp(-2.326617-0.289528*x+4.143284/x))/(0.000000+TMath::Exp(-9.265172-0.131857*x+87.146198/x-184.762029/(x*x)))"),
      std::pair<std::string, std::string>("nonpromptX", "(0.000001+TMath::Exp(-3.713476-0.267920*x+8.699709/x))/(0.000000+TMath::Exp(-9.995278-0.098159*x+74.733928/x-176.118111/(x*x)))"),
    };
    void Draw();
  private:
    const float ycut = 0.75;
    std::string fdir;
    void init();
    std::map<std::string, TF1*> ff;
  };
}

void ppref::ppATLAS::init()
{
  for(auto& tt : types)
    {
      ff[tt] = new TF1(Form("ff_%s", tt.c_str()), formula[tt].c_str());
      std::ifstream inf(fdir+"/"+tt+".dat");
      std::string line;
      std::vector<float> xx, yy, xxel, xxeh, xxelsyst, xxehsyst, yyel_stat, yyeh_stat, yyel_syst, yyeh_syst;
      std::vector<float> xbin, yhist, yehist_stat, yehist_syst;
      while(std::getline(inf, line))
        {
          if(line.empty() || line.length() <= 1) break;
          if(line.front() == '#') continue;

          float x, y, xel, xeh, yel_stat, yeh_stat, yel_syst, yeh_syst;
          std::istringstream ss(line);
          ss >> x >> xeh >> xel >> y >> yeh_stat >> yel_stat >> yeh_syst >> yel_syst;
          if(tt.find("Ratio")==std::string::npos) y *= 1.e-3; yel_stat*=1.e-3; yeh_stat*=1.e-3; yel_syst*=1.e-3; yeh_syst*=1.e-3;

          if(xbin.empty()) xbin.push_back(xel);
          xbin.push_back(xeh);
          yhist.push_back(y);
          yehist_stat.push_back(TMath::Abs(yeh_stat));
          yehist_syst.push_back(TMath::Abs(yeh_syst));

          xx.push_back(x);
          yy.push_back(y);
          xxel.push_back(x-xel);
          xxeh.push_back(xeh-x);
          xxelsyst.push_back(0.5);
          xxehsyst.push_back(0.5);
          yyel_stat.push_back(TMath::Abs(yel_stat));
          yyeh_stat.push_back(TMath::Abs(yeh_stat));
          yyel_syst.push_back(TMath::Abs(yel_syst));
          yyeh_syst.push_back(TMath::Abs(yeh_syst));
        }
      gg[tt]["stat"] = new TGraphAsymmErrors(xx.size(), xx.data(), yy.data(), xxel.data(), xxeh.data(), yyel_stat.data(), yyeh_stat.data());
      gg[tt]["stat"]->SetName(Form("grATLAS_%s_%s", tt.c_str(), "stat"));
      gg[tt]["syst"] = new TGraphAsymmErrors(xx.size(), xx.data(), yy.data(), xxelsyst.data(), xxehsyst.data(), yyel_syst.data(), yyeh_syst.data());
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

float ppref::ppATLAS::getweight(float pt, std::string type)
{
  if(pt < 15 || pt > 50) { return 0; }
  else { return ff[type]->Eval(pt); }
}

void ppref::ppATLAS::Draw()
{
  xjjroot::setthgrstyle(gg["promptRatio"]["syst"], xjjroot::mycolor_middle["azure"], 20, 1.2, 0, 1, 1, xjjroot::mycolor_middle["azure"], 0.3, 1001);
  xjjroot::setthgrstyle(gg["promptRatio"]["stat"], xjjroot::mycolor_middle["azure"], 20, 1.2, xjjroot::mycolor_middle["azure"], 1, 1, xjjroot::mycolor_middle["azure"], 0.3, 1001);
  xjjroot::setthgrstyle(gg["nonpromptRatio"]["syst"], xjjroot::mycolor_middle["azure"], 24, 1.2, 0, 1, 1, xjjroot::mycolor_middle["azure"], 0.3, 1001);
  xjjroot::setthgrstyle(gg["nonpromptRatio"]["stat"], xjjroot::mycolor_middle["azure"], 24, 1.2, xjjroot::mycolor_middle["azure"], 1, 1, xjjroot::mycolor_middle["azure"], 0.3, 1001);

  gg["promptRatio"]["syst"]->Draw("2same");
  gg["promptRatio"]["stat"]->Draw("pesame");
  gg["nonpromptRatio"]["syst"]->Draw("2same");
  gg["nonpromptRatio"]["stat"]->Draw("pesame");
}
