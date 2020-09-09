#ifndef __DATAMC_VAR_H_
#define __DATAMC_VAR_H_

#include <vector>
#include <map>
#include "xjjrootuti.h"
#include "xjjcuti.h"

#include "fitX.h"

namespace datamc
{
  class var
  {
  public:
    var(std::string _type) : ftype(_type) { fvalid = init(); }
    std::vector<double> vars() { return fvars; }
    std::string type() { return ftype; }
    std::string formula() { return fformula; }
    std::string title() { return ftitle; }
    std::string unit() { return funit; }
    bool gt() { return fgt; }
    bool valid() { return fvalid; }
    int n() { return fvars.size(); }
    int nfine() { return fnfine; }
    double minfine() { return fminfine; }
    double maxfine() { return fmaxfine; }
  private:
    bool init();
    bool fvalid;
    bool fgt;
    std::vector<double> fvars;
    std::string ftype;
    std::string fformula;
    std::string ftitle;
    std::string funit;
    int fnfine;
    double fminfine, fmaxfine;
  };
}

bool datamc::var::init()
{
  if(ftype=="BDT")
    {
      fformula = "BDT";
      fvars = std::vector<double>({0.06, 0.08, 0.10, 0.12, 0.14});
      fgt = true;
      ftitle = "BDT";
      funit = "";
      fnfine = 40;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="Qvalue")
    {
      fformula = "(Bmass-3.096916-Btktkmass)";
      // fformula = "(Bmass-Bmumumass-Btktkmass)";
      fvars = std::vector<double>({0., 0.04, 0.08, 0.2});
      fgt = false;
      ftitle = "Q = m_{#mu#mu#pi#pi}-m_{#mu#mu}-m_{#pi#pi}";
      funit = "(GeV/c^{2})";
      fnfine = 40;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="pt")
    {
      fformula = "Bpt";
      // fvars = std::vector<double>({15., 20., 25., 50.});
      fvars = std::vector<double>({15., 17, 20., 24., 28., 50.});
      // fvars = std::vector<double>({15., 18, 20., 23., 27., 50.});
      // fvars = std::vector<double>({15., 18, 21., 24., 28., 50.});
      fgt = true;
      ftitle = "p_{T}";
      funit = "(GeV/c)";
      fnfine = 35;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="absy")
    {
      // fformula = "TMath::Abs(By)";
      fformula = "fabs(By)";
      fvars = std::vector<double>({0, 0.4, 0.8, 1.6});
      fgt = false;
      ftitle = "|y|";
      funit = "";
      fnfine = 32;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="chi2cl")
    {
      fformula = "Bchi2cl";
      fvars = std::vector<double>({0.1, 0.4, 0.7, 1});
      fgt = true;
      ftitle = "#chi^{2} prob";
      funit = "";
      fnfine = 30;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="trk1pt")
    {
      fformula = "Btrk1Pt";
      fvars = std::vector<double>({1, 2, 3, 4});
      fgt = true;
      ftitle = "trk1 p_{T}";
      funit = "(GeV/c)";
      fnfine = 30;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="trk2pt")
    {
      fformula = "Btrk2Pt";
      fvars = std::vector<double>({1, 2, 3, 4, 5});
      fgt = true;
      ftitle = "trk2 p_{T}";
      funit = "(GeV/c)";
      fnfine = 40;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="trkptimb")
    {
      fformula = "fabs(Btrk1Pt-Btrk2Pt)/(Btrk1Pt+Btrk2Pt)";
      fvars = std::vector<double>({0, 0.2, 0.5, 1});
      fgt = true;
      ftitle = "|p_{T,#pi1}-p_{T,#pi2}|/(p_{T,#pi1}+p_{T,#pi2})";
      funit = "";
      fnfine = 50;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="dRtrk1")
    {
      fformula = "sqrt(pow(acos(cos(Bujphi-Btrk1Phi)),2) + pow(Bujeta-Btrk1Eta,2))";
      fvars = std::vector<double>({0, 0.1, 0.2, 0.4});
      fgt = false;
      ftitle = "#DeltaR (trk1)";
      funit = "";
      fnfine = 40;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="dRtrk2")
    {
      fformula = "sqrt(pow(acos(cos(Bujphi-Btrk2Phi)),2) + pow(Bujeta-Btrk2Eta,2))";
      fvars = std::vector<double>({0, 0.1, 0.2, 0.4});
      fgt = false;
      ftitle = "#DeltaR (trk2)";
      funit = "";
      fnfine = 40;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="BDTnoQ")
    {
      fformula = "BDT";
      fvars = std::vector<double>({-0.4, -0.3, -0.2, -0.1, -0.02, 0.02, 0.06, 0.10, 0.16, 0.30});
      fgt = true;
      ftitle = "BDT";
      funit = "";
      fnfine = 60;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else if(ftype=="lxy")
    {
      fformula = "10*Blxy*Bmass/Bpt";
      fvars = std::vector<double>({-0.1, -0.04, 0, 0.04, 0.14, 0.3, 0.5, 0.7, 1.1});
      fgt = true;
      ftitle = "l_{xy}";
      funit = "(mm)";
      fnfine = 55;
      fminfine = fvars.front();
      fmaxfine = fvars.back();
    }
  else { return false; }
  return true;
}


#endif
