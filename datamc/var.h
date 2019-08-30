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
  private:
    bool init();
    bool fvalid;
    bool fgt;
    std::vector<double> fvars;
    std::string ftype;
    std::string fformula;
    std::string ftitle;
    std::string funit;
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
    }
  else if(ftype=="Qvalue")
    {
      fformula = "(Bmass-3.096916-Btktkmass)";
      // fformula = "(Bmass-Bmumumass-Btktkmass)";
      fvars = std::vector<double>({0., 0.04, 0.08, 0.2});
      fgt = false;
      ftitle = "Q = m_{#mu#mu#pi#pi}-m_{#mu#mu}-m_{#pi#pi}";
      funit = "(GeV/c^{2})";
    }
  else if(ftype=="pt")
    {
      fformula = "Bpt";
      fvars = std::vector<double>({15., 20., 50.});
      fgt = true;
      ftitle = "p_{T}";
      funit = "(GeV/c)";
    }

  else if(ftype=="absy")
    {
      // fformula = "TMath::Abs(By)";
      fformula = "fabs(By)";
      fvars = std::vector<double>({0, 0.4, 0.8, 1.6});
      fgt = false;
      ftitle = "|y|";
      funit = "";
    }
  else { return false; }
  return true;
}


#endif
