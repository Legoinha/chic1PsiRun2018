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
    std::vector<float> vars() { return fvars; }
    std::string type() { return ftype; }
    std::string formula() { return fformula; }
    std::string title() { return ftitle; }
    bool gt() { return fgt; }
    bool valid() { return fvalid; }
    int n() { return fvars.size(); }
  private:
    bool init();
    bool fvalid;
    bool fgt;
    std::vector<float> fvars;
    std::string ftype;
    std::string fformula;
    std::string ftitle;
  };
}

bool datamc::var::init()
{
  if(ftype=="BDT")
    {
      fformula = "BDT";
      fvars = std::vector<float>({0.06, 0.08, 0.10, 0.12, 0.14});
      fgt = true;
      ftitle = "BDT";
    }
  else if(ftype=="Qvalue")
    {
      fformula = "(Bmass-3.096916-Btktkmass)";
      fvars = std::vector<float>({0., 0.07, 0.10, 0.2});
      fgt = false;
      ftitle = "m_{#mu#mu#pi#pi}-m_{#mu#mu}-m_{#pi#pi} (GeV/c^{2})";
    }
  else if(ftype=="pt")
    {
      fformula = "Bpt";
      fvars = std::vector<float>({15., 20., 50.});
      fgt = true;
      ftitle = "p_{T} (GeV/c)";
    }

  else if(ftype=="absy")
    {
      fformula = "TMath::Abs(By)";
      fvars = std::vector<float>({0, 0.4, 0.8, 1.6});
      fgt = false;
      ftitle = "|y|";
    }
  else { return false; }
  return true;
}


#endif
