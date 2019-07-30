#include <vector>
#include <iostream>
#include <string>
#include <TMath.h>

//============================================================
// type: psi'(type=0), X(type=1), X/psi'(type=2)
// opt: upper("u"), lower("d")
//============================================================

namespace syst
{
  float getsyst(int type, std::string opt);
}

namespace syst
{
  std::vector<float> syst_fit = {5.9, 15.5, 16.6};
  std::vector<float> syst_eff = {32.2, 35.5, 48.0};
  std::vector<float> syst_ptshape = {15.2, 2.4, 15.4};
  std::vector<float> syst_tnp_u = {5.4, 5.3, 0.1};
  std::vector<float> syst_tnp_d = {4.9, 4.8, 0.1};
  std::vector<float> syst_fprompt = {15.7, 3.7, 16.2};
}

float syst::getsyst(int type, std::string opt)
{
  float syst = 0;
  syst += syst_fit[type]*syst_fit[type];
  syst += syst_eff[type]*syst_eff[type];
  syst += syst_ptshape[type]*syst_ptshape[type];
  syst += syst_fprompt[type]*syst_fprompt[type];
  if(opt.find("u") != std::string::npos)
    syst += syst_tnp_u[type]*syst_tnp_u[type];
  else if(opt.find("d") != std::string::npos)
    syst += syst_tnp_d[type]*syst_tnp_d[type];
  
  syst = TMath::Sqrt(syst)*1.e-2;
  std::cout<<"==> "<<type<<" "<<opt<<" "<<syst<<std::endl;
  return syst;
}
