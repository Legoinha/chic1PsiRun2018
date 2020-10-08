#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <TString.h>
#include "xjjcuti.h"

#include "getdata.h"

int rebin()
{
  ppref::getdata g_fprompt("CMSpsi2Spp7TeV/rebin_fNonprompt_Psi2S.dat", true);
  ppref::getdata g_xsec("CMSpsi2Spp7TeV/rebin_xsec_promptPsi2S.dat", false);
  // ppref::getdata g_fprompt_all("CMSpsi2Spp7TeV/rebin_fNonprompt_Psi2S.dat", true);
  // ppref::getdata g_xsec_all("CMSpsi2Spp7TeV/rebin_xsec_promptPsi2S.dat", false);
  std::cout<<std::endl;
  
  std::map<std::string, std::vector<float>> rebin_fprompt, rebin_xsec;
  int n = g_fprompt.n();
  float sum_dpt_xsec_p=0, sum_dpt_xsec=0, sum_fprompt_stat=0, sum_fprompt_syst=0, sum_xsec_stat=0, sum_xsec_syst=0;
  for(int i=0; i<n; i++)
    {
      float dpt_xsec_p = g_xsec["center"][i]*g_xsec.binwidth(i);
      float dpt_xsec = dpt_xsec_p / g_fprompt["center"][i];

      sum_dpt_xsec_p += dpt_xsec_p;
      sum_dpt_xsec += dpt_xsec;
      sum_fprompt_stat += dpt_xsec * g_fprompt["stat"][i];
      sum_fprompt_syst += dpt_xsec * g_fprompt["syst"][i];
      sum_xsec_stat += dpt_xsec_p * g_xsec["stat"][i];
      sum_xsec_syst += dpt_xsec_p * g_xsec["syst"][i];
    }

  rebin_fprompt["ptmin"].push_back(g_fprompt["ptmin"].front());
  rebin_fprompt["ptmax"].push_back(g_fprompt["ptmax"].back());
  rebin_xsec["ptmin"].push_back(g_xsec["ptmin"].front());
  rebin_xsec["ptmax"].push_back(g_xsec["ptmax"].back());
  rebin_fprompt["center"].push_back(sum_dpt_xsec_p / sum_dpt_xsec);
  rebin_xsec["center"].push_back(sum_dpt_xsec_p / (g_xsec["ptmax"].back() - g_xsec["ptmin"].front()));
  rebin_fprompt["stat"].push_back(sum_fprompt_stat / sum_dpt_xsec);
  rebin_fprompt["syst"].push_back(sum_fprompt_syst / sum_dpt_xsec);
  rebin_xsec["stat"].push_back(sum_xsec_stat / sum_dpt_xsec_p);
  rebin_xsec["syst"].push_back(sum_xsec_syst / sum_dpt_xsec_p);

  ppref::getdata g_rebin_fprompt(rebin_fprompt, 1, "rebin_fprompt");
  ppref::getdata g_rebin_xsec(rebin_xsec, 1, "rebin_xsec");

  // print
  std::cout<<"\e[36;1m"; g_fprompt.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[36;1m"; g_xsec.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[33;1m"; g_rebin_fprompt.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[33;1m"; g_rebin_xsec.print(); std::cout<<"\e[0m"<<std::endl;

}

int main()
{
  return rebin();
}
