#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <TString.h>
#include "xjjcuti.h"

#include "getdata.h"

int getCMSRprompt()
{
  ppref::getdata g_fprompt_X("CMSX3872pp7TeV/fNonprompt_X.dat", true);
  ppref::getdata g_R("CMSX3872pp7TeV/R.dat", false);
  ppref::getdata g_fprompt_Psi2S("CMSpsi2Spp7TeV/fNonprompt_Psi2S.dat", true);
  ppref::getdata g_xsec_promptPsi2S("CMSpsi2Spp7TeV/xsec_promptPsi2S.dat", false, 2.4); // |y| < 1.2

  ppref::getdata g_xsec_promptX("CMSX3872pp7TeV/xsec_promptX.dat", false);
  float BRjpsipipi = 0.34, BRee = 7.82e-3;
  float relerrBRjpsipipi = 0.004/BRjpsipipi, relerrBRee = 0.17e-3/BRee;

  std::cout<<std::endl;

  int n = g_R.n();
  std::map<std::string, std::vector<float>> Rprompt, xsec_promptX_my;

  for(int i=0; i<n; i++)
    {
      xsec_promptX_my["ptmin"].push_back(g_R["ptmin"][i]);
      xsec_promptX_my["ptmax"].push_back(g_R["ptmax"][i]);
      Rprompt["ptmin"].push_back(g_R["ptmin"][i]);
      Rprompt["ptmax"].push_back(g_R["ptmax"][i]);

      // Rprompt
      float center = (g_fprompt_X["center"][i] / g_fprompt_Psi2S["center"][i]) * g_R["center"][i];
      Rprompt["center"].push_back(center);

      float rel_stat_sq = 
        pow(g_fprompt_X["stat_rel"][i], 2) +
        pow(g_R["stat_rel"][i], 2) +
        pow(g_fprompt_Psi2S["stat_rel"][i], 2);
      float stat = sqrt(rel_stat_sq) * center;
      Rprompt["stat"].push_back(stat);

      float rel_syst_sq = 
        pow(g_fprompt_X["syst_rel"][i], 2) +
        pow(g_R["syst_rel"][i], 2) +
        pow(g_fprompt_Psi2S["syst_rel"][i], 2);
      float syst = sqrt(rel_syst_sq) * center;
      Rprompt["syst"].push_back(syst);

      // prompt X xsec
      center = center * g_xsec_promptPsi2S["center"][i] * BRjpsipipi / BRee;
      xsec_promptX_my["center"].push_back(center);

      rel_stat_sq += pow(g_xsec_promptPsi2S["stat_rel"][i], 2); //
      stat = sqrt(rel_stat_sq) * center;
      xsec_promptX_my["stat"].push_back(stat);

      rel_syst_sq = rel_syst_sq + 
        pow(g_xsec_promptPsi2S["syst_rel"][i], 2) + //
        pow(relerrBRjpsipipi, 2) +
        pow(relerrBRee, 2);
      syst = sqrt(rel_syst_sq) * center;
      xsec_promptX_my["syst"].push_back(syst);

    }
  ppref::getdata g_xsec_promptX_my(xsec_promptX_my, n, "xsec_promptX_my");
  ppref::getdata g_Rprompt(Rprompt, n, "Rprompt");

  // print
  std::cout<<"\e[36;1m"; g_fprompt_X.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[36;1m"; g_fprompt_Psi2S.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[36;1m"; g_R.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[36;1m"; g_xsec_promptPsi2S.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[33;1m"; g_xsec_promptX_my.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[32;1m"; g_xsec_promptX.print(); std::cout<<"\e[0m"<<std::endl;
  std::cout<<"\e[33;1m"; g_Rprompt.print(); std::cout<<"\e[0m"<<std::endl;
  return 0;
}

int main()
{
  return getCMSRprompt();
}
