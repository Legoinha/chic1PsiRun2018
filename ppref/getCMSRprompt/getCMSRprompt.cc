#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

int getCMSRprompt()
{
  std::ifstream getdata_psip_fprompt("CMSpsi2Spp7TeV/fNonprompt_Psi2S.dat");
  std::ifstream getdata_psip_xsec("CMSpsi2Spp7TeV/xSection_promptPsi2S.dat");
  std::vector<float> psip_fprompt_ptmin, psip_fprompt_ptmax, psip_fprompt_center, psip_fprompt_stat, psip_fprompt_syst;
  std::vector<float> psip_xsec_ptmin, psip_xsec_ptmax, psip_xsec_center, psip_xsec_stat, psip_xsec_syst, psi_xsec_syst_lumi;
  std::vector<float> psip_fpromptrebin_ptmin, psip_fpromptrebin_ptmax, psip_fpromptrebin_center, psip_fpromptrebin_stat, psip_fpromptrebin_syst, psi_fpromptrebin_syst_lumi;

  while(true)
    {
      float ptmin, ptmax, center, stat, syst;
      getdata_psip_fprompt >> ptmin;
      if(getdata_psip_fprompt.eof()) break;
      getdata_psip_fprompt >> ptmax >> center >> stat >> syst;
      psip_fprompt_ptmin.push_back(ptmin);
      psip_fprompt_ptmax.push_back(ptmax);
      psip_fprompt_center.push_back(center);
      psip_fprompt_stat.push_back(stat);
      psip_fprompt_syst.push_back(syst);
    }

  while(true)
    {
      float ptmin, ptmax, center, stat, syst, syst_lumi;
      getdata_psip_xsec >> ptmin;
      if(getdata_psip_xsec.eof()) break;
      getdata_psip_xsec >> ptmax >> center >> stat >> syst >> syst_lumi;
      psip_xsec_ptmin.push_back(ptmin);
      psip_xsec_ptmax.push_back(ptmax);
      psip_xsec_center.push_back(center);
      psip_xsec_stat.push_back(stat);
      psip_xsec_syst.push_back(syst);
      psip_xsec_syst.push_back(syst_lumi);
    }

  float rebin_ptmin = 10., rebin_ptmax = 13.5;
  float rebin_fprompt_num = 0, rebin_fprompt_den = 0, rebin_stat_num = 0, rebin_syst_num = 0;
  for(int i=0; i<3; i++)
    {
      float Np_i = psip_xsec_center[i]*(psip_xsec_ptmax[i]-psip_xsec_ptmin[i]);
      rebin_fprompt_num += Np_i;
      rebin_fprompt_den += Np_i/psip_fprompt_center[i];
      rebin_stat_num += pow(Np_i, 2) * pow(psip_fprompt_stat[i], 2) / pow(psip_fprompt_center[i], 4);
      rebin_syst_num += pow(Np_i, 2) * pow(psip_fprompt_syst[i], 2) / pow(psip_fprompt_center[i], 4);
    }
  float rebin_center = rebin_fprompt_num / rebin_fprompt_den;
  float rebin_stat = rebin_center * sqrt(rebin_stat_num) / rebin_fprompt_den;
  float rebin_syst = rebin_center * sqrt(rebin_syst_num) / rebin_fprompt_den;
  for(int i=0; i<psip_fprompt_ptmin.size(); i++)
    {
      if(i==1 || i==2) continue;
      if(i==0)
        {
          psip_fpromptrebin_ptmin.push_back(rebin_ptmin);          
          psip_fpromptrebin_ptmax.push_back(rebin_ptmax);          
          psip_fpromptrebin_center.push_back(1-rebin_center);          
          psip_fpromptrebin_stat.push_back(rebin_stat);          
          psip_fpromptrebin_syst.push_back(rebin_syst);          
        }
      else
        {
          psip_fpromptrebin_ptmin.push_back(psip_fprompt_ptmin[i]);
          psip_fpromptrebin_ptmax.push_back(psip_fprompt_ptmax[i]);
          psip_fpromptrebin_center.push_back(1-psip_fprompt_center[i]);
          psip_fpromptrebin_stat.push_back(psip_fprompt_stat[i]);
          psip_fpromptrebin_syst.push_back(psip_fprompt_syst[i]);
        }
    }

  for(int i=0; i<psip_fpromptrebin_ptmin.size(); i++)
    {
      std::cout << std::left 
                << std::setw(10) << psip_fpromptrebin_ptmin[i]
                << std::setw(10) << psip_fpromptrebin_ptmax[i]
                << std::setw(10) << psip_fpromptrebin_center[i]
                << std::setw(10) << psip_fpromptrebin_stat[i]
                << std::setw(10) << psip_fpromptrebin_syst[i]
                << std::endl;
    }

  std::ifstream getdata_xp_fprompt("CMSX3872pp7TeV/fNonprompt_X.dat");
  std::ifstream getdata_xp_R("CMSX3872pp7TeV/R.dat");
  std::vector<float> xp_fprompt_ptmin, xp_fprompt_ptmax, xp_fprompt_center, xp_fprompt_stat, xp_fprompt_syst;
  std::vector<float> xp_R_ptmin, xp_R_ptmax, xp_R_center, xp_R_stat, xp_R_syst;
  std::vector<float> xp_Rprompt_ptmin, xp_Rprompt_ptmax, xp_Rprompt_center, xp_Rprompt_stat, xp_Rprompt_syst;

  while(true)
    {
      float ptmin, ptmax, center, stat, syst;
      getdata_xp_fprompt >> ptmin;
      if(getdata_xp_fprompt.eof()) break;
      getdata_xp_fprompt >> ptmax >> center >> stat >> syst;
      xp_fprompt_ptmin.push_back(ptmin);
      xp_fprompt_ptmax.push_back(ptmax);
      xp_fprompt_center.push_back(1-center);
      xp_fprompt_stat.push_back(stat);
      xp_fprompt_syst.push_back(syst);
    }

  while(true)
    {
      float ptmin, ptmax, center, stat, syst;
      getdata_xp_R >> ptmin;
      if(getdata_xp_R.eof()) break;
      getdata_xp_R >> ptmax >> center >> stat >> syst;
      xp_R_ptmin.push_back(ptmin);
      xp_R_ptmax.push_back(ptmax);
      xp_R_center.push_back(center);
      xp_R_stat.push_back(stat);
      xp_R_syst.push_back(syst);
    }

  for(int i=0; i<psip_fpromptrebin_ptmin.size(); i++)
    {
      xp_Rprompt_ptmin.push_back(xp_R_ptmin[i]);
      xp_Rprompt_ptmax.push_back(xp_R_ptmax[i]);
      xp_Rprompt_center.push_back(xp_R_center[i]*xp_fprompt_center[i]/psip_fpromptrebin_center[i]);
      xp_Rprompt_stat.push_back(sqrt(pow(xp_R_stat[i]/xp_R_center[i], 2) + 
                                     pow(xp_fprompt_stat[i]/xp_fprompt_center[i], 2) + 
                                     pow(psip_fpromptrebin_stat[i]/psip_fpromptrebin_center[i], 2)) * xp_Rprompt_center[i]);
      xp_Rprompt_syst.push_back(sqrt(pow(xp_R_syst[i]/xp_R_center[i], 2) + 
                                     pow(xp_fprompt_syst[i]/xp_fprompt_center[i], 2) + 
                                     pow(psip_fpromptrebin_syst[i]/psip_fpromptrebin_center[i], 2)) * xp_Rprompt_center[i]);      
    }

  std::cout<<std::endl;
  for(int i=0; i<xp_Rprompt_ptmin.size(); i++)
    {
      std::cout << std::left << std::setprecision(4)
                << std::setw(10) << xp_Rprompt_ptmin[i]
                << std::setw(10) << xp_Rprompt_ptmax[i]
                << std::setw(10) << xp_Rprompt_center[i]
                << std::setw(10) << xp_Rprompt_stat[i]
                << std::setw(10) << xp_Rprompt_syst[i]
                << std::setw(10) << xp_R_center[i]
                << std::setw(10) << xp_R_stat[i]
                << std::setw(10) << xp_R_syst[i]
                << std::endl;
    }



  return 0;
}

int main()
{
  return getCMSRprompt();
}
