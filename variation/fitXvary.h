#ifndef __FITX_VARYCUT_H_
#define __FITX_VARYCUT_H_

#include <vector>
#include <string>
#include <TH1F.h>
#include <TEfficiency.h>

#include "packtree.h"
#include "ntuple.h"
#include "xjjcuti.h"

#include "fitX.h"
#include "lxydis.h"

namespace fitX
{
  template<typename T>
  class varycut
  {
  public:
    std::vector<TH1F*> hdata, hdataBenr; // mass
    std::vector<TH1F*> hmc_a, hmc_b; // MC mass
    std::vector<TH1F*> hlxymcnp_a, hlxymcnp_b, hlxymcp_a, hlxymcp_b; // MC np/p lxy dis
    TH1F *heffmc_a, *heffmc_b;
    TH1F *heffgen_a, *heffgen_b;
    TH1F *heff_a, *heff_b;
    TEfficiency *greff_a, *greff_b;
    TH1F *hsideband_a, *hsideband_b; // sideband
    TH1F *hmcdisp_a, *hmcdisp_b, *hmcdisnp_a, *hmcdisnp_b;
    varycut(std::vector<T> vvary, std::string titlevary, int iscut) : fvv(vvary), fvt(titlevary), fiscutordis(iscut)
    {
      if(fiscutordis) fnv = vvary.size();
      else fnv = vvary.size() - 1;
      hdata.resize(fnv);
      hmc_a.resize(fnv);
      hmc_b.resize(fnv);
      hdataBenr.resize(fnv);
      hlxymcnp_a.resize(fnv);
      hlxymcnp_b.resize(fnv);
      hlxymcp_a.resize(fnv);
      hlxymcp_b.resize(fnv);
      lxyxbins = lxydis::setupbins();
    }
    const float masswinL = 0.2, masswinH = 0.4;

    int producehist();
    int loop_vary_data(xjjroot::packtree* pt, int nentries);
    int loop_vary_mcprompt(xjjroot::packtree* pt, int nentries, std::string name);
    int loop_vary_mcnonprompt(xjjroot::packtree* pt, int nentries, std::string name);
    int getnv() { return fnv; }

    int produceeff();

  private:
    int fnv;
    std::vector<T> fvv;
    std::string fvt;
    int fiscutordis;
    std::vector<TH1F*> fhmc, fhlxymcnp, fhlxymcp;
    TH1F *fheffmc, *fheffgen, *fheff;
    TEfficiency *fgreff;
    TH1F *fhmcdisp, *fhmcdisnp;
    std::map<std::string, std::vector<float>> lxyxbins;
    int setaorb(std::string name) 
    {
      if(xjjc::str_contains(name, "_a"))
        {
          fhmc = hmc_a;
          fheffmc = heffmc_a;
          fheffgen = heffgen_a;
          fheff = heff_a;
          fgreff = greff_a;
          fhlxymcnp = hlxymcnp_a;
          fhlxymcp = hlxymcp_a;
          fhmcdisp = hmcdisp_a;
          fhmcdisnp = hmcdisnp_a;
          return 0;
        }
      if(xjjc::str_contains(name, "_b")) 
        {
          fhmc = hmc_b;
          fheffmc = heffmc_b;
          fheffgen = heffgen_b;
          fheff = heff_b;
          fgreff = greff_b;
          fhlxymcnp = hlxymcnp_b;
          fhlxymcp = hlxymcp_b;
          fhmcdisp = hmcdisp_b;
          fhmcdisnp = hmcdisnp_b;
          return 1;
        }
      return -1;
    }
  };
}

template<typename T>
int fitX::varycut<T>::loop_vary_data(xjjroot::packtree* pt, int nentries)
{
  mytmva::ntuple* ntp = pt->ntp;
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          
          if(!(ntp->Bpt[j] > fitX::ptcut && TMath::Abs(ntp->By[j]) < fitX::ycut)) continue;
          if(!xjjc::str_contains(fvt, "BDT")) { if(ntp->BDTG[j]<0.7) { continue; } }

          for(int l=0; l<fnv; l++)
            {
              if((fiscutordis && __VARIABLE__ > fvv[l]) ||
                 (!fiscutordis && __VARIABLE__ >= fvv[l] && __VARIABLE__ < fvv[l+1]))
                {
                  hdata[l]->Fill(ntp->Bmass[j]);
                  if(TMath::Abs(ntp->Bmass[j]-MASS_PSI2S) > masswinL && TMath::Abs(ntp->Bmass[j]-MASS_PSI2S) < masswinH) { hsideband_a->Fill(hsideband_a->GetBinCenter(l+1)); }
                  if(TMath::Abs(ntp->Bmass[j]-MASS_X) > masswinL && TMath::Abs(ntp->Bmass[j]-MASS_X) < masswinH) { hsideband_b->Fill(hsideband_b->GetBinCenter(l+1)); }
                  if(ntp->Blxy[j] > 0.1)
                    {
                      hdataBenr[l]->Fill(ntp->Bmass[j]);
                    }
                  if(!fiscutordis) break;
                }
              else if(fiscutordis) { break; }
            }
        }
    }
  xjjc::progressbar_summary(nentries);
}

template<typename T>
int fitX::varycut<T>::loop_vary_mcprompt(xjjroot::packtree* pt, int nentries, std::string name)
{
  int isam = setaorb(name);
  if(isam < 0) return 1;
  mytmva::ntuple* ntp = pt->ntp;

  //
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!(ntp->hiBin >=0 && ntp->hiBin<180)) continue;
      float weight = ntp->pthatweight * ntp->Ncoll;

      for(int j=0; j<ntp->Gsize; j++)
        {
          if(!(ntp->GisSignal[j]==7 && ntp->GcollisionId[j]==0)) continue;
          if(!(ntp->Gpt[j] > fitX::ptcut && TMath::Abs(ntp->Gy[j]) < fitX::ycut)) continue;
          for(int l=0; l<fnv; l++)
            {
              fheffgen->Fill(fheffgen->GetBinCenter(l+1), weight);
            }
        }

      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;

          if(!(ntp->Bpt[j] > fitX::ptcut && TMath::Abs(ntp->By[j]) < fitX::ycut)) continue;
          if(!xjjc::str_contains(fvt, "BDT")) { if(ntp->BDTG[j]<0.7) { continue; } }
          for(int l=0; l<fnv; l++)
            {
              if((fiscutordis && __VARIABLE__ > fvv[l]) ||
                 (!fiscutordis && __VARIABLE__ >= fvv[l] && __VARIABLE__ < fvv[l+1]))
                {
                  fhmc[l]->Fill(ntp->Bmass[j]);
                  fheffmc->Fill(fheffmc->GetBinCenter(l+1), weight); // weight!
                  fhlxymcp[l]->Fill(ntp->Blxy[j], weight); // weight!
                  fhmcdisp->Fill(__VARIABLE__, weight);
                  if(!fiscutordis) break;
                }
              else if(fiscutordis) { break; }
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  return 0;
}

template<typename T>
int fitX::varycut<T>::loop_vary_mcnonprompt(xjjroot::packtree* pt, int nentries, std::string name)
{
  int isam = setaorb(name);
  if(isam < 0) return 1;
  mytmva::ntuple* ntp = pt->ntp;

  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!(ntp->hiBin >=0 && ntp->hiBin<180)) continue;
      float weight = ntp->pthatweight * ntp->Ncoll;
      if(!ntp->passedevtfil()) continue;

      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bpt[j] > fitX::ptcut && TMath::Abs(ntp->By[j]) < fitX::ycut)) continue;
          if(!xjjc::str_contains(fvt, "BDT")) { if(ntp->BDTG[j]<0.7) { continue; } }
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;
          for(int l=0; l<fnv; l++)
            {
              if((fiscutordis && __VARIABLE__ > fvv[l]) ||
                 (!fiscutordis && __VARIABLE__ >= fvv[l] && __VARIABLE__ < fvv[l+1]))
                {
                  fhlxymcnp[l]->Fill(ntp->Blxy[j], weight); // weight!
                  fhmcdisnp->Fill(__VARIABLE__, weight);
                  if(!fiscutordis) break;
                }
              else if(fiscutordis) { break; }
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  return 0;
}

template<typename T>
int fitX::varycut<T>::producehist()
{
  for(int l=0; l<fnv; l++)
    {
      hdata[l] = new TH1F(Form("hdata_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
      hdata[l]->Sumw2();
      hdataBenr[l] = new TH1F(Form("hdataBenr_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
      hdataBenr[l]->Sumw2();
      hmc_a[l] = new TH1F(Form("hmc_a_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L);
      hmc_a[l]->Sumw2();
      hmc_b[l] = new TH1F(Form("hmc_b_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H);
      hmc_b[l]->Sumw2();
      hlxymcnp_a[l] = new TH1F(Form("hlxymcnp_a_%d", l), ";l_{xy} (mm);Probability", lxyxbins["lxynonprompt"].size()-1, lxyxbins["lxynonprompt"].data());
      hlxymcnp_a[l]->Sumw2();
      hlxymcnp_b[l] = new TH1F(Form("hlxymcnp_b_%d", l), ";l_{xy} (mm);Probability", lxyxbins["lxynonprompt"].size()-1, lxyxbins["lxynonprompt"].data());
      hlxymcnp_b[l]->Sumw2();
      hlxymcp_a[l] = new TH1F(Form("hlxymcp_a_%d", l), ";l_{xy} (mm);Probability", lxyxbins["lxyprompt"].size()-1, lxyxbins["lxyprompt"].data());
      hlxymcp_a[l]->Sumw2();
      hlxymcp_b[l] = new TH1F(Form("hlxymcp_b_%d", l), ";l_{xy} (mm);Probability", lxyxbins["lxyprompt"].size()-1, lxyxbins["lxyprompt"].data());
      hlxymcp_b[l]->Sumw2();
    }
  heffmc_a = new TH1F("heffmc_a", Form(";%s;#alpha #times #epsilon", fvt.c_str()), fvv.size()-1, fvv.data());
  heffmc_b = new TH1F("heffmc_b", Form(";%s;#alpha #times #epsilon", fvt.c_str()), fvv.size()-1, fvv.data());
  heffgen_a = new TH1F("heffgen_a", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  heffgen_b = new TH1F("heffgen_b", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hsideband_a = new TH1F("hsideband_a", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hsideband_b = new TH1F("hsideband_b", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hmcdisp_a = new TH1F("hmcdisp_a", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hmcdisp_b = new TH1F("hmcdisp_b", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hmcdisnp_a = new TH1F("hmcdisnp_a", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hmcdisnp_b = new TH1F("hmcdisnp_b", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
}

template<typename T>
int fitX::varycut<T>::produceeff()
{
  heff_a = (TH1F*)heffmc_a->Clone("heff_a");
  heff_b = (TH1F*)heffmc_b->Clone("heff_b");
  greff_a = new TEfficiency(*heffmc_a, *heffgen_a); greff_a->SetName("greff_a");
  greff_b = new TEfficiency(*heffmc_b, *heffgen_b); greff_b->SetName("greff_b");
  for(int i=0;i<heff_a->GetNbinsX();i++)
    {
      heff_a->SetBinContent(i+1, greff_a->GetEfficiency(i+1));
      heff_a->SetBinError(i+1, greff_a->GetEfficiencyErrorUp(i+1));
      heff_b->SetBinContent(i+1, greff_b->GetEfficiency(i+1));
      heff_b->SetBinError(i+1, greff_b->GetEfficiencyErrorUp(i+1));
    }
  return 0;
}

#endif
