#include <vector>
#include <string>
#include <TH1F.h>
#include <TEfficiency.h>

#include "packtree.h"
#include "ntuple.h"
#include "xjjrootuti.h"
#include "xjjcuti.h"

// std::vector<float> bdtg = {-1, -0.9, -0.8, 0, 0.4, 0.6, 0.70, 0.76, 0.80}; int nbdtg = bdtg.size();
std::vector<float> bdtg = {-1   , -0.8 , -0.6 , -0.4, -0.2 , -0.1 , 0   , 0.1  , 0.2  , 0.3  , 0.4 , 
                           0.5  , 0.55 , 0.6  , 0.65 , 0.70 , 0.75 , 0.80, 0.85 , 0.9, 1.0};
std::vector<bool> pbdtg = {true , true , true, false, false, false, true, false, false, false, true, 
                           false, false, true , false, false, true , true, false, false, false};
int nbdtg = bdtg.size();
std::vector<float> dls = {0, 0.8}; int ndls = dls.size();

void drawalltext()
{
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.042, "X(3872)", 0.038, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.90, 0.84, "p_{T} > 15 GeV/c", 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.042, "|y| < 1.5", 0.038, 32, 62);
}


//
std::vector<TH1F*> hdata(nbdtg*ndls), hmc_a(nbdtg*ndls), hmc_b(nbdtg*ndls), hdatagt(nbdtg*ndls);
std::vector<TH1F*> heffmc_a(ndls), heffmc_b(ndls);
std::vector<TH1F*> hlxymcnp_a(nbdtg*ndls), hlxymcnp_b(nbdtg*ndls);
TH1F* heffgen_a;
TH1F* heffgen_b;
std::vector<TH1F*> heff_a(ndls), heff_b(ndls);
std::vector<TEfficiency*> greff_a(ndls), greff_b(ndls);
//
std::vector<TF1*> ff(nbdtg*ndls, 0), ffgt(nbdtg*ndls, 0);

namespace fitX
{
  void loop_vary_data(xjjroot::packtree* pt, mytmva::ntuple* ntp, int nentries);
  void loop_vary_mcprompt(xjjroot::packtree* pt, mytmva::ntuple* ntp, int nentries, std::string name);
  void loop_vary_mcnonprompt(xjjroot::packtree* pt, mytmva::ntuple* ntp, int nentries, std::string name);
}

void fitX::loop_vary_data(xjjroot::packtree* pt, mytmva::ntuple* ntp, int nentries)
{
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bpt[j] > 15 && TMath::Abs(ntp->By[j]) < 1.5)) continue;
          for(int l=0; l<nbdtg; l++)
            {
              if(ntp->BDTG[j] > bdtg[l])
                {
                  for(int k=0; k<ndls; k++)
                    {
                      if(ntp->BsvpvDistance[j]/ntp->BsvpvDisErr[j] > dls[k])
                        {
                          int idx = k*nbdtg + l;
                          hdata[idx]->Fill(ntp->Bmass[j]);
                          if(ntp->Blxy[j] > 0.1)
                            {
                              hdatagt[idx]->Fill(ntp->Bmass[j]);
                            }
                        }
                    }
                }
              else break;
            }
        }
    }
  xjjc::progressbar_summary(nentries);
}

void fitX::loop_vary_mcprompt(xjjroot::packtree* pt, mytmva::ntuple* ntp, int nentries, std::string name)
{
  std::vector<std::vector<TH1F*>> hmcs     = {hmc_a,     hmc_b};
  std::vector<std::vector<TH1F*>> heffmcs  = {heffmc_a,  heffmc_b};
  std::vector<TH1F*>              heffgens = {heffgen_a, heffgen_b};
  int isam = xjjc::str_contains(name, "_a")?0:1;
  std::vector<TH1F*> hmc = hmcs[isam], heffmc = heffmcs[isam];
  TH1F* heffgen = heffgens[isam];
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
          if(!(ntp->Gpt[j] > 15 && TMath::Abs(ntp->Gy[j]) < 1.5)) continue;
          for(int l=0; l<nbdtg; l++)
            {
              heffgen->Fill(heffgen->GetBinCenter(l+1), weight);
            }
        }

      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bpt[j] > 15 && TMath::Abs(ntp->By[j]) < 1.5)) continue;
          if(!(ntp->Bgen[j]==23333 && ntp->BgencollisionId[j]==0)) continue;
          for(int l=0; l<nbdtg; l++)
            {
              if(ntp->BDTG[j] > bdtg[l])
                {
                  for(int k=0; k<ndls; k++)
                    {
                      if(ntp->BsvpvDistance[j]/ntp->BsvpvDisErr[j] > dls[k])
                        {
                          int idx = k*nbdtg + l;
                          hmc[idx]->Fill(ntp->Bmass[j]);
                          heffmc[k]->Fill(heffmc[k]->GetBinCenter(l+1), weight); // weight!
                        }
                    }
                }
              else break;
            }
        }
    }
  xjjc::progressbar_summary(nentries);
}
  
void fitX::loop_vary_mcnonprompt(xjjroot::packtree* pt, mytmva::ntuple* ntp, int nentries, std::string name)
{
  std::vector<std::vector<TH1F*>> hlxymcnps = {hlxymcnp_a, hlxymcnp_b};
  int isam = xjjc::str_contains(name, "_a")?0:1;
  std::vector<TH1F*> hlxymcnp = hlxymcnps[isam];
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
          if(!(ntp->Bpt[j] > 15 && TMath::Abs(ntp->By[j]) < 1.5)) continue;
          if(!(ntp->Bgen[j]==23333 && ntp->BgencollisionId[j]==0)) continue;
          for(int l=0; l<nbdtg; l++)
            {
              if(ntp->BDTG[j] > bdtg[l])
                {
                  for(int k=0; k<ndls; k++)
                    {
                      if(ntp->BsvpvDistance[j]/ntp->BsvpvDisErr[j] > dls[k])
                        {
                          int idx = k*nbdtg + l;
                          hlxymcnp[idx]->Fill(ntp->Blxy[j], weight);
                        }
                    }
                }
              else break;
            }
        }
    }
  xjjc::progressbar_summary(nentries);
}
