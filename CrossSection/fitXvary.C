#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <string>

#include "packtree.h"
#include "ntuple.h"
#include "xjjcuti.h"

#include "fitXvary.h"

void fitXvary(std::string inputdata, std::string inputmc_a, std::string inputmc_b, 
              std::string output)
{

  std::vector<TH1F*> hdata(nbdtg*ndls), hmc_a(nbdtg*ndls), hmc_b(nbdtg*ndls);
  for(int k=0; k<ndls; k++)
    {
      for(int l=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hdata[idx] = new TH1F(Form("hdata_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); 
          hdata[idx]->Sumw2();
          hmc_a[idx] = new TH1F(Form("hmc_a_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); 
          hmc_a[idx]->Sumw2();
          hmc_b[idx] = new TH1F(Form("hmc_b_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); 
          hmc_b[idx]->Sumw2();
        }
    }

  TFile* infdata = TFile::Open(inputdata.c_str());
  if(!infdata->IsOpen()) return;
  TFile* infmc_a = TFile::Open(inputmc_a.c_str());
  if(!infmc_a->IsOpen()) return;
  TFile* infmc_b = TFile::Open(inputmc_b.c_str());
  if(!infmc_b->IsOpen()) return;

  // xjjc::packtree* ptmc_a = new xjjc::packtree(infmc_a, "Bfinder/ntmix", "mc_a");
  // mytmva::ntuple* ntpmc_a = ptmc_a->ntp;
  // xjjc::packtree* ptmc_b = new xjjc::packtree(infmc_b, "Bfinder/ntmix", "mc_b");
  // mytmva::ntuple* ntpmc_b = ptmc_b->ntp;

  xjjc::packtree* pt;
  mytmva::ntuple* ntp;
  int nentries;

  pt = new xjjc::packtree(infdata, "Bfinder/ntmix", "data");
  ntp = pt->ntp;
  nentries = pt->getentries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%100 == 0) { xjjc::progressbar(i, nentries); }

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
                        }
                    }
                }
              else break;
            }
        }
    }


  pt = new xjjc::packtree(infmc_a, "Bfinder/ntmix", "mc_a");
  ntp = pt->ntp;
  nentries = pt->getentries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%100 == 0) { xjjc::progressbar(i, nentries); }

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
                          hmc_a[idx]->Fill(ntp->Bmass[j]);
                        }
                    }
                }
              else break;
            }
        }
    }


  pt = new xjjc::packtree(infmc_b, "Bfinder/ntmix", "mc_b");
  ntp = pt->ntp;
  nentries = pt->getentries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%100 == 0) { xjjc::progressbar(i, nentries); }

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
                          hmc_b[idx]->Fill(ntp->Bmass[j]);
                        }
                    }
                }
              else break;
            }
        }
    }

  // int nentriesmc_a = ptmc_a->getentries();
  // int nentriesmc_b = ptmc_b->getentries();

  TFile* outf = new TFile(Form("rootfiles/root_fitXvary_%s.root", output.c_str()), "recreate");
  outf->cd();
  for(auto& hh : hdata) { hh->Write(); }
  for(auto& hh : hmc_a) { hh->Write(); }
  for(auto& hh : hmc_b) { hh->Write(); }
  outf->Close();
  
}

int main(int argc, char* argv[])
{
  if(argc==5) { fitXvary(argv[1], argv[2], argv[3], argv[4]); return 0; }
  return 1;
}
