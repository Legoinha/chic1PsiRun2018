#ifndef __OPTCUT__
#define __OPTCUT__ ntp->BDTG[j]>0.7

#include <iostream>
#include <map>
#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>

#include "tnp_weight_lowptPbPb.h"
#include "packtree.h"
#include "ntuple.h"
#include "fitX.h"
#include "xjjcuti.h"

#include "tnpcc_tmp.h"

void tnp_converter(std::string inputname, std::string dirname, int nevt=-1)
{
  std::string outputname = "rootfiles/"+dirname+"/tnp_"+(xjjc::str_divide(inputname, "/").back());
  gSystem->Exec(Form("mkdir -p rootfiles/%s", dirname.c_str()));

  std::cout<<"==> Opening files"<<std::endl;
  std::cout<<"input: "<<inputname<<std::endl;
  std::cout<<"output: "<<outputname<<std::endl;

  std::cout<<"==> Building histograms"<<std::endl;
  std::map<std::string, std::map<std::string, TH1D*>> hh;
  // std::map<std::string, std::map<std::string, TH1D*>> hscales;
  std::map<std::string, std::map<std::string, double>> scales;
  for(auto& tt : tnpcc::types)
    {
      for(auto& idxk : tnpcc::idxname)
        {
          hh[tt][idxk.second] = new TH1D(Form("htnp_%s_%s", tt.c_str(), idxk.second.c_str()), ";p_{T} (GeV/c);", tnpcc::nptbins, tnpcc::ptbins);
          hh[tt][idxk.second]->Sumw2();
          // hhscale[tt][idxk.second] = new TH1D(Form("hscaletnp_%s_%s", tt.c_str(), idxk.second.c_str()), ";#beta;", tnpcc::scalemin, tnpcc::scalemax);
          // hhscale[tt][idxk.second]->Sumw2();
        }
    }
  
  TFile* inf = TFile::Open(inputname.c_str());
  xjjroot::packtree* pt = new xjjroot::packtree(inf, "Bfinder/ntmix", "mcp_tnp");
  mytmva::ntuple* ntp = pt->ntp;

  std::cout<<"==> Scaning file"<<std::endl;
  int nentries = (nevt>0&&nevt<ntp->getnt()->GetEntries())?nevt:ntp->getnt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i); 
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      float weight = ntp->pthatweight * ntp->Ncoll;

      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;

          if(!(TMath::Abs(ntp->By[j]) < fitX::ycut)) continue;
          // ==>
          if(!(__OPTCUT__)) { continue; } 
          // <==
          // ==>
          for(auto& idxk : tnpcc::idxname)
            {
              scales["nominal"][idxk.second] = 1.;
              scales["muid"][idxk.second] = tnp_weight_muid_pbpb(ntp->Bmu1pt[j], ntp->Bmu1eta[j], idxk.first)                 * tnp_weight_muid_pbpb(ntp->Bmu2pt[j], ntp->Bmu2eta[j], idxk.first); 
              scales["trk"][idxk.second]  = tnp_weight_trk_pbpb(ntp->Bmu1eta[j], idxk.first)                                  * tnp_weight_trk_pbpb(ntp->Bmu2eta[j], idxk.first); 
              scales["trg"][idxk.second]  = tnp_weight_trg_pbpb(ntp->Bmu1pt[j], ntp->Bmu1eta[j], tnpcc::filterId, idxk.first) * tnp_weight_trg_pbpb(ntp->Bmu2pt[j], ntp->Bmu2eta[j], tnpcc::filterId, idxk.first); 
              scales["total"][idxk.second] = scales["muid"][idxk.second] * scales["trk"][idxk.second] * scales["trg"][idxk.second];
              //
              for(auto& tt : tnpcc::types)
                {
                  hh[tt][idxk.second]->Fill(ntp->Bpt[j], weight*scales[tt][idxk.second]);
                }
            }
          // <==
        }
    }
  xjjc::progressbar_summary(nentries);

  std::cout<<"==> Writing into output file"<<std::endl;
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  for(auto& ht : hh)
    {
      for(auto& hk : ht.second)
        {
          std::cout<<"\e[2m"<<"writing \e[0m"<<hk.second->GetName()<<" \e[36;1m("<<(int)hk.second->GetEntries()<<")\e[0m"<<std::endl;
          hk.second->Write();
        }
    }
  outf->Close();
  std::cout<<"==> Output file"<<std::endl;
  std::cout<<outputname<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==4) { tnp_converter(argv[1], argv[2], atoi(argv[3])); return 0; }
  if(argc==3) { tnp_converter(argv[1], argv[2]); return 0; }
  return 1;
}

#endif

// +++++++++++++++++++++++++++++++++++++++
// - Trigger: (tnp_weight_trg_pbpb)
//   * filterId = 0: Jpsi L2 filter
//   * filterId = 1: Jpsi L3 filter
//   * filterId = 2: Upsi L2 filter
//   * filterId = 3: Upsi L3 filter
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - MuID: (tnp_weight_muid_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - Inner tracking: (tnp_weight_trk_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// +++++++++++++++++++++++++++++++++++++++
