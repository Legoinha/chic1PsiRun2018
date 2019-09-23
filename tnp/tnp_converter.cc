#ifndef __OPTCUT__
#define __OPTCUT__ __CUTINPUT__ //ntp->BDTG[j]>0.7

#include <iostream>
#include <map>
#include <string>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TSystem.h>

#include "fitX.h"
#include "tnp_weight_lowptPbPb.h"
#include "packtree.h"
#include "ntuple.h"
#include "xjjcuti.h"

#include "tnpcc_tmp.h"

int n = 5;
void converter(std::string inputname, std::string outputname, std::string name, std::string ptweight="1", std::string wegihtname="no weight");

void tnp_converter(std::string inputname, std::string dirname, std::string name, std::string weightfilename)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* weightfile = TFile::Open(weightfilename.c_str());
  if(weightfile)
    {
      TTree* ntf = (TTree*)weightfile->Get("ntf");
      std::vector<TString*> pars(n), funs(n);
      for(int i=0; i<n; i++)
        {
          pars[i] = 0;
          funs[i] = 0;
          ntf->SetBranchAddress(Form("par%s-%d", name.c_str(), i+1), &(pars[i]));
          ntf->SetBranchAddress(Form("fun%s-%d", name.c_str(), i+1), &(funs[i]));
        }
      ntf->GetEntry(0);
      for(int i=0; i<n; i++)
        {
          std::string weight = xjjc::str_replaceall(xjjc::str_tolower(pars[i]->Data()), "tmath::", "");
          weight = xjjc::str_replaceall(weight, "x[0]", "x");
          std::string outputname = "rootfiles/"+dirname+fitX::tagname()+"/funs/fun-"+std::string(Form("%d", i+1))+"/"+"tnp"+name+".root";
          converter(inputname, outputname, name, weight, funs[i]->Data());
        }
    }
  else
    {
      std::string outputname = "rootfiles/"+dirname+fitX::tagname()+"/tnp"+name+".root";
      converter(inputname, outputname, name);
    }

}

void converter(std::string inputname, std::string outputname, std::string name, std::string ptweight, std::string weightname)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  std::cout<<"==> Checking pt weight"<<std::endl;
  std::cout<<ptweight<<std::endl;
  TF1* fptweight = new TF1("fptweight", ptweight.c_str(), fitX::ptmincut, fitX::ptmaxcut);

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
  TH2F* hptweight = new TH2F("hptweight", Form(";p_{T} (GeV/c);%s", weightname.c_str()), 50, tnpcc::ptbins[0], tnpcc::ptbins[tnpcc::nptbins-1], 50, 0, 5);
  
  TFile* inf = TFile::Open(inputname.c_str());
  xjjroot::packtree* pt = new xjjroot::packtree(inf, "Bfinder/ntmix", "mcp_tnp");
  mytmva::ntuple* ntp = pt->ntp();

  std::cout<<"==> Scaning file"<<std::endl;
  int nentries = ntp->nt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i); 
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      float weight = ntp->pthatweight * ntp->Ncoll;

      if(!(ntp->hiBin >= fitX::centmincut && ntp->hiBin <= fitX::centmaxcut)) continue;
      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;

          if(!(TMath::Abs(ntp->By[j])>=fitX::ymincut && TMath::Abs(ntp->By[j])<fitX::ymaxcut)) continue;
          // ==>
          if(!(__OPTCUT__)) { continue; } 
          // <==
          float wptweight = fptweight->Eval(ntp->Bgenpt[j]);
          hptweight->Fill(ntp->Bgenpt[j], wptweight, weight);
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
                  hh[tt][idxk.second]->Fill(ntp->Bpt[j], weight*scales[tt][idxk.second]*wptweight);
                }
            }
          // <==
        }
    }
  xjjc::progressbar_summary(nentries);

  std::cout<<"==> Writing into output file"<<std::endl;
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  for(auto& ht : hh)
    {
      for(auto& hk : ht.second)
        {
          xjjroot::printhist(hk.second, 20);
          hk.second->Write();
        }
    }
  hptweight->Write();
  fitX::write();
  outf->Close();
  std::cout<<"==> Output file"<<std::endl;
  std::cout<<outputname<<std::endl;
  std::cout<<std::endl;

  gROOT->cd();
  for(auto& ht : hh)
    { for(auto& hk : ht.second) { delete hk.second; } }
  delete hptweight;
}

int main(int argc, char* argv[])
{
  if(argc==11) { 
    fitX::init(atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]));
    tnp_converter(argv[1], argv[2], argv[3], argv[4]); return 0; }
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

