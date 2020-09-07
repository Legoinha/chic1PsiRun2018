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
void converter(std::string inputname, std::string outputname, std::string name, std::string ptweight="1");

void checkL2L3(std::string inputname, std::string dirname, std::string name, std::string ptweightfile)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::string outputname = "rootfiles/"+dirname+fitX::tagname()+"/muL2L3"+name+".root";
  converter(inputname, outputname, name);
}

void converter(std::string inputname, std::string outputname, std::string name, std::string ptweight)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  std::cout<<"==> Checking pt weight"<<std::endl;
  std::cout<<ptweight<<std::endl;
  TF1* fptweight = new TF1("fptweight", ptweight.c_str(), fitX::ptmincut, fitX::ptmaxcut);

  std::cout<<"==> Opening files"<<std::endl;
  std::cout<<"input: "<<inputname<<std::endl;
  std::cout<<"output: "<<outputname<<std::endl;

  std::cout<<"==> Building histograms"<<std::endl;
  TH1D* hL2L3 = new TH1D("hL2L3", ";;", 3, 0, 3);
  TH2D* hL2L3pair = new TH2D("hL2L3pair", ";#mu^{(1)}; #mu^{(2)}", 3, 0, 3, 3, 0, 3);
  TH2D* hL2L3ver = new TH2D("hL2L3ver", ";Aug 2019;Aug 2020", 2, 0, 2, 2, 0, 2);
  
  TFile* inf = TFile::Open(inputname.c_str());
  xjjroot::packtree* pt = new xjjroot::packtree(inf, "Bfinder/ntmix", "mcp_tnp");
  mytmva::ntuple* ntp = pt->ntp();

  std::cout<<"==> Scaning file"<<std::endl;
  int nentries = ntp->nt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i); 
      if(i%1000 == 0 || i==nentries-1) { xjjc::progressbar(i, nentries); }

      float weight = ntp->pthatweight * ntp->Ncoll;

      if(!(ntp->hiBin >= fitX::centmincut && ntp->hiBin <= fitX::centmaxcut)) continue;
      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          // if(!ntp->mvapref[j]) continue;
          if(!ntp->passedpre(j)) continue;
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;

          if(!(TMath::Abs(ntp->By[j])>=fitX::ymincut && TMath::Abs(ntp->By[j])<fitX::ymaxcut)) continue;
          // ==>
          if(!(__OPTCUT__)) { continue; } 
          // <==
          
          int mu1ver2019 = 0, mu1ver2020 = 0;
          if(ntp->Bmu1isTriggered[j]) mu1ver2019 = 1;
          if(ntp->mu1hltL2L3lastflt[j] >= 0) mu1ver2020 = 1;
          hL2L3ver->Fill(mu1ver2019, mu1ver2020, weight);
          int mu2ver2019 = 0, mu2ver2020 = 0;
          if(ntp->Bmu2isTriggered[j]) mu2ver2019 = 1;
          if(ntp->mu2hltL2L3lastflt[j] >= 0) mu2ver2020 = 1;
          hL2L3ver->Fill(mu2ver2019, mu2ver2020, weight);

          // variation !!
          if(!(ntp->mu1hltL2L3lastflt[j] >= 0 && ntp->mu2hltL2L3lastflt[j] >= 0)) continue;

          int mu1type = 0;
          if(ntp->mu1hltL3flt[j] >= 0) mu1type = tnpcc::filterId["jpsiL3"]+1;
          else if(ntp->mu1hltL2flt[j] >= 0) mu1type = tnpcc::filterId["jpsiL2"]+1;
          hL2L3->Fill(mu1type, weight);

          int mu2type = 0;
          if(ntp->mu2hltL3flt[j] >= 0) mu2type = tnpcc::filterId["jpsiL3"]+1;
          else if(ntp->mu2hltL2flt[j] >= 0) mu2type = tnpcc::filterId["jpsiL2"]+1;
          hL2L3->Fill(mu2type, weight);

          hL2L3pair->Fill(mu1type, mu2type, weight);
        }
    }
  xjjc::progressbar_summary(nentries);

  std::cout<<"==> Writing into output file"<<std::endl;
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hL2L3->Write();
  hL2L3pair->Write();
  hL2L3ver->Write();
  fitX::write();
  outf->Close();
  std::cout<<"==> Output file"<<std::endl;
  std::cout<<outputname<<std::endl;
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==11) { 
    fitX::init(atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]));
    checkL2L3(argv[1], argv[2], argv[3], argv[4]); return 0; }
  return 1;
}

#endif

