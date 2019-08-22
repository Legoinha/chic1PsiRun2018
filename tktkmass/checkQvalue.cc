#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCut.h>

#include <string>
#include <iostream>

#include "project.h"
#include "Qvalue.h"
#include "fitX.h"
#include "xjjrootuti.h"

void checkQvalue(std::string inputmc, int itype, std::string name)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TH1F* hq = new TH1F("hq", ";Q value = m_{#mu#mu#pi#pi}-m_{#mu#mu}-m_{#pi#pi} (GeV/c^{2});", Qvalue::NBIN, Qvalue::BIN_MIN, Qvalue::BIN_MAX);
  TH1F* hptimb = new TH1F("hptimb", ";|p_{T,#pi1}-p_{T,#pi2}|/(p_{T,#pi1}+p_{T,#pi2});", ptimb::NBIN, ptimb::BIN_MIN, ptimb::BIN_MAX);
  TTree* ntmixmc = fitX::getnt(inputmc, "Bfinder/ntmix"); if(!ntmixmc) { return; }
  std::string cut = "Bgen>=23333 && BgencollisionId==0 && Bpt > 15 && Bpt < 50 && TMath::Abs(By)>=0 && TMath::Abs(By)<1.6 && hiBin>=0 && hiBin<=190";
  std::string mcweight = "pthatweight * Ncoll";

  gROOT->cd();
  std::cout<<" == input ==>"<<std::endl;
  std::cout<<cut<<std::endl;
  std::cout<<Qvalue::Q.c_str()<<std::endl;
  std::cout<<" == fill ==>"<<std::endl;
  ntmixmc->Project("hq", Qvalue::Q.c_str(), TCut(mcweight.c_str())*TCut(cut.c_str()));
  fitX::printhist(hq);
  ntmixmc->Project("hptimb", "TMath::Abs(Btrk1Pt-Btrk2Pt)/(Btrk1Pt+Btrk2Pt)", TCut(mcweight.c_str())*TCut(cut.c_str()));
  fitX::printhist(hptimb);
  
  std::string outputname(Form("rootfiles/%s/Qvalue_%s.root", name.c_str(), Qvalue::types[itype].c_str()));
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = TFile::Open(outputname.c_str(), "recreate");
  hq->Write();
  hptimb->Write();
  outf->Close();
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==4) { checkQvalue(argv[1], atoi(argv[2]), argv[3]); return 0; }
  return 1;
}
