#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>

#include <string>

float weight_ss = 1./0.681;

void fitXmc(std::string input, std::string input_ss, std::string inputmc_a, std::string inputmc_b, 
            std::string cut, std::string output)
{

  std::cout<<cut<<std::endl;

  bool is_ss = false;

  std::cout<<" -- data >>>>"<<std::endl;
  TFile* inf = TFile::Open(input.c_str());
  if(!inf->IsOpen()) return;
  TH1F* h = new TH1F("h", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h->Sumw2();
  TTree* ntmix = (TTree*)inf->Get("Bfinder/ntmix");
  ntmix->AddFriend("hiEvtAnalyzer/HiTree");
  ntmix->AddFriend("hltanalysis/HltTree");
  ntmix->AddFriend("skimanalysis/HltTree");
  ntmix->AddFriend("dataset/mva");
  ntmix->Project("h", "Bmass", TCut(Form("(%s)", cut.c_str())));
  std::cout<<h->GetEntries()<<std::endl;
  std::cout<<" -- same-sign >>>>"<<std::endl;
  // TFile* inf_ss = TFile::Open(input_ss.c_str());
  TH1F* h_ss = new TH1F("h_ss", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h_ss->Sumw2();
  // if(inf_ss->IsOpen()) 
  //   {  
  //     is_ss = true;
  //     TTree* ntmix_ss = (TTree*)inf_ss->Get("Bfinder/ntmix");
  //     ntmix_ss->AddFriend("hiEvtAnalyzer/HiTree");
  //     ntmix_ss->AddFriend("hltanalysis/HltTree");
  //     ntmix_ss->AddFriend("skimanalysis/HltTree");
  //     ntmix_ss->AddFriend("dataset/mva");
  //     ntmix_ss->Project("h_ss", "Bmass", TCut(Form("(%s)", cut.c_str())));
  //     std::cout<<h_ss->GetEntries()<<std::endl;
  //     h_ss->Scale(weight_ss);
  //   }

  std::cout<<" -- mc_a >>>>"<<std::endl;
  TFile* infmc_a = TFile::Open(inputmc_a.c_str());
  if(!infmc_a->IsOpen()) return;
  TH1F* hmc_a = new TH1F("hmc_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); hmc_a->Sumw2();
  TTree* ntmixmc_a = (TTree*)infmc_a->Get("Bfinder/ntmix");
  ntmixmc_a->AddFriend("hiEvtAnalyzer/HiTree");
  ntmixmc_a->AddFriend("hltanalysis/HltTree");
  ntmixmc_a->AddFriend("skimanalysis/HltTree");
  ntmixmc_a->AddFriend("dataset/mva");
  ntmixmc_a->Project("hmc_a", "Bmass", TCut("1")*TCut(Form("(%s) && Bgen==23333 && BgencollisionId==0", cut.c_str())));
  std::cout<<hmc_a->GetEntries()<<std::endl;
  hmc_a->Scale(hmc_a->GetEntries()/hmc_a->Integral());

  std::cout<<" -- mc_b >>>>"<<std::endl;
  TFile* infmc_b = TFile::Open(inputmc_b.c_str());
  if(!infmc_b->IsOpen()) return;
  TH1F* hmc_b = new TH1F("hmc_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); hmc_b->Sumw2();
  TTree* ntmixmc_b = (TTree*)infmc_b->Get("Bfinder/ntmix");
  ntmixmc_b->AddFriend("hiEvtAnalyzer/HiTree");
  ntmixmc_b->AddFriend("hltanalysis/HltTree");
  ntmixmc_b->AddFriend("skimanalysis/HltTree");
  ntmixmc_b->AddFriend("dataset/mva");
  ntmixmc_b->Project("hmc_b", "Bmass", TCut("1")*TCut(Form("(%s) && Bgen==23333 && BgencollisionId==0", cut.c_str())));
  std::cout<<hmc_b->GetEntries()<<std::endl;
  hmc_b->Scale(hmc_b->GetEntries()/hmc_b->Integral());

  fitX::fit(h, (is_ss?h_ss:0), hmc_a, hmc_b, output);

  TFile* outf = new TFile(Form("rootfiles/root_fitXmc_%s.root", output.c_str()), "recreate");
  outf->cd();
  h->Write();
  h_ss->Write();
  hmc_a->Write();
  hmc_b->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==7) { fitXmc(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); return 0; }
  return 1;
}
