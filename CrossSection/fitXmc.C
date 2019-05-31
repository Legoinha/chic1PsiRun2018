#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>

#include <string>

void fitXmc(std::string input, std::string inputmc_a, std::string inputmc_b, std::string cut,
            std::string output)
{
  TFile* inf = TFile::Open(input.c_str());
  TTree* ntmix = (TTree*)inf->Get("Bfinder/ntmix");
  ntmix->AddFriend("hiEvtAnalyzer/HiTree");
  ntmix->AddFriend("dataset/mva");

  TFile* infmc_a = TFile::Open(inputmc_a.c_str());
  TTree* ntmixmc_a = (TTree*)infmc_a->Get("Bfinder/ntmix");
  ntmixmc_a->AddFriend("hiEvtAnalyzer/HiTree");
  ntmixmc_a->AddFriend("dataset/mva");

  TFile* infmc_b = TFile::Open(inputmc_b.c_str());
  TTree* ntmixmc_b = (TTree*)infmc_b->Get("Bfinder/ntmix");
  ntmixmc_b->AddFriend("hiEvtAnalyzer/HiTree");
  ntmixmc_b->AddFriend("dataset/mva");

  TH1F* h = new TH1F("h", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
  ntmix->Project("h", "Bmass", TCut(cut.c_str())&&"hiBin<=200");
  TH1F* hmc_a = new TH1F("hmc_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L);
  ntmixmc_a->Project("hmc_a", "Bmass", TCut("1")*(TCut(cut.c_str())&&TCut("Bgen==23333&&hiBin<=200")));
  hmc_a->Scale(hmc_a->GetEntries()/hmc_a->Integral());
  TH1F* hmc_b = new TH1F("hmc_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H);
  ntmixmc_b->Project("hmc_b", "Bmass", TCut("1")*(TCut(cut.c_str())&&TCut("Bgen==23333&&hiBin<=200")));
  hmc_b->Scale(hmc_b->GetEntries()/hmc_b->Integral());

  fitX::fit(h, hmc_a, hmc_b, output);
}

int main(int argc, char* argv[])
{
  if(argc==6) { fitXmc(argv[1], argv[2], argv[3], argv[4], argv[5]); return 0; }
  return 1;
}
