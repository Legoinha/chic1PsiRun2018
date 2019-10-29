#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCut.h>
// #include <TSystem.h>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

void fillpthatweight(std::string inputname)
{
  std::string outputname = xjjc::str_divide(inputname, "/").back();
  outputname = xjjc::str_replaceall(outputname, ".root", "");
  outputname = xjjc::str_replaceallspecial(outputname);
  outputname = "rootfiles/pthatweight_"+outputname+".root";

  std::cout<<"==> "<<__FUNCTION__<<": "<<inputname<<std::endl;
  std::cout<<"<== "<<__FUNCTION__<<": "<<outputname<<std::endl;

  TFile* inf = TFile::Open(inputname.c_str());
  // TTree* hi = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  TTree* gen = (TTree*)inf->Get("Bfinder/ntGen");
  gen->AddFriend("hiEvtAnalyzer/HiTree");
  
  TH1F* hpthat = new TH1F("hpthat", ";#hat{p}_{T} (GeV/c);", 80, 0, 80); hpthat->Sumw2();
  TH1F* hpthatweight = new TH1F("hpthatweight", ";#hat{p}_{T} (GeV/c);", 80, 0, 80); hpthatweight->Sumw2();
  TH1F* hgpt = new TH1F("hgpt", ";p_{T} (GeV/c);", 70, 15, 85); hgpt->Sumw2();
  TH1F* hgptweight = new TH1F("hgptweight", ";p_{T} (GeV/c);", 70, 15, 85); hgptweight->Sumw2();

  gen->Project("hpthat", "pthat");
  xjjroot::printhist(hpthat, 14);
  gen->Project("hpthatweight", "pthat", "pthatweight");
  xjjroot::printhist(hpthatweight, 14);

  gen->Project("hgpt", "Gpt", "GisSignal==7");
  xjjroot::printhist(hgpt, 14);
  gen->Project("hgptweight", "Gpt", TCut("pthatweight")*TCut("GisSignal==7"));
  xjjroot::printhist(hgptweight, 14);

  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  hpthat->Write();
  hpthatweight->Write();
  hgpt->Write();
  hgptweight->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==2) { fillpthatweight(argv[1]); return 0; }
  return 1;
}
