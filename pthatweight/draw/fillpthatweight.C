#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TSystem.h>
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
  TTree* hi = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  
  TH1F* hpthat = new TH1F("hpthat", ";#hat{p}_{T} (GeV/c);", 80, 0, 80); hpthat->Sumw2();
  TH1F* hpthatweight = new TH1F("hpthatweight", ";#hat{p}_{T} (GeV/c);", 80, 0, 80); hpthatweight->Sumw2();

  hi->Project("hpthat", "pthat");
  hi->Project("hpthatweight", "pthat", "pthatweight");

  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  hpthat->Write();
  hpthatweight->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==2) { fillpthatweight(argv[1]); return 0; }
  return 1;
}
