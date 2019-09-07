#include <TFile.h>
#include <TCut.h>
#include <string>
#include "tree_flatten.h"

void fitX_flatten(std::string inputname, std::string outputname="", std::string cut="1")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if(outputname=="")
    outputname = xjjc::str_replaceall(inputname, ".root", "_flatten.root");
  std::cout<<"output: "<<outputname<<std::endl;
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  fitX::tree_flatten* ntf = new fitX::tree_flatten(inputname, "ntmix", outf);
  ntf->flatten();
  // TTree* ntskim = (TTree*)ntf->outnt->CopyTree(TCut(cut.c_str()));
  // ntskim->SetName("ntmix_skim");
  ntf->outnt()->Write();
  xjjroot::printhist(ntf->outnt());
  // ntskim->Write();
  // std::cout<<ntskim->GetEntries()<<std::endl;;  
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==4) { fitX_flatten(argv[1], argv[2], argv[3]); return 0; }
  if(argc==3) { fitX_flatten(argv[1], argv[2]); return 0; }
  if(argc==2) { fitX_flatten(argv[1]); return 0; }
  return 1;
}
