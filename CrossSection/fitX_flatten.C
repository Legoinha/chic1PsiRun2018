#include "TFile.h"
#include <string>
#include "tree_flatten.h"

void fitX_flatten(std::string inputname, std::string outputname="")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if(outputname=="")
    outputname = xjjc::str_replaceall(inputname, ".root", "_flatten.root");
  fitX::tree_flatten ntf(inputname);
  ntf.flatten();
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  ntf.outnt->Write();
  std::cout<<"output: "<<otuputname<<std::endl;
  std::cout<<ntf.outnt->GetEntries()<<std::endl;;
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX_flatten(argv[1], argv[2]); return 0; }
  if(argc==2) { fitX_flatten(argv[1]); return 0; }
  return 1;
}
