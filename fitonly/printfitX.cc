#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>
#include "xjjcuti.h"

void printfitX(std::string outputdir, std::string fitopts="default")
{
  std::vector<std::string> options = xjjc::str_divide(fitopts, ",");
  int nopt = options.size();

  std::cout<<"\e[36;1m"<<std::endl;
  UInt_t stdwid = 15;
  std::cout<<std::string(stdwid*5+1, '-')<<std::endl;
  std::cout<<std::left<<std::setw(stdwid)<<"| Title"<<std::setw(stdwid)<<"| Sig (PLL)"<<std::setw(stdwid)<<"| Sig (Null)"<<std::setw(stdwid)<<"| p value"<<std::setw(stdwid)<<"| Yield"<<"|"<<std::endl;
  std::cout<<std::string(stdwid*5+1, '-')<<std::endl;
  for(int i=0; i<nopt; i++)
    {
      std::string outputname = "rootfiles/" + outputdir + "/" + options[i] + "/paramnll.root";
      TFile* inf = TFile::Open(outputname.c_str());
      TTree* param = (TTree*)inf->Get("param");
      float pllsig; param->SetBranchAddress("pllsig", &pllsig);
      float nullsig; param->SetBranchAddress("nullsig", &nullsig);
      float nullpvalue; param->SetBranchAddress("nullpvalue", &nullpvalue);
      float xmeanpll; param->SetBranchAddress("xmeanpll", &xmeanpll);
      param->GetEntry(0);
      std::cout<<std::left<<std::setw(stdwid)<<Form("| %s", options[i].c_str())<<std::setw(stdwid)<<Form("| %.4f", pllsig)<<std::setw(stdwid)<<Form("| %.4f", nullsig)<<std::setw(stdwid)<<Form("| %.4fe-5", nullpvalue*1.e+5)<<std::setw(stdwid)<<Form("| %.1f", xmeanpll)<<"|"<<std::endl;
      std::cout<<std::string(stdwid*5+1, '-')<<std::endl;
    }
  std::cout<<"\e[0m"<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==3) { printfitX(argv[1], argv[2]); return 0; }
  return 1;
}
