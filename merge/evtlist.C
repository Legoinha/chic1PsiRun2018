#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>

#include "xjjcuti.h"

void evtlist(std::string inputname, std::string outputname)
{
  TFile* inf = TFile::Open(inputname.c_str());
  TTree* hi = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  UInt_t run; hi->SetBranchAddress("run", &run);
  ULong64_t evt; hi->SetBranchAddress("evt", &evt);
  UInt_t lumi;hi->SetBranchAddress("lumi", &lumi);
  
  int nentries = hi->GetEntries();
  std::ofstream evtlistfile(outputname.c_str());
  
  for(int i=0; i<nentries; i++)
    {
      hi->GetEntry(i);
      if(i%10000==0) xjjc::progressbar(i, nentries);

      evtlistfile<<run<<" "<<lumi<<" "<<evt<<std::endl;
    }
  xjjc::progressbar_summary(nentries);
}

int main(int argc, char* argv[])
{
  if(argc==3) { evtlist(argv[1], argv[2]); return 0; }
  return 1;
}
