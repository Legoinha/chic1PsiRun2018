#include <TFile.h>
#include <iostream>
#include <TTree.h>

namespace fitX
{
  TTree* getnt(std::string inputfilename, std::string treename = "Bfinder/ntmix", bool addfriend = true)
  {
    std::cout<<__FUNCTION__<<": opening "<<inputfilename<<std::endl;
    TFile* inf = TFile::Open(inputfilename.c_str());
    if(!inf->IsOpen()) return 0;
    TTree* nt = (TTree*)inf->Get(Form("%s", treename.c_str()));
    if(addfriend)
      {
        nt->AddFriend("hiEvtAnalyzer/HiTree");
        nt->AddFriend("hltanalysis/HltTree");
        nt->AddFriend("skimanalysis/HltTree");
        nt->AddFriend("dataset/mva");
      }
    return nt;
  }
}
