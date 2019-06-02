#ifndef _XJJC_PACKTREE_H_
#define _XJJC_PACKTREE_H_

#include "ntuple.h"

#include <TFile.h>
#include <TTree.h>
#include <string>

namespace xjjc
{
  typedef std::vector<TTree*> trlist;
  class packtree
  {
  public:
    packtree(TFile* inf, std::string treename, std::string type);
    mytmva::ntuple* ntp;
    int getentries() { return fentries; }
    std::string getname() { return name; }
    void getentry(int i);
    // TTree* getnt() { return nt; }
    // TTree* gethlt() { return hlt; }
    // TTree* getskim() { return skim; }
    // TTree* gethi() { return hi; }
  private:
    std::string name;
    TTree* nt;
    TTree* hlt;
    TTree* skim;
    TTree* hi;
    TTree* mva;
    int fentries;
  };
}

void xjjc::packtree::getentry(int i)
{
  nt->GetEntry(i);
  hlt->GetEntry(i);
  hi->GetEntry(i);
  skim->GetEntry(i);
  mva->GetEntry(i);
}

xjjc::packtree::packtree(TFile* inf, std::string treename, std::string type)
{
  name = type;

  nt = (TTree*)inf->Get(treename.c_str());
  fentries = nt->GetEntries();

  hlt = (TTree*)inf->Get("hltanalysis/HltTree");
  hlt->SetBranchStatus("*", 0);
  hlt->SetBranchStatus("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1", 1);
  nt->AddFriend(hlt);

  skim = (TTree*)inf->Get("skimanalysis/HltTree");
  nt->AddFriend(skim);

  hi = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  nt->AddFriend(hi);

  mva = (TTree*)inf->Get("dataset/mva");
  nt->AddFriend(mva);

  ntp = new mytmva::ntuple(nt);
}

#endif
