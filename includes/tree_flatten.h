#ifndef _FITX_TREE_FLATTEN__
#define _FITX_TREE_FLATTEN__

#include <TFile.h>
#include <TTree.h>
#include <string>
#include "xjjcuti.h"

#include "fitX.h"
#include "project.h"

#ifndef MAX_XB
#define MAX_XB      20000
#endif

namespace fitX
{
  template <typename T>
  class br_element
  {
  public:
    br_element(std::string brname, bool isarray=true) : exist(false), name(brname), fisarray(isarray) { ; }
    std::string name;
    T value;
    T array[MAX_XB];
    bool exist;
    bool isarray() { return fisarray; }
  private:
    bool fisarray;
  };

  class tree_flatten
  {
  public:
    tree_flatten(std::string inputname, std::string treename="ntmix");
    void flatten();
    TTree* outnt;

  private:
    std::string ffname, fntname;
    TTree* fnt;
    int Bsize;
  
    std::vector<br_element<float>*> branchesf = {
      new br_element<float>("Bmass"),
      new br_element<float>("Btktkmass"),
      new br_element<float>("Bpt"), 
      new br_element<float>("Bgen"), 
      new br_element<float>("By"),
      new br_element<float>("Blxy"),
      new br_element<float>("BDT"),
      new br_element<float>("BDTG"),
      new br_element<float>("BDTD"),
      new br_element<float>("BDTF"),
      new br_element<float>("pthatweight", false)
    };
    std::vector<br_element<int>*> branchesi = {
      new br_element<int>("BgencollisionId"),
      new br_element<int>("hiBin", false),
      new br_element<int>("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1", false),
      new br_element<int>("pprimaryVertexFilter", false),
      new br_element<int>("phfCoincFilter2Th4", false),
      new br_element<int>("pclusterCompatibilityFilter", false),
    };
    std::vector<br_element<bool>*> brancheso = {
      new br_element<bool>("mvapref"),
    };
  };

  tree_flatten::tree_flatten(std::string inputname, std::string treename) : 
    ffname(inputname), 
    fntname(treename)
  {
    fnt = (TTree*)fitX::getnt(ffname, Form("Bfinder/%s", fntname.c_str()));
    gDirectory->cd("root:/");
    outnt = new TTree(Form("%s_flatten", fnt->GetName()), "");
    fnt->SetBranchAddress("Bsize", &Bsize);
    for(auto& br : branchesf)
      {
        int exitcode = -5;
        if(br->isarray()) exitcode = fnt->SetBranchAddress(br->name.c_str(), br->array);
        else exitcode = fnt->SetBranchAddress(br->name.c_str(), &br->value);
        if(exitcode == -5) { std::cout<<__FUNCTION__<<": warning: no branch \e[4m"<<br->name<<"\e[0m"<<std::endl; }
        else 
          { 
            br->exist = true; 
            outnt->Branch(br->name.c_str(), &br->value);
          }
      }
    for(auto& br : branchesi)
      {
        int exitcode = -5;
        if(br->isarray()) exitcode = fnt->SetBranchAddress(br->name.c_str(), br->array);
        else exitcode = fnt->SetBranchAddress(br->name.c_str(), &br->value);
        if(exitcode == -5) { std::cout<<__FUNCTION__<<": warning: no branch \e[4m"<<br->name<<"\e[0m"<<std::endl; }
        else 
          { 
            br->exist = true; 
            outnt->Branch(br->name.c_str(), &br->value);
          }
      }
    for(auto& br : brancheso)
      {
        int exitcode = -5;
        if(br->isarray()) exitcode = fnt->SetBranchAddress(br->name.c_str(), br->array);
        else exitcode = fnt->SetBranchAddress(br->name.c_str(), &br->value);
        if(exitcode == -5) { std::cout<<__FUNCTION__<<": warning: no branch \e[4m"<<br->name<<"\e[0m"<<std::endl; }
        else 
          { 
            br->exist = true; 
            outnt->Branch(br->name.c_str(), &br->value);
          }
      }
  }

  void tree_flatten::flatten()
  {
    int nentries = fnt->GetEntries();
    // outnt->cd();
    for(int i=0; i<nentries; i++)
      {
        fnt->GetEntry(i);

        if(i%1000==0) xjjc::progressbar(i, nentries);

        for(int j=0; j<Bsize; j++)
          {
            for(auto& br : branchesf) { if(br->isarray()) { br->value = br->array[j]; } }
            for(auto& br : branchesi) { if(br->isarray()) { br->value = br->array[j]; } }
            for(auto& br : brancheso) { if(br->isarray()) { br->value = br->array[j]; } }
            outnt->Fill();
          }
      }
    // outnt->Write();
    xjjc::progressbar_summary(nentries);
  }

}

#endif
