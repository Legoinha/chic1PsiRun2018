#ifndef _FITX_TREE_FLATTEN__
#define _FITX_TREE_FLATTEN__

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <string>
#include "xjjcuti.h"

#include "fitX.h"
#include "project.h"

#ifndef MAX_XB
#define MAX_XB      20000
#endif

namespace fitX
{
  template <typename T, size_t N>
  class br_element
  {
  public:
    br_element(std::string brname, bool isarray=true) : exist(false), name(brname), fisarray(isarray) { ; }
    std::string name;
    T value;
    T array[N];
    bool exist;
    bool isarray() { return fisarray; }
  private:
    bool fisarray;
  };

  class tree_flatten
  {
  public:
    tree_flatten(std::string inputname, std::string treename="ntmix", TDirectory* outdir=gDirectory);
    void flatten();
    TTree* outnt() { return foutnt; }

  private:
    TTree* foutnt;
    std::string ffname, fntname;
    TTree* fnt;
    int Bsize;
  
    std::vector<br_element<float, MAX_XB>*> branchesf = {
      new br_element<float, MAX_XB>("Bmass"),
      new br_element<float, MAX_XB>("Btktkmass"),
      new br_element<float, MAX_XB>("Bmumumass"),
      new br_element<float, MAX_XB>("Bujmass"),
      new br_element<float, MAX_XB>("Bpt"), 
      new br_element<float, MAX_XB>("Bgen"), 
      new br_element<float, MAX_XB>("By"),
      new br_element<float, MAX_XB>("Blxy"),
      new br_element<float, MAX_XB>("BDT"),
      new br_element<float, MAX_XB>("BDTG"),
      new br_element<float, MAX_XB>("BDTD"),
      new br_element<float, MAX_XB>("BDTF"),
      new br_element<float, MAX_XB>("BsvpvDistance"),
      new br_element<float, MAX_XB>("BsvpvDisErr"),
      new br_element<float, MAX_XB>("Balpha"),
      new br_element<float, MAX_XB>("Beta"),
      new br_element<float, MAX_XB>("pthatweight", false)
    };
    std::vector<br_element<int, MAX_XB>*> branchesi = {
      new br_element<int, MAX_XB>("BgencollisionId"),
      new br_element<int, MAX_XB>("hiBin", false),
      new br_element<int, MAX_XB>("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1", false),
      new br_element<int, MAX_XB>("pprimaryVertexFilter", false),
      new br_element<int, MAX_XB>("phfCoincFilter2Th4", false),
      new br_element<int, MAX_XB>("pclusterCompatibilityFilter", false),
    };
    std::vector<br_element<bool, MAX_XB>*> brancheso = {
      new br_element<bool, MAX_XB>("mvapref"),
    };
  };

  tree_flatten::tree_flatten(std::string inputname, std::string treename, TDirectory* outdir) : 
    ffname(inputname), 
    fntname(treename)
  {
    fnt = (TTree*)fitX::getnt(ffname, Form("Bfinder/%s", fntname.c_str()));
    // gDirectory->cd("root:/");
    outdir->cd();
    foutnt = new TTree(Form("%s_flatten", fnt->GetName()), "");
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
            foutnt->Branch(br->name.c_str(), &br->value);
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
            foutnt->Branch(br->name.c_str(), &br->value);
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
            foutnt->Branch(br->name.c_str(), &br->value);
          }
      }
  }

  void tree_flatten::flatten()
  {
    int nentries = fnt->GetEntries();
    // foutnt->cd();
    for(int i=0; i<nentries; i++)
      {
        fnt->GetEntry(i);

        if(i%1000==0) xjjc::progressbar(i, nentries);

        for(int j=0; j<Bsize; j++)
          {
            for(auto& br : branchesf) { if(br->isarray()) { br->value = br->array[j]; } }
            for(auto& br : branchesi) { if(br->isarray()) { br->value = br->array[j]; } }
            for(auto& br : brancheso) { if(br->isarray()) { br->value = br->array[j]; } }
            foutnt->Fill();
          }
      }
    // foutnt->Write();
    xjjc::progressbar_summary(nentries);
  }

}

#endif
