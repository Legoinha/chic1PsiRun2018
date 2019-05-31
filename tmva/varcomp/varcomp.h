#ifndef _VARCOMP_H_
#define _VARCOMP_H_

#include "ntuple.h"
#include "TMVAClassification.h"
#include <TMath.h>
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>

#include <vector>
#include <map>
#include <string>
#include <iostream>

#ifndef MASS_JPSI
#define MASS_JPSI  3.096916
#endif

namespace varcomp
{
  class variable
  {
  public:
    int nvar;
    std::vector<TH1F*> hist;
    const int nbin = 50;

    variable(mytmva::ntuple* nt_, std::string name_);
    variable(TFile* histfile, std::string name_);
    ~variable() { for(auto& h : hist) { delete h; } } 
    void fillhist(int j);
    void savehist() { for(auto& h : hist) { std::cout<<"--- "<<__FUNCTION__<<":  \e[2mTH1F\e[0m "<<h->GetName()<<" \e[2m("<<h->GetEntries()<<")\e[0m"<<" Write."<<std::endl; h->Write(); } }
    void sethist(int k);
    void clearnt() { fnt = 0; }
    bool isvalid() { return fvalid; }

  private:
    const std::vector<mytmva::tmvavar> varlist = mytmva::varlist;
    std::string name;
    mytmva::ntuple* fnt; //~
    TFile* finf; //~
    mytmva::varval* fvarval;
    bool fvalid;

  };

  std::vector<Color_t> lcolor = {kRed, kBlue-3, kGreen+2};
  std::vector<Color_t> fcolor = {kRed, kBlue-3, kGreen+2};
  std::vector<Style_t> fstyle = {3005, 1001, 1001};
  std::vector<float>   falpha = {1, 0.4, 0.4};

  variable::variable(mytmva::ntuple* nt_, std::string name_) : fnt(nt_), name(name_), finf(0), fvalid(true)
  { 
    nvar = varlist.size();
    fvarval = new mytmva::varval(nt_);
    hist.resize(nvar, 0);
    for(int i=0; i<nvar; i++)
      {
        hist[i] = new TH1F(Form("h_%s_%s_%d", name.c_str(), varlist[i].varname.c_str(), i), Form(";%s;", varlist[i].vartex.c_str()), nbin, varlist[i].varmin, varlist[i].varmax); 
        hist[i]->Sumw2();
      }
    std::cout<<"--- "<<__FUNCTION__<<":  \e[1m"<<name<<"\e[0m ("<<nvar<<") created. => pthatweight "<<(fnt->isweight()?"applied.":"\e[31;1mnot\e[0m applied.")<<std::endl;
  }
  variable::variable(TFile* histfile, std::string name_) : fnt(0), name(name_), finf(histfile), fvalid(true), fvarval(0)
  {
    nvar = varlist.size();
    hist.resize(nvar, 0);
    for(int i=0; i<nvar; i++)
      {
        hist[i] = (TH1F*)finf->Get(Form("h_%s_%s_%d", name.c_str(), varlist[i].varname.c_str(), i));
        if(hist[i]==0) { fvalid = false; break; }
        std::cout<<"--- "<<__FUNCTION__<<":  \e[2mTH1F\e[0m "<<hist[i]->GetName()<<" \e[2m("<<hist[i]->GetEntries()<<")\e[0m"<<" attached."<<std::endl;
      }
    if(fvalid) std::cout<<"--- "<<__FUNCTION__<<":  \e[1m"<<name<<"\e[0m created ("<<nvar<<")."<<std::endl;    
    else std::cout<<"--- "<<__FUNCTION__<<":  \e[1m"<<name<<"\e[0m fails."<<std::endl;    
  }

  void variable::fillhist(int j)
  {
    float weight = fnt->isweight()?fnt->pthatweight:1.;
    for(int i=0; i<nvar; i++) { hist[i]->Fill(fvarval->getval(varlist[i].varname, j), weight); }
  }
  void variable::sethist(int k)
  {
    for(auto& h : hist)
      {
        h->Scale(1./h->Integral());
        h->SetLineColor(lcolor[k]);
        h->SetFillStyle(fstyle[k]);
        h->SetFillColor(fcolor[k]);
        h->SetFillColorAlpha(fcolor[k], falpha[k]);
        h->SetLineWidth(2);
      }
  }

  //
  

}

#endif
