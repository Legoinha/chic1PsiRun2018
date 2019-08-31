#include "xjjrootuti.h"

#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include <string>

namespace MCeff
{
  std::vector<float> ptBins = {15, 20, 25, 30, 35, 40, 45, 50};
  // std::vector<float> ptBins = {15, 25, 35, 45, 55};
  const int nPtBins = ptBins.size() - 1;

  std::vector<float> ptBins_incl = {15, 50};
  const int nPtBins_incl = ptBins_incl.size() - 1;

  class MCefficiency
  {
  public:
    MCefficiency(std::string name, int whichincl=0, std::vector<float> _ptbins = MCeff::ptBins);
    MCefficiency(TFile* inf, std::string name, int whichincl=0);

    TH1F* heffmc() { return fheffmc; }
    TH1F* heffmc_incl() { return fheffmc_incl; }
    TH1F* heffgen() { return fheffgen; }
    TH1F* heffgen_incl() { return fheffgen_incl; }
    TEfficiency* greff() { return fgreff; }
    TEfficiency* greff_incl() { return fgreff_incl; }
    TH1F* heff() { return fheff; }
    TH1F* heff_incl() { return fheff_incl; }

    void calceff();
    void setstyle(Color_t color, Style_t mstyle=20, Style_t lstyle=2);
    std::vector<float> ptbins() { return fptbins; }
    int nbins() { return fnptbins; }
  private:
    int fincl;
    std::string fname;
    void createhist();
    void readhist(TFile* inf);
    std::vector<float> fptbins;
    int fnptbins;
    TH1F* fheffmc;
    TH1F* fheffmc_incl;
    TH1F* fheffgen;
    TH1F* fheffgen_incl;
    TEfficiency* fgreff;
    TEfficiency* fgreff_incl;
    TH1F* fheff;
    TH1F* fheff_incl;
    std::vector<float> fptbins_incl;
    int fnptbins_incl;
  };
}

MCeff::MCefficiency::MCefficiency(std::string name, int whichincl, std::vector<float> _ptbins) : fname(name), fincl(whichincl), fptbins(_ptbins), fnptbins(_ptbins.size()-1)
{
  createhist();
  fptbins_incl = MCeff::ptBins_incl;
  fnptbins_incl = MCeff::nPtBins_incl;
}

MCeff::MCefficiency::MCefficiency(TFile* inf, std::string name, int whichincl) : fname(name), fincl(whichincl)
{
  readhist(inf);
  fptbins_incl = MCeff::ptBins_incl;
  fnptbins_incl = MCeff::nPtBins_incl;
}

void MCeff::MCefficiency::createhist()
{
  fheffmc = new TH1F(Form("heffmc%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data()); fheffmc->Sumw2();
  if(fincl==0) { fheffmc_incl = new TH1F(Form("heffmc_incl%s", fname.c_str()), "", 5, 0, 5); fheffmc_incl->Sumw2(); }
  if(fincl==1) { fheffmc_incl = new TH1F(Form("heffmc_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffmc_incl->Sumw2(); }
  fheffgen = new TH1F(Form("heffgen%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data());
  if(fincl==0) { fheffgen_incl = new TH1F(Form("heffgen_incl%s", fname.c_str()), "", 5, 0, 5); fheffgen_incl->Sumw2(); }
  if(fincl==1) { fheffgen_incl = new TH1F(Form("heffgen_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffgen_incl->Sumw2(); }
}

void MCeff::MCefficiency::readhist(TFile* inf)
{
  fheffmc = (TH1F*)inf->Get(Form("heffmc%s", fname.c_str()));
  fheffmc_incl = (TH1F*)inf->Get(Form("heffmc_incl%s", fname.c_str()));
  fheffgen = (TH1F*)inf->Get(Form("heffgen%s", fname.c_str()));
  fheffgen_incl = (TH1F*)inf->Get(Form("heffgen_incl%s", fname.c_str()));
  fnptbins = fheffmc->GetXaxis()->GetNbins();
  fptbins.clear();
  for(int i=0; i<fnptbins+1; i++)
    {
      fptbins.push_back((*(fheffmc->GetXaxis()->GetXbins()))[i]);
    }
}

void MCeff::MCefficiency::calceff()
{
  fgreff = new TEfficiency(*fheffmc, *fheffgen); fgreff->SetName(Form("greff%s", fname.c_str()));
  fgreff_incl = new TEfficiency(*fheffmc_incl, *fheffgen_incl); fgreff_incl->SetName(Form("greff_incl%s", fname.c_str()));
  fheff = (TH1F*)fheffmc->Clone(Form("heff%s", fname.c_str()));
  fheff->Divide(fheffgen);
  fheff_incl = (TH1F*)fheffmc_incl->Clone(Form("heff_incl%s", fname.c_str()));
  fheff_incl->Divide(fheffgen_incl);
}

void MCeff::MCefficiency::setstyle(Color_t color, Style_t mstyle/*=20*/, Style_t lstyle/*=2*/)
{
  xjjroot::setthgrstyle(fheff, 0, 0, 0, color, lstyle, 3, color, 0.1, 1001);
  xjjroot::setthgrstyle(fheff_incl, color, mstyle, 1.2, color, lstyle, 3, color, 0.1, 1001);
  xjjroot::setthgrstyle(fgreff, 0, 0, 0, color, lstyle, 3, color, 0.1, 1001);
  xjjroot::setthgrstyle(fgreff_incl, color, mstyle, 1.2, color, lstyle, 3, color, 0.1, 1001);
}
