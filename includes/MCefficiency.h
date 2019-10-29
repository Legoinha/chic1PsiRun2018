#include "xjjrootuti.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
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
    TH1F* heffgenacc() { return fheffgenacc; }
    TH1F* heffgenacc_incl() { return fheffgenacc_incl; }
    TH1F* heffmcpre() { return fheffmcpre; }
    TH1F* heffmcpre_incl() { return fheffmcpre_incl; }
    TH1F* heffmcbdt() { return fheffmcbdt; }
    TH1F* heffmcbdt_incl() { return fheffmcbdt_incl; }
    TEfficiency* greff() { return fgreff; }
    TEfficiency* greff_incl() { return fgreff_incl; }
    TEfficiency* gracc() { return fgracc; }
    TEfficiency* gracc_incl() { return fgracc_incl; }
    TEfficiency* greffpre() { return fgreffpre; }
    TEfficiency* greffpre_incl() { return fgreffpre_incl; }
    TEfficiency* greffbdt() { return fgreffbdt; }
    TEfficiency* greffbdt_incl() { return fgreffbdt_incl; }
    TEfficiency* greffqvl() { return fgreffqvl; }
    TEfficiency* greffqvl_incl() { return fgreffqvl_incl; }
    TEfficiency* greffcut() { return fgreffcut; }
    TEfficiency* greffcut_incl() { return fgreffcut_incl; }
    TH1F* heff() { return fheff; }
    TH1F* heff_incl() { return fheff_incl; }

    void calceff();
    void calcacc();
    void setstyle(Color_t color, Style_t mstyle=20, Style_t lstyle=2, Width_t lwidth=3);
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
    TH1F* fheffgenacc;
    TH1F* fheffgenacc_incl;
    TH1F* fheffmcpre;
    TH1F* fheffmcpre_incl;
    TH1F* fheffmcbdt;
    TH1F* fheffmcbdt_incl;
    TEfficiency* fgreff;
    TEfficiency* fgreff_incl;
    TEfficiency* fgracc;
    TEfficiency* fgracc_incl;
    TEfficiency* fgreffpre;
    TEfficiency* fgreffpre_incl;
    TEfficiency* fgreffbdt;
    TEfficiency* fgreffbdt_incl;
    TEfficiency* fgreffqvl;
    TEfficiency* fgreffqvl_incl;
    TEfficiency* fgreffcut;
    TEfficiency* fgreffcut_incl;
    TH1F* fheff;
    TH1F* fheff_incl;
    std::vector<float> fptbins_incl;
    int fnptbins_incl;
  };

  TH2F* createhempty(std::string name, std::string ytitle, float ymax=1.);
  TH2F* createhempty_incl(std::string name, std::string ytitle, float ymax=1.);
}

MCeff::MCefficiency::MCefficiency(std::string name, int whichincl, std::vector<float> _ptbins) : fname(name), fincl(whichincl), fptbins(_ptbins), fnptbins(_ptbins.size()-1)
{
  fheffmc = 0;
  fheffmc_incl = 0;
  fheffgen = 0;
  fheffgen_incl = 0;
  fheffgenacc = 0;
  fheffgenacc_incl = 0;
  fheffmcpre = 0;
  fheffmcpre_incl = 0;
  fheffmcbdt = 0;
  fheffmcbdt_incl = 0;
  fgreff = 0;
  fgreff_incl = 0;
  fgracc = 0;
  fgracc_incl = 0;
  fgreffpre = 0;
  fgreffpre_incl = 0;
  fgreffbdt = 0;
  fgreffbdt_incl = 0;
  fgreffqvl = 0;
  fgreffqvl_incl = 0;
  fgreffcut = 0;
  fgreffcut_incl = 0;
  fheff = 0;
  fheff_incl = 0;
  fptbins_incl = MCeff::ptBins_incl;
  fnptbins_incl = MCeff::nPtBins_incl;
  createhist();
}

MCeff::MCefficiency::MCefficiency(TFile* inf, std::string name, int whichincl) : fname(name), fincl(whichincl)
{
  fheffmc = 0;
  fheffmc_incl = 0;
  fheffgen = 0;
  fheffgen_incl = 0;
  fheffgenacc = 0;
  fheffgenacc_incl = 0;
  fheffmcpre = 0;
  fheffmcpre_incl = 0;
  fheffmcbdt = 0;
  fheffmcbdt_incl = 0;
  fgreff = 0;
  fgreff_incl = 0;
  fgracc = 0;
  fgracc_incl = 0;
  fgreffpre = 0;
  fgreffpre_incl = 0;
  fgreffbdt = 0;
  fgreffbdt_incl = 0;
  fgreffqvl = 0;
  fgreffqvl_incl = 0;
  fgreffcut = 0;
  fgreffcut_incl = 0;
  fheff = 0;
  fheff_incl = 0;
  fptbins_incl = MCeff::ptBins_incl;
  fnptbins_incl = MCeff::nPtBins_incl;
  readhist(inf);
}

void MCeff::MCefficiency::createhist()
{
  fheffmc = new TH1F(Form("heffmc%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data()); fheffmc->Sumw2();
  if(fincl==0) { fheffmc_incl = new TH1F(Form("heffmc_incl%s", fname.c_str()), "", 5, 0, 5); fheffmc_incl->Sumw2(); }
  if(fincl==1) { fheffmc_incl = new TH1F(Form("heffmc_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffmc_incl->Sumw2(); }
  fheffgen = new TH1F(Form("heffgen%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data());
  if(fincl==0) { fheffgen_incl = new TH1F(Form("heffgen_incl%s", fname.c_str()), "", 5, 0, 5); fheffgen_incl->Sumw2(); }
  if(fincl==1) { fheffgen_incl = new TH1F(Form("heffgen_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffgen_incl->Sumw2(); }
  fheffgenacc = new TH1F(Form("heffgenacc%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data());
  if(fincl==0) { fheffgenacc_incl = new TH1F(Form("heffgenacc_incl%s", fname.c_str()), "", 5, 0, 5); fheffgenacc_incl->Sumw2(); }
  if(fincl==1) { fheffgenacc_incl = new TH1F(Form("heffgenacc_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffgenacc_incl->Sumw2(); }
  fheffmcpre = new TH1F(Form("heffmcpre%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data()); fheffmcpre->Sumw2();
  if(fincl==0) { fheffmcpre_incl = new TH1F(Form("heffmcpre_incl%s", fname.c_str()), "", 5, 0, 5); fheffmcpre_incl->Sumw2(); }
  if(fincl==1) { fheffmcpre_incl = new TH1F(Form("heffmcpre_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffmcpre_incl->Sumw2(); }
  fheffmcbdt = new TH1F(Form("heffmcbdt%s", fname.c_str()), ";p_{T} (GeV/c);", fnptbins, fptbins.data()); fheffmcbdt->Sumw2();
  if(fincl==0) { fheffmcbdt_incl = new TH1F(Form("heffmcbdt_incl%s", fname.c_str()), "", 5, 0, 5); fheffmcbdt_incl->Sumw2(); }
  if(fincl==1) { fheffmcbdt_incl = new TH1F(Form("heffmcbdt_incl%s", fname.c_str()), "", fnptbins_incl, fptbins_incl.data()); fheffmcbdt_incl->Sumw2(); }
}

void MCeff::MCefficiency::readhist(TFile* inf)
{
  fheffmc          = (TH1F*)inf->Get(Form("heffmc%s", fname.c_str()));
  fheffmc_incl     = (TH1F*)inf->Get(Form("heffmc_incl%s", fname.c_str()));
  fheffgen         = (TH1F*)inf->Get(Form("heffgen%s", fname.c_str()));
  fheffgen_incl    = (TH1F*)inf->Get(Form("heffgen_incl%s", fname.c_str()));
  fheffgenacc      = (TH1F*)inf->Get(Form("heffgenacc%s", fname.c_str()));
  fheffgenacc_incl = (TH1F*)inf->Get(Form("heffgenacc_incl%s", fname.c_str()));
  fheffmcpre       = (TH1F*)inf->Get(Form("heffmcpre%s", fname.c_str()));
  fheffmcpre_incl  = (TH1F*)inf->Get(Form("heffmcpre_incl%s", fname.c_str()));
  fheffmcbdt       = (TH1F*)inf->Get(Form("heffmcbdt%s", fname.c_str()));
  fheffmcbdt_incl  = (TH1F*)inf->Get(Form("heffmcbdt_incl%s", fname.c_str()));
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

void MCeff::MCefficiency::calcacc()
{
  fgracc         = new TEfficiency(*fheffgenacc,      *fheffgen);         fgracc->SetName(Form("gracc%s",                 fname.c_str()));
  fgracc_incl    = new TEfficiency(*fheffgenacc_incl, *fheffgen_incl);    fgracc_incl->SetName(Form("gracc_incl%s",       fname.c_str()));
  fgreffpre      = new TEfficiency(*fheffmcpre,       *fheffgenacc);      fgreffpre->SetName(Form("greffpre%s",           fname.c_str()));
  fgreffpre_incl = new TEfficiency(*fheffmcpre_incl,  *fheffgenacc_incl); fgreffpre_incl->SetName(Form("greffpre_incl%s", fname.c_str()));
  fgreffcut      = new TEfficiency(*fheffmc,          *fheffmcpre);       fgreffcut->SetName(Form("greffcut%s",           fname.c_str()));
  fgreffcut_incl = new TEfficiency(*fheffmc_incl,     *fheffmcpre_incl);  fgreffcut_incl->SetName(Form("greffcut_incl%s", fname.c_str()));
  fgreffbdt      = new TEfficiency(*fheffmcbdt,       *fheffmcpre);       fgreffbdt->SetName(Form("greffbdt%s",           fname.c_str()));
  fgreffbdt_incl = new TEfficiency(*fheffmcbdt_incl,  *fheffmcpre_incl);  fgreffbdt_incl->SetName(Form("greffbdt_incl%s", fname.c_str()));
  fgreffqvl      = new TEfficiency(*fheffmc,          *fheffmcbdt);       fgreffqvl->SetName(Form("greffqvl%s",           fname.c_str()));
  fgreffqvl_incl = new TEfficiency(*fheffmc_incl,     *fheffmcbdt_incl);  fgreffqvl_incl->SetName(Form("greffqvl_incl%s", fname.c_str()));
}

void MCeff::MCefficiency::setstyle(Color_t color, Style_t mstyle/*=20*/, Style_t lstyle/*=2*/, Width_t lwidth/*=3*/)
{
  if(fheff)       { xjjroot::setthgrstyle(fheff      , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fheff_incl)  { xjjroot::setthgrstyle(fheff_incl , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreff)      { xjjroot::setthgrstyle(fgreff     , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreff_incl) { xjjroot::setthgrstyle(fgreff_incl, color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgracc)      { xjjroot::setthgrstyle(fgracc     , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgracc_incl) { xjjroot::setthgrstyle(fgracc_incl, color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffpre)      { xjjroot::setthgrstyle(fgreffpre     , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffpre_incl) { xjjroot::setthgrstyle(fgreffpre_incl, color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffcut)      { xjjroot::setthgrstyle(fgreffcut     , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffcut_incl) { xjjroot::setthgrstyle(fgreffcut_incl, color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffbdt)      { xjjroot::setthgrstyle(fgreffbdt     , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffbdt_incl) { xjjroot::setthgrstyle(fgreffbdt_incl, color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffqvl)      { xjjroot::setthgrstyle(fgreffqvl     , color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
  if(fgreffqvl_incl) { xjjroot::setthgrstyle(fgreffqvl_incl, color, mstyle, 1.2, color, lstyle, lwidth, color, 0.2, 1001); }
}

TH2F* MCeff::createhempty(std::string name, std::string ytitle, float ymax)
{
  TH2F* hemptyeff = new TH2F(name.c_str(), Form(";p_{T} (GeV/c);%s", ytitle.c_str()), 10, MCeff::ptBins[0], MCeff::ptBins[MCeff::nPtBins], 10, 0, ymax);
  xjjroot::sethempty(hemptyeff, 0, 0.3);
  return hemptyeff;
}

TH2F* MCeff::createhempty_incl(std::string name, std::string ytitle, float ymax)
{
  TH2F* hemptyeff_incl = new TH2F(name.c_str(), Form(";;%s", ytitle.c_str()), 5, 0, 5, 10, 0, ymax);
  xjjroot::sethempty(hemptyeff_incl, 0, 0.3);
  hemptyeff_incl->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hemptyeff_incl->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hemptyeff_incl->GetXaxis()->SetLabelSize(hemptyeff_incl->GetXaxis()->GetLabelSize()*1.5);
  return hemptyeff_incl;
}
