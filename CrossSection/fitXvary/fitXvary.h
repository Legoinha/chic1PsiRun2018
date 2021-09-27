#ifndef __FITX_FITXVARY_H_
#define __FITX_FITXVARY_H_

#include <vector>
#include <string>
#include <TH1F.h>
#include <TEfficiency.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include "packtree.h"
#include "ntuple.h"
#include "xjjrootuti.h"
#include "xjjcuti.h"

#include "fit.h"

namespace fitX
{
  float Qvaluecut = 0.13;
  class varymva
  {
  public:
    varymva() { ; }
    int n() { return fmva.size(); }
    std::string type() { return ftype; }
    float thatval() { return fthatval; }
    float minval() { return fminval; }
    std::vector<float> mva() { return fmva; }
    std::vector<bool> ifdraw() { return fifdraw; }
    float binmin() { return fbinmin; }
    float binmax() { return fbinmax; }
    void set(std::string _type, std::vector<float> _mva, std::vector<bool> _ifdraw, float _thatval, float _minval) 
    {
      ftype = _type;
      fmva = _mva;
      fifdraw = _ifdraw;
      fthatval = _thatval;
      fminval = _minval;
      fbinmax = fmva.back();
      float delta = fmva[fmva.size()-1] - fmva[fmva.size()-2];
      fbinmin = fmva.front()-delta;
    }
  private:
    std::string ftype;
    std::vector<float> fmva;
    std::vector<bool> fifdraw;
    float fbinmin, fbinmax;
    float fthatval;
    float fminval;
  };
}

namespace fitX
{
  fitX::varymva* initvarycut(std::string type);
}

namespace fitX
{
  class varycut
  {
  public:
    std::vector<TH1F*> hdata, hdataBenr; // mass
    std::vector<TH1F*> hmc_a, hmc_b; // MC mass
    std::vector<RooDataSet*> dshdata, dshdataBenr; // mass
    std::vector<RooDataSet*> dshmc_a, dshmc_b; // MC mass
    std::vector<TH1F*> hlxymcnp_a, hlxymcnp_b, hlxymcp_a, hlxymcp_b; // MC np/p lxy dis
    TH1F *heffmc_a, *heffmc_b;
    TH1F *heffgen_a, *heffgen_b;
    TH1F *heff_a, *heff_b;
    TEfficiency *greff_a, *greff_b;
    TH1F *hsideband_a, *hsideband_b; // sideband
    varycut(std::vector<float> vvary, std::string titlevary) : fvv(vvary), fvt(titlevary) 
    {
      fnv = vvary.size();
      hdata.resize(fnv);
      hmc_a.resize(fnv);
      hmc_b.resize(fnv);
      hdataBenr.resize(fnv);
      dshdata.resize(fnv);
      dshmc_a.resize(fnv);
      dshmc_b.resize(fnv);
      dshdataBenr.resize(fnv);
      hlxymcnp_a.resize(fnv);
      hlxymcnp_b.resize(fnv);
      hlxymcp_a.resize(fnv);
      hlxymcp_b.resize(fnv);
      mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
      massBenr = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
      massmc_a = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_L, fitX::BIN_MAX_L);
      massmc_b = new RooRealVar("Bmass", "massmc_b", fitX::BIN_MIN_H, fitX::BIN_MAX_H);
      pthatweight = new RooRealVar("pthatweight", "pthatweight", 0., 10.); //
    }
    const float masswinL = 0.03, masswinH = 0.06;

    int producehist();
    int loop_vary_data(xjjroot::packtree* pt, int nentries);
    int loop_vary_mcprompt(xjjroot::packtree* pt, int nentries, std::string name);
    int loop_vary_mcnonprompt(xjjroot::packtree* pt, int nentries, std::string name);

    int produceeff();

  private:
    RooRealVar *mass, *massBenr, *massmc_a, *massmc_b;
    RooRealVar *pthatweight;
    float MVA(mytmva::ntuple* ntp, int j);
    int fnv;
    std::vector<float> fvv;
    std::string fvt;
    std::vector<TH1F*> fhmc, fhlxymcnp, fhlxymcp;
    std::vector<RooDataSet*> fdshmc;
    TH1F *fheffmc, *fheffgen, *fheff;
    TEfficiency *fgreff;
    RooRealVar* fmassmc;
    int setaorb(std::string name) 
    {
      if(xjjc::str_contains(name, "_a"))
        {
          fdshmc = dshmc_a;
          fmassmc = massmc_a;
          fhmc = hmc_a;
          fheffmc = heffmc_a;
          fheffgen = heffgen_a;
          fheff = heff_a;
          fgreff = greff_a;
          fhlxymcnp = hlxymcnp_a;
          fhlxymcp = hlxymcp_a;
          return 0;
        }
      if(xjjc::str_contains(name, "_b")) 
        {
          fdshmc = dshmc_b;
          fmassmc = massmc_b;
          fhmc = hmc_b;
          fheffmc = heffmc_b;
          fheffgen = heffgen_b;
          fheff = heff_b;
          fgreff = greff_b;
          fhlxymcnp = hlxymcnp_b;
          fhlxymcp = hlxymcp_b;
          return 1;
        }
      return -1;
    }
  };
}

float fitX::varycut::MVA(mytmva::ntuple* ntp, int j)
{
  if(fvt=="BDTG") { return ntp->BDTG[j]; }
  if(fvt=="BDT") { return ntp->BDT[j]; }
  if(fvt=="BDTD") { return ntp->BDTD[j]; }
  if(fvt=="BDTF") { return ntp->BDTF[j]; }
  return -10;
}

int fitX::varycut::loop_vary_data(xjjroot::packtree* pt, int nentries)
{
  mytmva::ntuple* ntp = pt->ntp();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!ntp->passedevtfil()) continue;
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          
          if(!(ntp->Bpt[j] > fitX::ptmincut && ntp->Bpt[j] < fitX::ptmaxcut && 
               TMath::Abs(ntp->By[j]) >= fitX::ymincut && TMath::Abs(ntp->By[j]) < fitX::ymaxcut && 
               ntp->hiBin >= fitX::centmincut && ntp->hiBin <= fitX::centmaxcut)) continue;
          if(!((ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < fitX::Qvaluecut)) continue; // !!
          if(ntp->Bmass[j] < fitX::BIN_MIN || ntp->Bmass[j] > fitX::BIN_MAX) continue;
          mass->setVal(ntp->Bmass[j]);
          massBenr->setVal(ntp->Bmass[j]);
          for(int l=0; l<fnv; l++)
            {
              if(MVA(ntp, j) > fvv[l])
                {
                  hdata[l]->Fill(ntp->Bmass[j]);
                  dshdata[l]->add(*mass);
                  if(TMath::Abs(ntp->Bmass[j]-MASS_PSI2S) > masswinL && TMath::Abs(ntp->Bmass[j]-MASS_PSI2S) < masswinH) { hsideband_a->Fill(hsideband_a->GetBinCenter(l+1)); }
                  if(TMath::Abs(ntp->Bmass[j]-MASS_X) > masswinL && TMath::Abs(ntp->Bmass[j]-MASS_X) < masswinH) { hsideband_b->Fill(hsideband_b->GetBinCenter(l+1)); }
                  if((10*ntp->Blxy[j]*ntp->Bmass[j]/ntp->Bpt[j]) > 0.1)
                    {
                      hdataBenr[l]->Fill(ntp->Bmass[j]);
                      dshdataBenr[l]->add(*massBenr);
                    }
                }
              else break;
            }
        }
    }
  xjjc::progressbar_summary(nentries);
}

int fitX::varycut::loop_vary_mcprompt(xjjroot::packtree* pt, int nentries, std::string name)
{
  int isam = setaorb(name);
  if(isam < 0) return 1;
  mytmva::ntuple* ntp = pt->ntp();

  //
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!(ntp->hiBin >=0 && ntp->hiBin<180)) continue;
      float weight = ntp->pthatweight * ntp->Ncoll;

      for(int j=0; j<ntp->Gsize; j++)
        {
          if(!(ntp->GisSignal[j]==7 && ntp->GcollisionId[j]==0)) continue;
          if(!(ntp->Gpt[j] > fitX::ptmincut && ntp->Gpt[j] < fitX::ptmaxcut && TMath::Abs(ntp->Gy[j]) >= fitX::ymincut && TMath::Abs(ntp->Gy[j]) < fitX::ymaxcut && ntp->hiBin >= fitX::centmincut && ntp->hiBin <= fitX::centmaxcut)) continue;
          for(int l=0; l<fnv; l++)
            {
              fheffgen->Fill(fheffgen->GetBinCenter(l+1), weight);
            }
        }

      if(!ntp->passedevtfil()) continue;
      pthatweight->setVal(ntp->pthatweight);
      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;

          if(!(ntp->Bpt[j] > fitX::ptmincut && ntp->Bpt[j] < fitX::ptmaxcut && TMath::Abs(ntp->By[j]) >= fitX::ymincut && TMath::Abs(ntp->By[j]) < fitX::ymaxcut && ntp->hiBin >= fitX::centmincut && ntp->hiBin <= fitX::centmaxcut)) continue;
          if(!((ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < fitX::Qvaluecut)) continue;
          for(int l=0; l<fnv; l++)
            {
              if(MVA(ntp, j) > fvv[l])
                {
                  fhmc[l]->Fill(ntp->Bmass[j], ntp->pthatweight);
                  fmassmc->setVal(ntp->Bmass[j]);
                  fdshmc[l]->add(*fmassmc, ntp->pthatweight);
                  fheffmc->Fill(fheffmc->GetBinCenter(l+1), weight); // weight!
                  fhlxymcp[l]->Fill((10*ntp->Blxy[j]*ntp->Bmass[j]/ntp->Bpt[j]), weight); // weight!
                }
              else break;
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  return 0;
}
  
int fitX::varycut::loop_vary_mcnonprompt(xjjroot::packtree* pt, int nentries, std::string name)
{
  int isam = setaorb(name);
  if(isam < 0) return 1;
  mytmva::ntuple* ntp = pt->ntp();

  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      if(!(ntp->hiBin >=0 && ntp->hiBin<180)) continue;
      float weight = ntp->pthatweight * ntp->Ncoll;
      if(!ntp->passedevtfil()) continue;

      for(int j=0; j<ntp->Bsize; j++)
        {
          if(!ntp->mvapref[j]) continue;
          if(!(ntp->Bpt[j] > fitX::ptmincut && ntp->Bpt[j] < fitX::ptmaxcut && TMath::Abs(ntp->By[j]) >= fitX::ymincut && TMath::Abs(ntp->By[j]) < fitX::ymaxcut && ntp->hiBin >= fitX::centmincut && ntp->hiBin <= fitX::centmaxcut)) continue;
          if(!(ntp->Bgen[j]>=23333 && ntp->BgencollisionId[j]==0)) continue;
          if(!((ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < fitX::Qvaluecut)) continue;
          for(int l=0; l<fnv; l++)
            {
              if(MVA(ntp, j) > fvv[l])
                {
                  fhlxymcnp[l]->Fill((10*ntp->Blxy[j]*ntp->Bmass[j]/ntp->Bpt[j]), weight); // weight!
                }
              else break;
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  return 0;
}

int fitX::varycut::producehist()
{
  for(int l=0; l<fnv; l++)
    {
      hdata[l] = new TH1F(Form("hdata_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
      hdata[l]->Sumw2();
      hdataBenr[l] = new TH1F(Form("hdataBenr_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
      hdataBenr[l]->Sumw2();
      hmc_a[l] = new TH1F(Form("hmc_a_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L);
      hmc_a[l]->Sumw2();
      hmc_b[l] = new TH1F(Form("hmc_b_%d", l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H);
      hmc_b[l]->Sumw2();

      dshdata[l] = new RooDataSet(Form("dshdata_%d", l), "", RooArgSet(*mass));
      dshdataBenr[l] = new RooDataSet(Form("dshdataBenr_%d", l), "", RooArgSet(*massBenr));
      dshmc_a[l] = new RooDataSet(Form("dshmc_a_%d", l), "", RooArgSet(*massmc_a, *pthatweight));
      dshmc_b[l] = new RooDataSet(Form("dshmc_b_%d", l), "", RooArgSet(*massmc_b, *pthatweight));
    }
  heffmc_a = new TH1F("heffmc_a", Form(";%s;#alpha #times #epsilon", fvt.c_str()), fvv.size()-1, fvv.data());
  heffmc_b = new TH1F("heffmc_b", Form(";%s;#alpha #times #epsilon", fvt.c_str()), fvv.size()-1, fvv.data());
  heffgen_a = new TH1F("heffgen_a", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  heffgen_b = new TH1F("heffgen_b", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hsideband_a = new TH1F("hsideband_a", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
  hsideband_b = new TH1F("hsideband_b", Form(";%s;", fvt.c_str()), fvv.size()-1, fvv.data());
}

int fitX::varycut::produceeff()
{
  heff_a = (TH1F*)heffmc_a->Clone("heff_a");
  heff_b = (TH1F*)heffmc_b->Clone("heff_b");
  greff_a = new TEfficiency(*heffmc_a, *heffgen_a); greff_a->SetName("greff_a");
  greff_b = new TEfficiency(*heffmc_b, *heffgen_b); greff_b->SetName("greff_b");
  for(int i=0;i<heff_a->GetNbinsX();i++)
    {
      heff_a->SetBinContent(i+1, greff_a->GetEfficiency(i+1));
      heff_a->SetBinError(i+1, greff_a->GetEfficiencyErrorUp(i+1));
      heff_b->SetBinContent(i+1, greff_b->GetEfficiency(i+1));
      heff_b->SetBinError(i+1, greff_b->GetEfficiencyErrorUp(i+1));
    }
  return 0;
}

void drawalltext()
{
  xjjroot::drawCMSleft("Internal");
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.04, "X(3872)", 0.038, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.90, 0.84, fitX::pttag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 62);
}

void drawvline(TH2F* hempty)
{
  xjjroot::drawline(0.06, hempty->GetYaxis()->GetXmin(), 0.06, hempty->GetYaxis()->GetXmax(), kGray+1, 2, 2);
  xjjroot::drawline(0.04, hempty->GetYaxis()->GetXmin(), 0.04, hempty->GetYaxis()->GetXmax(), kGray+1, 2, 2);
}

void drawvlineonly(TH2F* hempty)
{
  xjjroot::drawline(0.06, hempty->GetYaxis()->GetXmin(), 0.06, hempty->GetYaxis()->GetXmax(), kGray+1, 2, 2);
  // xjjroot::drawline(0.04, hempty->GetYaxis()->GetXmin(), 0.04, hempty->GetYaxis()->GetXmax(), kGray+1, 2, 2);
}

fitX::varymva* fitX::initvarycut(std::string type)
{
  fitX::varymva* varbin = new fitX::varymva();
  if(type=="BDTG")
    {
      varbin->set("BDTG",
                  std::vector<float>({-1 , -0.8 , -0.6 , -0.4 , -0.2 , -0.1 , 0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.55 , 0.6 , 0.65 , 0.70 , 0.72, 0.75 , 0.80, 0.85 , 1.0}),
                  std::vector<bool>({false, false, false, false, false, false, true, false, true , false, true, false, false, true , false, true, false, true , true, false, false}),
                  0.7,
                  0.5
                  );
    }
  else if(type=="BDT")
    {
      varbin->set("BDT",
                  std::vector<float>({0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  0.1,  0.11,  0.12, 0.13}),
                  std::vector<bool>({true,  true, true, true, true, false, true, false, true, false}),
                  0.06,
                  0
                  );
      // varbin->set("BDT",
      //             std::vector<float>({-0.2, -0.1 , -0.04, -0.02, 0    , 0.02 , 0.04, 0.05 , 0.06, 0.07, 0.08, 0.09 ,  0.1,  0.11,  0.12, 0.13 , 0.15 , 0.20 , 0.30}) ,
      //             std::vector<bool>({false, false, false, false, false, false, true,  true, true, true, true, false, true, false, true , false, false, false, false}),
      //             0.06,
      //             0
      //             );
    }
  else if(type=="BDTF")
    {
      varbin->set("BDTF",
                  std::vector<float>({-0.5, 0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4}),
                  std::vector<bool>({true, true, true, true, true, true, true, true}),
                  0.3,
                  0.1
                  );
    }
  else if(type=="BDTD")
    {
      varbin->set("BDTD",
                  std::vector<float>({-0.5, 0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4}),
                  std::vector<bool>({true, true, true, true, true, true, true, true}),
                  0.12,
                  0.1
                  );
    }
  else { return 0; }
  return varbin;
}


#endif
