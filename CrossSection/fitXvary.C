#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <string>

#include "packtree.h"
#include "ntuple.h"

#include "fitXvary.h"
#include "lxydis.h"

void fitXvary(std::string inputdata, std::string inputmc_a, std::string inputmc_b, std::string inputmcnp_a, std::string inputmcnp_b,
              std::string output)
{
  // init
  std::map<std::string, std::vector<float>> xbins = lxydis::setupbins();
  
  for(int k=0; k<ndls; k++)
    {
      for(int l=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hdata[idx] = new TH1F(Form("hdata_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); 
          hdata[idx]->Sumw2();
          hdatagt[idx] = new TH1F(Form("hdatagt_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); 
          hdatagt[idx]->Sumw2();
          hmc_a[idx] = new TH1F(Form("hmc_a_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); 
          hmc_a[idx]->Sumw2();
          hmc_b[idx] = new TH1F(Form("hmc_b_%d_%d", k, l), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); 
          hmc_b[idx]->Sumw2();
          hlxymcnp_a[idx] = new TH1F(Form("hlxymcnp_a_%d_%d", k, l), ";l_{xy} (mm);Probability", xbins["lxynonprompt"].size()-1, xbins["lxynonprompt"].data());
          hlxymcnp_a[idx]->Sumw2();
          hlxymcnp_b[idx] = new TH1F(Form("hlxymcnp_b_%d_%d", k, l), ";l_{xy} (mm);Probability", xbins["lxynonprompt"].size()-1, xbins["lxynonprompt"].data());
          hlxymcnp_b[idx]->Sumw2();
        }
      heffmc_a[k] = new TH1F(Form("heffmc_a_%d", k), ";BDTG;#alpha #times #epsilon", bdtg.size()-1, bdtg.data());
      heffmc_b[k] = new TH1F(Form("heffmc_b_%d", k), ";BDTG;#alpha #times #epsilon", bdtg.size()-1, bdtg.data());
    }
  heffgen_a = new TH1F("heffgen_a", ";BDTG;", bdtg.size()-1, bdtg.data());
  heffgen_b = new TH1F("heffgen_b", ";BDTG;", bdtg.size()-1, bdtg.data());

  // open files
  TFile* infdata = TFile::Open(inputdata.c_str());
  if(!infdata->IsOpen()) return;
  TFile* infmc_a = TFile::Open(inputmc_a.c_str());
  if(!infmc_a->IsOpen()) return;
  TFile* infmc_b = TFile::Open(inputmc_b.c_str());
  if(!infmc_b->IsOpen()) return;
  TFile* infmcnp_a = TFile::Open(inputmcnp_a.c_str());
  if(!infmcnp_a->IsOpen()) return;
  TFile* infmcnp_b = TFile::Open(inputmcnp_b.c_str());
  if(!infmcnp_b->IsOpen()) return;

  xjjroot::packtree* pt;
  mytmva::ntuple* ntp;
  int nentries;
  // ------------------------------

  // >>>> data <<<<
  pt = new xjjroot::packtree(infdata, "Bfinder/ntmix", "data");
  ntp = pt->ntp;
  nentries = pt->getentries();
  fitX::loop_vary_data(pt, ntp, nentries);

  // >>>> Prompt psi' MC <<<<
  pt = new xjjroot::packtree(infmc_a, "Bfinder/ntmix", "mc_a", "Bfinder/ntGen");
  ntp = pt->ntp;
  nentries = pt->getentries();
  fitX::loop_vary_mcprompt(pt, ntp, nentries, "mc_a");

  // >>>> Prompt X MC <<<<
  pt = new xjjroot::packtree(infmc_b, "Bfinder/ntmix", "mc_b", "Bfinder/ntGen");
  ntp = pt->ntp;
  nentries = pt->getentries();
  fitX::loop_vary_mcprompt(pt, ntp, nentries, "mc_b");

  // >>>> Nonprompt psi' MC <<<<
  pt = new xjjroot::packtree(infmcnp_a, "Bfinder/ntmix", "mcnp_a");
  ntp = pt->ntp;
  nentries = pt->getentries();
  fitX::loop_vary_mcnonprompt(pt, ntp, nentries, "mcnp_a");

  // >>>> Nonprompt X MC <<<<
  pt = new xjjroot::packtree(infmcnp_b, "Bfinder/ntmix", "mcnp_b");
  ntp = pt->ntp;
  nentries = pt->getentries();
  fitX::loop_vary_mcnonprompt(pt, ntp, nentries, "mcnp_b");

  // ------------------------------
  for(int k=0; k<ndls; k++)
    {
      heff_a[k] = (TH1F*)heffmc_a[k]->Clone(Form("heff_a_%d", k));
      heff_b[k] = (TH1F*)heffmc_b[k]->Clone(Form("heff_b_%d", k));
      greff_a[k] = new TEfficiency(*(heffmc_a[k]), *heffgen_a); greff_a[k]->SetName(Form("greff_a_%d", k));
      greff_b[k] = new TEfficiency(*(heffmc_b[k]), *heffgen_b); greff_b[k]->SetName(Form("greff_b_%d", k));
      for(int i=0;i<heff_a[k]->GetNbinsX();i++)
        {
          heff_a[k]->SetBinContent(i+1, greff_a[k]->GetEfficiency(i+1));
          heff_a[k]->SetBinError(i+1, greff_a[k]->GetEfficiencyErrorUp(i+1));
          heff_b[k]->SetBinContent(i+1, greff_b[k]->GetEfficiency(i+1));
          heff_b[k]->SetBinError(i+1, greff_b[k]->GetEfficiencyErrorUp(i+1));
        }
    }  

  //
  TFile* outf = new TFile(Form("rootfiles/root_fitXvary_%s.root", output.c_str()), "recreate");
  outf->cd();
  for(auto& hh : hdata) { hh->Write(); }
  for(auto& hh : hdatagt) { hh->Write(); }
  for(auto& hh : hmc_a) { hh->Write(); }
  for(auto& hh : hmc_b) { hh->Write(); }
  for(auto& hh : heff_a) { hh->Write(); }
  for(auto& hh : heff_b) { hh->Write(); }
  for(auto& hh : greff_a) { hh->Write(); }
  for(auto& hh : greff_b) { hh->Write(); }
  for(auto& hh : hlxymcnp_a) { hh->Write(); }
  for(auto& hh : hlxymcnp_b) { hh->Write(); }
  outf->Close();
  
}

int main(int argc, char* argv[])
{
  if(argc==7) { fitXvary(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); return 0; }
  return 1;
}
