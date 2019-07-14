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

void fitXvary(std::string inputdata,                            // ==> data
              std::string inputmc_a, std::string inputmc_b,     // ==> prompt MC
              std::string inputmcnp_a, std::string inputmcnp_b, // ==> nonprompt MC
              std::string output)
{
  // init
  std::map<std::string, std::vector<float>> xbins = lxydis::setupbins();

  fitX::varycut vc(bdtg, "BDTG");
  vc.producehist();
  for(int l=0; l<nbdtg; l++)
    {
      vc.hlxymcnp_a[l] = new TH1F(Form("hlxymcnp_a_%d", l), ";l_{xy} (mm);Probability", xbins["lxynonprompt"].size()-1, xbins["lxynonprompt"].data());
      vc.hlxymcnp_a[l]->Sumw2();
      vc.hlxymcnp_b[l] = new TH1F(Form("hlxymcnp_b_%d", l), ";l_{xy} (mm);Probability", xbins["lxynonprompt"].size()-1, xbins["lxynonprompt"].data());
      vc.hlxymcnp_b[l]->Sumw2();
      vc.hlxymcp_a[l] = new TH1F(Form("hlxymcp_a_%d", l), ";l_{xy} (mm);Probability", xbins["lxyprompt"].size()-1, xbins["lxyprompt"].data());
      vc.hlxymcp_a[l]->Sumw2();
      vc.hlxymcp_b[l] = new TH1F(Form("hlxymcp_b_%d", l), ";l_{xy} (mm);Probability", xbins["lxyprompt"].size()-1, xbins["lxyprompt"].data());
      vc.hlxymcp_b[l]->Sumw2();
    }

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
  // mytmva::ntuple* ntp;
  int nentries;
  
  // ------------------------------
  int ss = 0;
  // ==> data <==
  std::cout<<std::endl<<" ==> ("<<++ss<<"/5) ==> "<<"Processing data"<<std::endl;
  pt = new xjjroot::packtree(infdata, "Bfinder/ntmix", "data");
  nentries = pt->getentries();
  vc.loop_vary_data(pt, nentries);

  // ==> Prompt psi' MC <==
  std::cout<<std::endl<<" ==> ("<<++ss<<"/5) ==> "<<"Processing prompt psi' MC"<<std::endl;
  pt = new xjjroot::packtree(infmc_a, "Bfinder/ntmix", "mc_a", "Bfinder/ntGen");
  nentries = pt->getentries();
  vc.loop_vary_mcprompt(pt, nentries, "mc_a");

  // ==> Prompt X MC <==
  std::cout<<std::endl<<" ==> ("<<++ss<<"/5) ==> "<<"Processing prompt X MC"<<std::endl;
  pt = new xjjroot::packtree(infmc_b, "Bfinder/ntmix", "mc_b", "Bfinder/ntGen");
  nentries = pt->getentries();
  vc.loop_vary_mcprompt(pt, nentries, "mc_b");

  // ==> Nonprompt psi' MC <==
  std::cout<<std::endl<<" ==> ("<<++ss<<"/5) ==> "<<"Processing nonprompt psi' MC"<<std::endl;
  pt = new xjjroot::packtree(infmcnp_a, "Bfinder/ntmix", "mcnp_a");
  nentries = pt->getentries();
  vc.loop_vary_mcnonprompt(pt, nentries, "mcnp_a");

  // ==> Nonprompt X MC <==
  std::cout<<std::endl<<" ==> ("<<++ss<<"/5) ==> "<<"Processing nonprompt X MC"<<std::endl;
  pt = new xjjroot::packtree(infmcnp_b, "Bfinder/ntmix", "mcnp_b");
  nentries = pt->getentries();
  vc.loop_vary_mcnonprompt(pt, nentries, "mcnp_b");

  // ------------------------------
  std::cout<<std::endl;
  vc.produceeff();

  //
  TFile* outf = new TFile(Form("rootfiles/root_fitXvary_%s.root", output.c_str()), "recreate");
  outf->cd();
  for(auto& hh : vc.hdata) { hh->Write(); }
  for(auto& hh : vc.hdataBenr) { hh->Write(); }
  for(auto& hh : vc.hmc_a) { hh->Write(); }
  for(auto& hh : vc.hmc_b) { hh->Write(); }
  vc.heff_a->Write();
  vc.heff_b->Write();
  vc.greff_a->Write();
  vc.greff_b->Write();
  for(auto& hh : vc.hlxymcnp_a) { hh->Write(); }
  for(auto& hh : vc.hlxymcnp_b) { hh->Write(); }
  for(auto& hh : vc.hlxymcp_a) { hh->Write(); }
  for(auto& hh : vc.hlxymcp_b) { hh->Write(); }
  vc.hsideband_a->Write();
  vc.hsideband_b->Write();
  outf->Close();
  
}

int main(int argc, char* argv[])
{
  if(argc==7) { fitXvary(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); return 0; }
  return 1;
}
