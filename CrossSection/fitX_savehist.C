#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

#include "fitX.h"
#include "project.h"
#include "lxydis.h"
#include "MCefficiency.h"
#include "xjjcuti.h"

float weight_ss = 1./0.681;

void fitX_savehist(std::string input, std::string inputmcp_a, std::string inputmcp_b, std::string inputmcnp_a, std::string inputmcnp_b,
                   std::string cut, std::string cutgen, std::string output)
{

  std::cout<<cut<<std::endl;
  std::map<std::string, std::vector<float>> lxyxbins = lxydis::setupbins();
  
  std::string mcweight = "(pthatweight*Ncoll)";

  //
  TTree* ntmix = (TTree*)fitX::getnt(input, "ntmix"); if(!ntmix) { return; }
  TTree* ntmixmcp_a = fitX::getnt(inputmcp_a, "ntmix"); if(!ntmixmcp_a) { return; }
  TTree* ntmixmcp_b = fitX::getnt(inputmcp_b, "ntmix"); if(!ntmixmcp_b) { return; }
  TTree* ntGenmcp_a = fitX::getnt(inputmcp_a, "ntGen"); if(!ntGenmcp_a) { return; }
  TTree* ntGenmcp_b = fitX::getnt(inputmcp_b, "ntGen"); if(!ntGenmcp_b) { return; }
  TTree* ntmixmcnp_a = fitX::getnt(inputmcnp_a, "ntmix"); if(!ntmixmcnp_a) { return; }
  TTree* ntmixmcnp_b = fitX::getnt(inputmcnp_b, "ntmix"); if(!ntmixmcnp_b) { return; }

  // TH1 must be defined after TTree declaration (some tricky issue)
  TH1F* h = new TH1F("h", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h->Sumw2();
  TH1F* hBenr = new TH1F("hBenr", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); hBenr->Sumw2();
  TH1F* hmcp_a = new TH1F("hmcp_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); hmcp_a->Sumw2();
  TH1F* hmcp_b = new TH1F("hmcp_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); hmcp_b->Sumw2();
  TH1F* hlxymcnp_a = new TH1F("hlxymcnp_a", ";l_{xy} (mm);Probability", lxyxbins["lxynonprompt"].size()-1, lxyxbins["lxynonprompt"].data()); hlxymcnp_a->Sumw2();
  TH1F* hlxymcnp_b = new TH1F("hlxymcnp_b", ";l_{xy} (mm);Probability", lxyxbins["lxynonprompt"].size()-1, lxyxbins["lxynonprompt"].data()); hlxymcnp_b->Sumw2();
  TH1F* hlxymcp_a = new TH1F("hlxymcp_a", ";l_{xy} (mm);Probability", lxyxbins["lxyprompt"].size()-1, lxyxbins["lxyprompt"].data()); hlxymcp_a->Sumw2();
  TH1F* hlxymcp_b = new TH1F("hlxymcp_b", ";l_{xy} (mm);Probability", lxyxbins["lxyprompt"].size()-1, lxyxbins["lxyprompt"].data()); hlxymcp_b->Sumw2();
  MCeff::MCefficiency mceff_a("_a");
  MCeff::MCefficiency mceff_b("_b");

  std::string cutreco = Form("(%s) && Bpt>%f && Bpt<%f && TMath::Abs(By)<%f", cut.c_str(), fitX::ptmincut, fitX::ptmaxcut, fitX::ycut);
  std::string cutmcreco = Form("%s && Bgen>=23333 && BgencollisionId==0", cutreco.c_str());
  std::string cutmcgen = Form("(%s) && Gpt>%f && Gpt<%f && TMath::Abs(Gy)<%f && GisSignal==7 && GcollisionId==0", cutgen.c_str(), fitX::ptmincut, fitX::ptmaxcut, fitX::ycut);

  std::cout<<" == data ==>"<<std::endl;
  ntmix->Project("h", "Bmass", TCut(cutreco.c_str()));
  std::cout<<h->GetEntries()<<std::endl;
  ntmix->Project("hBenr", "Bmass", TCut(Form("%s && Blxy > 0.1", cutreco.c_str())));
  std::cout<<hBenr->GetEntries()<<std::endl;

  std::cout<<" == mcp_a ==>"<<std::endl;
  ntmixmcp_a->Project("hmcp_a", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // mass shape weight
  std::cout<<hmcp_a->GetEntries()<<std::endl;
  hmcp_a->Scale(hmcp_a->GetEntries()/hmcp_a->Integral());
  ntmixmcp_a->Project(mceff_a.heffmc->GetName(), "Bpt", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp_a->Project(mceff_a.heffmc_incl->GetName(), Form("%d", fitX::ibin_a-1), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp_a->Project("hlxymcp_a", "Blxy", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));

  std::cout<<" == mcp_b ==>"<<std::endl;
  ntmixmcp_b->Project("hmcp_b", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // mass shape weight
  std::cout<<hmcp_b->GetEntries()<<std::endl;
  hmcp_b->Scale(hmcp_b->GetEntries()/hmcp_b->Integral());
  ntmixmcp_b->Project(mceff_b.heffmc->GetName(), "Bpt", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp_b->Project(mceff_b.heffmc_incl->GetName(), Form("%d", fitX::ibin_b-1), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp_b->Project("hlxymcp_b", "Blxy", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));

  std::cout<<" == mcgenp_a ==>"<<std::endl;
  ntGenmcp_a->Project(mceff_a.heffgen->GetName(), "Gpt", TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));
  ntGenmcp_a->Project(mceff_a.heffgen_incl->GetName(), Form("%d", fitX::ibin_a-1), TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));

  std::cout<<" == mcgenp_b ==>"<<std::endl;
  ntGenmcp_b->Project(mceff_b.heffgen->GetName(), "Gpt", TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));
  ntGenmcp_b->Project(mceff_b.heffgen_incl->GetName(), Form("%d", fitX::ibin_b-1), TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));

  std::cout<<" == mcnp_a ==>"<<std::endl;
  ntmixmcnp_a->Project("hlxymcnp_a", "Blxy", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  std::cout<<hlxymcnp_a->GetEntries()<<std::endl;

  std::cout<<" == mcnp_b ==>"<<std::endl;
  ntmixmcnp_b->Project("hlxymcnp_b", "Blxy", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  std::cout<<hlxymcnp_b->GetEntries()<<std::endl;

  TFile* outf = new TFile(Form("%s.root", output.c_str()), "recreate");
  outf->cd();
  h->Write();
  hBenr->Write();
  hmcp_a->Write();
  hmcp_b->Write();
  hlxymcnp_a->Write();
  hlxymcnp_b->Write();
  hlxymcp_a->Write();
  hlxymcp_b->Write();
  mceff_a.heffmc->Write();
  mceff_a.heffmc_incl->Write();
  mceff_a.heffgen->Write();
  mceff_a.heffgen_incl->Write();
  mceff_b.heffmc->Write();
  mceff_b.heffmc_incl->Write();
  mceff_b.heffgen->Write();
  mceff_b.heffgen_incl->Write();
  outf->cd();
  TTree* info = new TTree("info", "cut info");
  info->Branch("input", &input);
  info->Branch("inputmcp_a", &inputmcp_a);
  info->Branch("inputmcp_b", &inputmcp_b);
  info->Branch("inputmcnp_a", &inputmcnp_a);
  info->Branch("inputmcnp_b", &inputmcnp_b);
  info->Branch("cutreco", &cutreco);
  info->Branch("cutmcreco", &cutmcreco);
  info->Branch("cutmcgen", &cutmcgen);
  info->Fill();
  info->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==9) { fitX_savehist(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]); return 0; }
  return 1;
}
