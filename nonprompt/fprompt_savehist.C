#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "fit.h"
#include "project.h"
#include "lxydis.h"
#include "xjjcuti.h"

void fprompt_savehist(std::string input, std::string inputmcp_a, std::string inputmcp_b, std::string inputmcnp_a, std::string inputmcnp_b,
                      std::string cut, std::string output, std::string lxyvar="lxy")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::cout<<cut<<std::endl;
  std::map<std::string, std::vector<float>> lxyxbins = lxydis::setupbins();
  std::string input_flatten = xjjc::str_replaceall(input, ".root", "_flatten.root");
  std::string inputmcp_a_flatten = xjjc::str_replaceall(inputmcp_a, ".root", "_flatten.root");
  std::string inputmcp_b_flatten = xjjc::str_replaceall(inputmcp_b, ".root", "_flatten.root");

  std::string mcweight = "(pthatweight*Ncoll)";
  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooRealVar* massmc_a = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_L, fitX::BIN_MAX_L);
  RooRealVar* massmc_b = new RooRealVar("Bmass", "massmc_b", fitX::BIN_MIN_H, fitX::BIN_MAX_H);
  RooRealVar* pthatweight = new RooRealVar("pthatweight", "pthatweight", 0, 10.); // pthatweight range!!

  // TH1 must be defined after TTree declaration (some tricky issue) if no `gDirectory->cd("root:/");`
  TH1F* h = new TH1F("h", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h->Sumw2();
  RooDataSet* dsh;
  std::vector<TH1F*> hBenr(lxydis::lxycut[lxyvar].size()); 
  for(int i=0;i<lxydis::lxycut[lxyvar].size();i++) 
    { hBenr[i] = new TH1F(Form("hBenr-%d", i), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); hBenr[i]->Sumw2(); }
  std::vector<RooDataSet*> dshBenr(lxydis::lxycut[lxyvar].size());
  TH1F* hmcp_a = new TH1F("hmcp_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); hmcp_a->Sumw2();
  RooDataSet* dshmcp_a;
  TH1F* hmcp_b = new TH1F("hmcp_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); hmcp_b->Sumw2();
  RooDataSet* dshmcp_b;
  TH1F* hlxymcnp_a = new TH1F("hlxymcnp_a", Form(";%s;Probability", lxydis::vars[lxyvar].c_str()), lxyxbins["lxynonprompt"].size()-1, lxyxbins["lxynonprompt"].data()); hlxymcnp_a->Sumw2();
  TH1F* hlxymcnp_b = new TH1F("hlxymcnp_b", Form(";%s;Probability", lxydis::vars[lxyvar].c_str()), lxyxbins["lxynonprompt"].size()-1, lxyxbins["lxynonprompt"].data()); hlxymcnp_b->Sumw2();
  TH1F* hlxymcp_a = new TH1F("hlxymcp_a", Form(";%s;Probability", lxydis::vars[lxyvar].c_str()), lxyxbins["lxyprompt"].size()-1, lxyxbins["lxyprompt"].data()); hlxymcp_a->Sumw2();
  TH1F* hlxymcp_b = new TH1F("hlxymcp_b", Form(";%s;Probability", lxydis::vars[lxyvar].c_str()), lxyxbins["lxyprompt"].size()-1, lxyxbins["lxyprompt"].data()); hlxymcp_b->Sumw2();

  //
  TTree* ntmix = fitX::getnt(input, "Bfinder/ntmix"); if(!ntmix) { return; }
  TTree* ntmixmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntmix"); if(!ntmixmcp_a) { return; }
  TTree* ntmixmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntmix"); if(!ntmixmcp_b) { return; }
  TTree* ntmixmcnp_a = fitX::getnt(inputmcnp_a, "Bfinder/ntmix"); if(!ntmixmcnp_a) { return; }
  TTree* ntmixmcnp_b = fitX::getnt(inputmcnp_b, "Bfinder/ntmix"); if(!ntmixmcnp_b) { return; }
  TTree* ntmix_flatten = fitX::getnt(input_flatten, "ntmix_flatten", false); if(!ntmix_flatten) { return; }
  TTree* ntmixmcp_a_flatten = fitX::getnt(inputmcp_a_flatten, "ntmix_flatten", false); if(!ntmixmcp_a_flatten) { return; }
  TTree* ntmixmcp_b_flatten = fitX::getnt(inputmcp_b_flatten, "ntmix_flatten", false); if(!ntmixmcp_b_flatten) { return; }
  
  gDirectory->cd("root:/");
  RooWorkspace* ww = new RooWorkspace("ww");

  std::string cutreco = Form("(%s) && Bmass >= %f && Bmass < %f && Bpt>%f && Bpt<%f && TMath::Abs(By)>=%f && TMath::Abs(By)<%f && hiBin>=%f && hiBin<=%f", cut.c_str(), 
                             fitX::BIN_MIN, fitX::BIN_MAX, 
                             fitX::ptmincut, fitX::ptmaxcut, 
                             fitX::ymincut, fitX::ymaxcut, 
                             fitX::centmincut*2, fitX::centmaxcut*2);
  std::string cutmcreco = Form("%s && Bgen>=23333 && BgencollisionId==0", cutreco.c_str());

  //
  std::cout<<" == data ==>"<<std::endl;
  fitX::printhist(ntmix_flatten);
  ntmix->Project("h", "Bmass", TCut(cutreco.c_str()));
  fitX::printhist(h);
  TTree* ntmix_skimh = (TTree*)ntmix_flatten->CopyTree(TCut(cutreco.c_str())); ntmix_skimh->SetName("ntmix_skimh");
  fitX::printhist(ntmix_skimh);
  dsh = new RooDataSet("dsh", "", ntmix_skimh, RooArgSet(*mass));
  ww->import(*dsh);
  for(int i=0; i<lxydis::lxycut[lxyvar].size(); i++)
    {
      ntmix->Project(Form("hBenr-%d", i), "Bmass", TCut(Form("%s && %s > %f", cutreco.c_str(), lxydis::formulas[lxyvar].c_str(), lxydis::lxycut[lxyvar][i])));
      fitX::printhist(hBenr[i]);
      TTree* ntmix_skimhBenr = (TTree*)ntmix_flatten->CopyTree(TCut(Form("%s && %s > %f", cutreco.c_str(), lxydis::formulas[lxyvar].c_str(), lxydis::lxycut[lxyvar][i]))); ntmix_skimhBenr->SetName(Form("ntmix_skimhBenr-%d", i));
      fitX::printhist(ntmix_skimhBenr);
      dshBenr[i] = new RooDataSet(Form("dshBenr-%d", i), "", ntmix_skimhBenr, RooArgSet(*mass));
      ww->import(*(dshBenr[i]));
    }
  //
  std::cout<<" == mcp_a ==>"<<std::endl;
  ntmixmcp_a->Project("hmcp_a", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // !! weight
  fitX::printhist(hmcp_a);
  TTree* ntmixmcp_a_skim = (TTree*)ntmixmcp_a_flatten->CopyTree(TCut(cutmcreco.c_str())); ntmixmcp_a_skim->SetName("ntmixmcp_a_skim");
  fitX::printhist(ntmixmcp_a_skim);
  dshmcp_a = new RooDataSet("dshmcp_a", "", ntmixmcp_a_skim, RooArgSet(*massmc_a, *pthatweight), 0, "pthatweight");
  ww->import(*dshmcp_a);
  ntmixmcp_a->Project("hlxymcp_a", lxydis::formulas[lxyvar].c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  fitX::printhist(hlxymcp_a);
  //
  std::cout<<" == mcp_b ==>"<<std::endl;
  ntmixmcp_b->Project("hmcp_b", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // !! weight
  fitX::printhist(hmcp_b);
  TTree* ntmixmcp_b_skim = (TTree*)ntmixmcp_b_flatten->CopyTree(TCut(cutmcreco.c_str())); ntmixmcp_b_skim->SetName("ntmixmcp_b_skim");
  fitX::printhist(ntmixmcp_b_skim);
  dshmcp_b = new RooDataSet("dshmcp_b", "", ntmixmcp_b_skim, RooArgSet(*massmc_b, *pthatweight), 0, "pthatweight");
  ww->import(*dshmcp_b);
  ntmixmcp_b->Project("hlxymcp_b", lxydis::formulas[lxyvar].c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  fitX::printhist(hlxymcp_b);

  //
  std::cout<<" == mcnp_a ==>"<<std::endl;
  ntmixmcnp_a->Project("hlxymcnp_a", lxydis::formulas[lxyvar].c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  fitX::printhist(hlxymcnp_a);

  std::cout<<" == mcnp_b ==>"<<std::endl;
  ntmixmcnp_b->Project("hlxymcnp_b", lxydis::formulas[lxyvar].c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  fitX::printhist(hlxymcnp_b);

  std::string outputname = "rootfiles/" + output + fitX::tagname() + "/fprompt_savehist.root";
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  h->Write();
  for(auto& hh : hBenr) hh->Write();
  hmcp_a->Write();
  hmcp_b->Write();
  hlxymcnp_a->Write();
  hlxymcnp_b->Write();
  hlxymcp_a->Write();
  hlxymcp_b->Write();
  outf->cd();
  gDirectory->Add(ww);
  ww->Write();
  ww->Print();
  outf->cd();
  TTree* info = new TTree("info", "cut info");
  info->Branch("input", &input);
  info->Branch("inputmcp_a", &inputmcp_a);
  info->Branch("inputmcp_b", &inputmcp_b);
  info->Branch("inputmcnp_a", &inputmcnp_a);
  info->Branch("inputmcnp_b", &inputmcnp_b);
  info->Branch("cutreco", &cutreco);
  info->Branch("cutmcreco", &cutmcreco);
  info->Fill();
  info->Write();
  fitX::write();
  outf->Close();
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==15) { 
    fitX::init(atof(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atof(argv[14]));
    fprompt_savehist(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]); return 0; }
  return 1;
}

