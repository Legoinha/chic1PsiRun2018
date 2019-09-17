#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TCut.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <string>
#include <vector>

#include "var.h"
#include "fit.h"
#include "project.h"

#include "xjjcuti.h"
#include "xjjrootuti.h"

void datamcmain(std::string input, std::string inputmcp_a, std::string inputmcp_b,
                std::string cut, std::string type, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::cout<<cut<<std::endl;
  datamc::var* vv = new datamc::var(type.c_str());
  if(!vv->valid()) return;

  output += (fitX::tagname()+"/"+type);

  std::string input_flatten = xjjc::str_replaceall(input, ".root", "_flatten.root");
  std::string inputmcp_a_flatten = xjjc::str_replaceall(inputmcp_a, ".root", "_flatten.root");
  std::string inputmcp_b_flatten = xjjc::str_replaceall(inputmcp_b, ".root", "_flatten.root");

  std::string mcweight = "(pthatweight*Ncoll)";
  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooRealVar* massmc_a = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_L, fitX::BIN_MAX_L);
  RooRealVar* massmc_b = new RooRealVar("Bmass", "massmc_b", fitX::BIN_MIN_H, fitX::BIN_MAX_H);
  RooRealVar* pthatweight = new RooRealVar("pthatweight", "pthatweight", 0, 10.); // pthatweight range!!
  RooRealVar* varr = new RooRealVar(vv->formula().c_str(), vv->title().c_str(), vv->vars().front(), vv->vars().back());

  std::vector<RooDataSet*> dsh(vv->n()-1);
  RooDataSet* dshmcp_a;
  RooDataSet* dshmcp_b;
  std::vector<TH1F*> h(vv->n()-1);
  for(int i=0; i<vv->n()-1; i++) { h[i] = new TH1F(Form("h%d", i), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h[i]->Sumw2(); }
  TH1F* hmcp_a = new TH1F("hmcp_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); hmcp_a->Sumw2();
  TH1F* hmcp_b = new TH1F("hmcp_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); hmcp_b->Sumw2();

  TH1F* hmcdis_a = new TH1F("hmcdis_a", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());
  TH1F* hmcdis_b = new TH1F("hmcdis_b", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());
  TH1F* hmcfinedis_a = new TH1F("hmcfinedis_a", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->nfine(), vv->minfine(), vv->maxfine());
  TH1F* hmcfinedis_b = new TH1F("hmcfinedis_b", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->nfine(), vv->minfine(), vv->maxfine());
  TH1F* hbkgdis_a = new TH1F("hbkgdis_a", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());
  TH1F* hbkgdis_b = new TH1F("hbkgdis_b", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());

  //
  TTree* ntmix = fitX::getnt(input, "Bfinder/ntmix"); if(!ntmix) { return; }
  TTree* ntmixmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntmix"); if(!ntmixmcp_a) { return; }
  TTree* ntmixmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntmix"); if(!ntmixmcp_b) { return; }
  TTree* ntmix_flatten = fitX::getnt(input_flatten, "ntmix_flatten", false); if(!ntmix_flatten) { return; }
  TTree* ntmixmcp_a_flatten = fitX::getnt(inputmcp_a_flatten, "ntmix_flatten", false); if(!ntmixmcp_a_flatten) { return; }
  TTree* ntmixmcp_b_flatten = fitX::getnt(inputmcp_b_flatten, "ntmix_flatten", false); if(!ntmixmcp_b_flatten) { return; }

  gROOT->cd();
  RooWorkspace* ww = new RooWorkspace("ww");

  std::string cutreco = Form("(%s) && Bmass >= %f && Bmass < %f && Bpt>%f && Bpt<%f && TMath::Abs(By)>=%f && TMath::Abs(By)<%f && hiBin>=%f && hiBin<=%f", cut.c_str(),
                             fitX::BIN_MIN, fitX::BIN_MAX,
                             fitX::ptmincut, fitX::ptmaxcut,
                             fitX::ymincut, fitX::ymaxcut,
                             fitX::centmincut*2, fitX::centmaxcut*2);
  std::string cutmcreco = Form("%s && Bgen>=23333 && BgencollisionId==0", cutreco.c_str());
  std::vector<std::string> cutrecoi(vv->n()-1);
  for(int i=0; i<vv->n()-1; i++)
    {
      int icut = vv->gt()?i:i+1;
      cutrecoi[i] = cutreco + std::string(Form("&& %s%s%f && %s%s%f", 
                                               vv->formula().c_str(), (vv->gt()?">":"<"), vv->vars()[icut], 
                                               vv->formula().c_str(), (vv->gt()?"<":">"), (vv->gt()?vv->vars().back():vv->vars().front())));
    }
  std::string cutbkgreco_a = cutreco + std::string(Form("&& fabs(Bmass-%f)>0.02 && fabs(Bmass-%f)<0.05", fitX::FIT_MASS_PSI2S, fitX::FIT_MASS_PSI2S));
  std::string cutbkgreco_b = cutreco + std::string(Form("&& fabs(Bmass-%f)>0.02 && fabs(Bmass-%f)<0.05", fitX::FIT_MASS_X, fitX::FIT_MASS_X));

  std::cout<<" == data ==>"<<std::endl;
  xjjroot::printhist(ntmix_flatten);
  for(int i=0; i<vv->n()-1; i++)
    {
      // std::cout<<cutrecoi[i]<<std::endl;
      ntmix->Project(Form("h%d", i), "Bmass", TCut(cutrecoi[i].c_str()));
      xjjroot::printhist(h[i]);
      TTree* ntmix_skimh = (TTree*)ntmix_flatten->CopyTree(TCut(cutrecoi[i].c_str())); ntmix_skimh->SetName(Form("ntmix_skimh%d", i));
      xjjroot::printhist(ntmix_skimh);
      dsh[i] = new RooDataSet(Form("dsh%d", i), "", ntmix_skimh, RooArgSet(*mass, *varr));
      ww->import(*(dsh[i]));
    }
  ntmix->Project(hbkgdis_a->GetName(), vv->formula().c_str(), TCut(cutbkgreco_a.c_str()));
  xjjroot::printhist(hbkgdis_a);
  ntmix->Project(hbkgdis_b->GetName(), vv->formula().c_str(), TCut(cutbkgreco_b.c_str()));
  xjjroot::printhist(hbkgdis_b);

  std::cout<<" == mcp_a ==>"<<std::endl;
  ntmixmcp_a->Project("hmcp_a", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // !! weight
  xjjroot::printhist(hmcp_a);
  TTree* ntmixmcp_a_skim = (TTree*)ntmixmcp_a_flatten->CopyTree(TCut(cutmcreco.c_str())); ntmixmcp_a_skim->SetName("ntmixmcp_a_skim");
  xjjroot::printhist(ntmixmcp_a_skim);
  dshmcp_a = new RooDataSet("dshmcp_a", "", ntmixmcp_a_skim, RooArgSet(*massmc_a, *pthatweight), "pthatweight");
  ww->import(*dshmcp_a);
  ntmixmcp_a->Project(hmcdis_a->GetName(), vv->formula().c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  xjjroot::printhist(hmcdis_a);
  ntmixmcp_a->Project(hmcfinedis_a->GetName(), vv->formula().c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  xjjroot::printhist(hmcfinedis_a);

  std::cout<<" == mcp_b ==>"<<std::endl;
  ntmixmcp_b->Project("hmcp_b", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // !! weight
  xjjroot::printhist(hmcp_b);
  TTree* ntmixmcp_b_skim = (TTree*)ntmixmcp_b_flatten->CopyTree(TCut(cutmcreco.c_str())); ntmixmcp_b_skim->SetName("ntmixmcp_b_skim");
  xjjroot::printhist(ntmixmcp_b_skim);
  dshmcp_b = new RooDataSet("dshmcp_b", "", ntmixmcp_b_skim, RooArgSet(*massmc_b, *pthatweight), "pthatweight");
  ww->import(*dshmcp_b);
  ntmixmcp_b->Project(hmcdis_b->GetName(), vv->formula().c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  xjjroot::printhist(hmcdis_b);
  ntmixmcp_b->Project(hmcfinedis_b->GetName(), vv->formula().c_str(), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  xjjroot::printhist(hmcfinedis_b);

  std::string outputname = "rootfiles/" + output + "/datamc_savehist.root";
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  for(auto& hh : h) { hh->Write(); }
  hmcp_a->Write();
  hmcp_b->Write();
  outf->cd();
  gDirectory->Add(ww);
  ww->Write();
  ww->Print();
  outf->cd();
  hmcdis_a->Write();
  hmcdis_b->Write();
  hmcfinedis_a->Write();
  hmcfinedis_b->Write();
  hbkgdis_a->Write();
  hbkgdis_b->Write();
  TTree* info = new TTree("info", "cut info");
  info->Branch("input", &input);
  info->Branch("inputmcp_a", &inputmcp_a);
  info->Branch("inputmcp_b", &inputmcp_b);
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
  if(argc==13)
    {
      fitX::init(atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]));
      datamcmain(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
      return 0;
    }
  return 1;
}

