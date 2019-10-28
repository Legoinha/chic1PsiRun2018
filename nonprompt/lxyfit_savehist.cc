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

int nsmear = 10;
void lxyfit_savehist(std::string input, std::string inputmcp_a, std::string inputmcp_b, std::string inputmcnp_a, std::string inputmcnp_b,
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
  std::string mcweight_a = "(pthatweight*Ncoll*(-1.502425+0.049572*Bgenpt+29.751642/Bgenpt))";
  std::string mcweight_b = "(pthatweight*Ncoll*(-10.797509+0.239998*Bgenpt+137.456058/Bgenpt))";
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

  std::vector<TH1F*> hmcpdis_a(nsmear); 
  for(int i=0;i<nsmear;i++) { hmcpdis_a[i] = new TH1F(Form("hmcpdis_a_%d", i), Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data()); }
  std::vector<TH1F*> hmcpdis_b(nsmear); 
  for(int i=0;i<nsmear;i++) { hmcpdis_b[i] = new TH1F(Form("hmcpdis_b_%d", i), Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data()); }
  std::vector<TH1F*> hmcnpdis_a(nsmear); 
  for(int i=0;i<nsmear;i++) { hmcnpdis_a[i] = new TH1F(Form("hmcnpdis_a_%d", i), Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data()); }
  std::vector<TH1F*> hmcnpdis_b(nsmear); 
  for(int i=0;i<nsmear;i++) { hmcnpdis_b[i] = new TH1F(Form("hmcnpdis_b_%d", i), Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data()); }

  TH1F* hbkgdis_a = new TH1F("hbkgdis_a", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());
  TH1F* hbkgdis_b = new TH1F("hbkgdis_b", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());

  //
  TTree* ntmix = fitX::getnt(input, "Bfinder/ntmix"); if(!ntmix) { return; }
  TTree* ntmixmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntmix"); if(!ntmixmcp_a) { return; }
  TTree* ntmixmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntmix"); if(!ntmixmcp_b) { return; }
  TTree* ntmixmcnp_a = fitX::getnt(inputmcnp_a, "Bfinder/ntmix"); if(!ntmixmcnp_a) { return; }
  TTree* ntmixmcnp_b = fitX::getnt(inputmcnp_b, "Bfinder/ntmix"); if(!ntmixmcnp_b) { return; }
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
      cutrecoi[i] = cutreco + std::string(Form("&& %s>%f && %s<%f", 
                                               vv->formula().c_str(), vv->vars()[i], 
                                               vv->formula().c_str(), vv->vars()[i+1]));
    }
  std::string cutbkgreco_a = cutreco + std::string(Form("&& fabs(Bmass-%f)>0.02 && fabs(Bmass-%f)<0.05", fitX::FIT_MASS_PSI2S, fitX::FIT_MASS_PSI2S));
  std::string cutbkgreco_b = cutreco + std::string(Form("&& fabs(Bmass-%f)>0.02 && fabs(Bmass-%f)<0.05", fitX::FIT_MASS_X, fitX::FIT_MASS_X));

  std::cout<<" == data ==>"<<std::endl;
  xjjroot::printhist(ntmix_flatten, 13);
  for(int i=0; i<vv->n()-1; i++)
    {
      // std::cout<<cutrecoi[i]<<std::endl;
      ntmix->Project(Form("h%d", i), "Bmass", TCut(cutrecoi[i].c_str()));
      xjjroot::printhist(h[i], 13);
      TTree* ntmix_skimh = (TTree*)ntmix_flatten->CopyTree(TCut(cutrecoi[i].c_str())); ntmix_skimh->SetName(Form("ntmix_skimh%d", i));
      xjjroot::printhist(ntmix_skimh, 13);
      dsh[i] = new RooDataSet(Form("dsh%d", i), "", ntmix_skimh, RooArgSet(*mass, *varr));
      ww->import(*(dsh[i]));
    }
  ntmix->Project(hbkgdis_a->GetName(), vv->formula().c_str(), TCut(cutbkgreco_a.c_str()));
  xjjroot::printhist(hbkgdis_a, 13);
  ntmix->Project(hbkgdis_b->GetName(), vv->formula().c_str(), TCut(cutbkgreco_b.c_str()));
  xjjroot::printhist(hbkgdis_b, 13);

  std::cout<<" == mcp_a ==>"<<std::endl;
  ntmixmcp_a->Project("hmcp_a", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // !! weight
  xjjroot::printhist(hmcp_a, 13);
  TTree* ntmixmcp_a_skim = (TTree*)ntmixmcp_a_flatten->CopyTree(TCut(cutmcreco.c_str())); ntmixmcp_a_skim->SetName("ntmixmcp_a_skim");
  xjjroot::printhist(ntmixmcp_a_skim, 13);
  dshmcp_a = new RooDataSet("dshmcp_a", "", ntmixmcp_a_skim, RooArgSet(*massmc_a, *pthatweight), "pthatweight");
  ww->import(*dshmcp_a);
  // ntmixmcp_a->Project(hmcpdis_a->GetName(), vv->formula().c_str(), TCut(mcweight_a.c_str())*TCut(cutmcreco.c_str()));
  // xjjroot::printhist(hmcpdis_a, 13);
  for(int i=0; i<nsmear; i++) 
    {
      ntmixmcp_a->Project(hmcpdis_a[i]->GetName(), Form("%s*%f", vv->formula().c_str(), 1.0+0.04*i), TCut(mcweight_a.c_str())*TCut(cutmcreco.c_str()));
      xjjroot::printhist(hmcpdis_a[i], 13);
    }

  std::cout<<" == mcp_b ==>"<<std::endl;
  ntmixmcp_b->Project("hmcp_b", "Bmass", TCut("pthatweight")*TCut(cutmcreco.c_str())); // !! weight
  xjjroot::printhist(hmcp_b, 13);
  TTree* ntmixmcp_b_skim = (TTree*)ntmixmcp_b_flatten->CopyTree(TCut(cutmcreco.c_str())); ntmixmcp_b_skim->SetName("ntmixmcp_b_skim");
  xjjroot::printhist(ntmixmcp_b_skim, 13);
  dshmcp_b = new RooDataSet("dshmcp_b", "", ntmixmcp_b_skim, RooArgSet(*massmc_b, *pthatweight), "pthatweight");
  ww->import(*dshmcp_b);
  // ntmixmcp_b->Project(hmcpdis_b->GetName(), vv->formula().c_str(), TCut(mcweight_b.c_str())*TCut(cutmcreco.c_str()));
  // xjjroot::printhist(hmcpdis_b, 13);
  for(int i=0; i<nsmear; i++) 
    {
      ntmixmcp_b->Project(hmcpdis_b[i]->GetName(), Form("%s*%f", vv->formula().c_str(), 1.0+0.04*i), TCut(mcweight_b.c_str())*TCut(cutmcreco.c_str()));
      xjjroot::printhist(hmcpdis_b[i], 13);
    }

  std::cout<<" == mcnp_a ==>"<<std::endl;
  // ntmixmcnp_a->Project(hmcnpdis_a->GetName(), vv->formula().c_str(), TCut(mcweight_a.c_str())*TCut(cutmcreco.c_str()));
  // xjjroot::printhist(hmcnpdis_a, 13);
  for(int i=0; i<nsmear; i++) 
    {
      ntmixmcnp_a->Project(hmcnpdis_a[i]->GetName(), Form("%s*%f", vv->formula().c_str(), 1.0+0.04*i), TCut(mcweight_a.c_str())*TCut(cutmcreco.c_str()));
      xjjroot::printhist(hmcnpdis_a[i], 13);
    }

  std::cout<<" == mcnp_b ==>"<<std::endl;
  // ntmixmcnp_b->Project(hmcnpdis_b->GetName(), vv->formula().c_str(), TCut(mcweight_b.c_str())*TCut(cutmcreco.c_str()));
  // xjjroot::printhist(hmcnpdis_b, 13);
  for(int i=0; i<nsmear; i++) 
    {
      ntmixmcnp_b->Project(hmcnpdis_b[i]->GetName(), Form("%s*%f", vv->formula().c_str(), 1.0+0.04*i), TCut(mcweight_b.c_str())*TCut(cutmcreco.c_str()));
      xjjroot::printhist(hmcnpdis_b[i], 13);
    }

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
  for(auto& hh : hmcpdis_a) hh->Write();
  for(auto& hh : hmcpdis_b) hh->Write();
  for(auto& hh : hmcnpdis_a) hh->Write();
  for(auto& hh : hmcnpdis_b) hh->Write();
  hbkgdis_a->Write();
  hbkgdis_b->Write();
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
  if(argc==15)
    {
      fitX::init(atof(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atof(argv[14]));
      lxyfit_savehist(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
      return 0;
    }
  return 1;
}

