#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH2F.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "fit.h"
#include "project.h"
#include "MCefficiency.h"
#include "xjjcuti.h"
#include "xjjrootuti.h"

void closure_savehist(std::string inputmcp_a, std::string inputmcp_b,
                      std::string cut, std::string cutgen, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::string inputmcp_a_flatten = xjjc::str_replaceall(inputmcp_a, ".root", "_flatten.root");
  std::string inputmcp_b_flatten = xjjc::str_replaceall(inputmcp_b, ".root", "_flatten.root");

  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooRealVar* massmc_a = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_L, fitX::BIN_MAX_L);
  RooRealVar* massmc_b = new RooRealVar("Bmass", "massmc_b", fitX::BIN_MIN_H, fitX::BIN_MAX_H);
  RooRealVar* pthatweight = new RooRealVar("pthatweight", "pthatweight", 0, 10.);

  MCeff::MCefficiency* mceff_a = new MCeff::MCefficiency("_a");
  MCeff::MCefficiency* mceff_b = new MCeff::MCefficiency("_b");
  std::vector<float> bins = mceff_a->ptbins();
  int nbins = mceff_a->nbins();
  std::vector<TH1F*> h_a(nbins, 0), h_b(nbins, 0), hmcp_a(nbins, 0), hmcp_b(nbins, 0);
  std::vector<RooDataSet*> dsh_a(nbins, 0), dsh_b(nbins, 0), dshmcp_a(nbins, 0), dshmcp_b(nbins, 0);
  for(int i=0; i<nbins; i++)
    {
      h_a[i] = new TH1F(Form("h_a_%d", i), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h_a[i]->Sumw2();
      h_b[i] = new TH1F(Form("h_b_%d", i), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h_b[i]->Sumw2();
      hmcp_a[i] = new TH1F(Form("hmcp_a_%d", i), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); hmcp_a[i]->Sumw2();
      hmcp_b[i] = new TH1F(Form("hmcp_b_%d", i), Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); hmcp_b[i]->Sumw2();
    }

  TTree* ntmixmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntmix"); if(!ntmixmcp_a) { return; }
  TTree* ntmixmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntmix"); if(!ntmixmcp_b) { return; }
  TTree* ntGenmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntGen"); if(!ntGenmcp_a) { return; }
  TTree* ntGenmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntGen"); if(!ntGenmcp_b) { return; }
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
  std::string cutmcgen = Form("(%s) && Gpt>%f && Gpt<%f && TMath::Abs(Gy)>=%f && TMath::Abs(Gy)<%f && hiBin>=%f && hiBin<=%f && GisSignal==7 && GcollisionId==0", cutgen.c_str(), fitX::ptmincut, fitX::ptmaxcut, fitX::ymincut, fitX::ymaxcut, fitX::centmincut*2, fitX::centmaxcut*2);

  std::cout<<" == \"data\" (_a) ==>"<<std::endl;
  xjjroot::printhist(ntmixmcp_a_flatten);
  for(int i=0; i<nbins; i++)
    {
      std::string cutrecoi = cutreco + std::string(Form("&& Bpt>%f && Bpt<%f", bins[i], bins[i+1]));
      std::string cutmcrecoi = cutmcreco + std::string(Form("&& Bpt>%f && Bpt<%f", bins[i], bins[i+1]));
      // "data"
      ntmixmcp_a->Project(Form("h_a_%d", i), "Bmass", TCut("pthatweight")*TCut(cutrecoi.c_str()));
      xjjroot::printhist(h_a[i]);
      TTree* ntmix_a_skim = (TTree*)ntmixmcp_a_flatten->CopyTree(TCut(cutrecoi.c_str())); ntmix_a_skim->SetName(Form("ntmix_a_skim_%d", i));
      xjjroot::printhist(ntmix_a_skim);
      dsh_a[i] = new RooDataSet(Form("dsh_a_%d", i), "", ntmix_a_skim, RooArgSet(*mass, *pthatweight), 0, "pthatweight");
      ww->import(*dsh_a[i]);
      // mc
      ntmixmcp_a->Project(Form("hmcp_a_%d", i), "Bmass", TCut("pthatweight")*TCut(cutmcrecoi.c_str()));
      xjjroot::printhist(hmcp_a[i]);
      TTree* ntmixmcp_a_skim = (TTree*)ntmixmcp_a_flatten->CopyTree(TCut(cutmcrecoi.c_str())); ntmixmcp_a_skim->SetName(Form("ntmixmcp_a_skim_%d", i));
      xjjroot::printhist(ntmixmcp_a_skim);
      dshmcp_a[i] = new RooDataSet(Form("dshmcp_a_%d", i), "", ntmixmcp_a_skim, RooArgSet(*massmc_a, *pthatweight), 0, "pthatweight");
      ww->import(*dshmcp_a[i]);
      delete ntmix_a_skim;
      delete ntmixmcp_a_skim;
    }

  std::cout<<" == \"data\" (_b) ==>"<<std::endl;
  xjjroot::printhist(ntmixmcp_b_flatten);
  for(int i=0; i<nbins; i++)
    {
      std::string cutrecoi = cutreco + std::string(Form("&& Bpt>%f && Bpt<%f", bins[i], bins[i+1]));
      std::string cutmcrecoi = cutmcreco + std::string(Form("&& Bpt>%f && Bpt<%f", bins[i], bins[i+1]));
      // "data"
      ntmixmcp_b->Project(Form("h_b_%d", i), "Bmass", TCut("pthatweight")*TCut(cutrecoi.c_str()));
      xjjroot::printhist(h_b[i]);
      TTree* ntmix_b_skim = (TTree*)ntmixmcp_b_flatten->CopyTree(TCut(cutrecoi.c_str())); ntmix_b_skim->SetName(Form("ntmix_b_skim_%d", i));
      xjjroot::printhist(ntmix_b_skim);
      dsh_b[i] = new RooDataSet(Form("dsh_b_%d", i), "", ntmix_b_skim, RooArgSet(*mass, *pthatweight), 0, "pthatweight");
      ww->import(*dsh_b[i]);
      // mc
      ntmixmcp_b->Project(Form("hmcp_b_%d", i), "Bmass", TCut("pthatweight")*TCut(cutmcrecoi.c_str()));
      xjjroot::printhist(hmcp_b[i]);
      TTree* ntmixmcp_b_skim = (TTree*)ntmixmcp_b_flatten->CopyTree(TCut(cutmcrecoi.c_str())); ntmixmcp_b_skim->SetName(Form("ntmixmcp_b_skim_%d", i));
      xjjroot::printhist(ntmixmcp_b_skim);
      dshmcp_b[i] = new RooDataSet(Form("dshmcp_b_%d", i), "", ntmixmcp_b_skim, RooArgSet(*massmc_b, *pthatweight), 0, "pthatweight");
      ww->import(*dshmcp_b[i]);
      delete ntmix_b_skim;
      delete ntmixmcp_b_skim;
    }

  std::cout<<" == mcp_a ==>"<<std::endl;
  ntmixmcp_a->Project(mceff_a->heffmc()->GetName(), "Bpt", TCut("pthatweight")*TCut(cutmcreco.c_str()));
  xjjroot::printhist(mceff_a->heffmc());
  // ntmixmcp_a->Project(mceff_a->heffmc_incl()->GetName(), Form("%d", fitX::ibin_a-1), TCut("pthatweight")*TCut(cutmcreco.c_str()));
  // xjjroot::printhist(mceff_a->heffmc_incl());

  std::cout<<" == mcp_b ==>"<<std::endl;
  ntmixmcp_b->Project(mceff_b->heffmc()->GetName(), "Bpt", TCut("pthatweight")*TCut(cutmcreco.c_str()));
  xjjroot::printhist(mceff_b->heffmc());
  // ntmixmcp_b->Project(mceff_b->heffmc_incl()->GetName(), Form("%d", fitX::ibin_b-1), TCut("pthatweight")*TCut(cutmcreco.c_str()));
  // xjjroot::printhist(mceff_b->heffmc_incl());

  std::cout<<" == mcgenp_a ==>"<<std::endl;
  ntGenmcp_a->Project(mceff_a->heffgen()->GetName(), "Gpt", TCut("pthatweight")*TCut(cutmcgen.c_str()));
  xjjroot::printhist(mceff_a->heffgen());
  // ntGenmcp_a->Project(mceff_a->heffgen_incl()->GetName(), Form("%d", fitX::ibin_a-1), TCut("pthatweight")*TCut(cutmcgen.c_str()));
  // xjjroot::printhist(mceff_a->heffgen_incl());

  std::cout<<" == mcgenp_b ==>"<<std::endl;
  ntGenmcp_b->Project(mceff_b->heffgen()->GetName(), "Gpt", TCut("pthatweight")*TCut(cutmcgen.c_str()));
  xjjroot::printhist(mceff_b->heffgen());
  // ntGenmcp_b->Project(mceff_b->heffgen_incl()->GetName(), Form("%d", fitX::ibin_b-1), TCut("pthatweight")*TCut(cutmcgen.c_str()));
  // xjjroot::printhist(mceff_b->heffgen_incl());

  std::string outputname = "rootfiles/" + output + fitX::tagname() + "/closure_savehist.root";
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  for(auto& hh : h_a) hh->Write();
  for(auto& hh : h_b) hh->Write();
  for(auto& hh : hmcp_a) hh->Write();
  for(auto& hh : hmcp_b) hh->Write();
  mceff_a->heffmc()->Write();
  mceff_a->heffmc_incl()->Write();
  mceff_a->heffgen()->Write();
  mceff_a->heffgen_incl()->Write();
  mceff_b->heffmc()->Write();
  mceff_b->heffmc_incl()->Write();
  mceff_b->heffgen()->Write();
  mceff_b->heffgen_incl()->Write();
  outf->cd();
  gDirectory->Add(ww);
  ww->Write();
  ww->Print();
  outf->cd();
  TTree* info = new TTree("info", "cut info");
  info->Branch("inputmcp_a", &inputmcp_a);
  info->Branch("inputmcp_b", &inputmcp_b);
  info->Branch("cutreco", &cutreco);
  info->Branch("cutmcreco", &cutmcreco);
  info->Branch("cutmcgen", &cutmcgen);
  info->Fill();
  info->Write();
  fitX::write();
  outf->Close();
  std::cout<<std::endl;

}

int main(int argc, char* argv[])
{
  if(argc==12) {
    fitX::init(atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]), atof(argv[11]));
    closure_savehist(argv[1], argv[2], argv[3], argv[4], argv[5]); return 0; }
  return 1;
}
