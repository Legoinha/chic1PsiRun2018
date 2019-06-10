#include "xjjrootuti.h"
#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include <string>

namespace MCeff
{
  std::vector<float> ptBins = {15, 20, 25, 30, 40, 60, 80, 100};
  const int nPtBins = ptBins.size() - 1;
}

void MCefficiency(std::string input_a, std::string input_b, std::string cut, std::string cutgen,
                  std::string output)
{
  TFile* inf_a = TFile::Open(input_a.c_str());
  TTree* ntmix_a = (TTree*)inf_a->Get("Bfinder/ntmix");
  ntmix_a->AddFriend("hiEvtAnalyzer/HiTree");
  ntmix_a->AddFriend("hltanalysis/HltTree");
  ntmix_a->AddFriend("skimanalysis/HltTree");
  ntmix_a->AddFriend("dataset/mva");
  TTree* ntGen_a = (TTree*)inf_a->Get("Bfinder/ntGen");
  ntGen_a->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen_a->AddFriend("dataset/mva");

  TFile* inf_b = TFile::Open(input_b.c_str());
  TTree* ntmix_b = (TTree*)inf_b->Get("Bfinder/ntmix");
  ntmix_b->AddFriend("hiEvtAnalyzer/HiTree");
  ntmix_b->AddFriend("hltanalysis/HltTree");
  ntmix_b->AddFriend("skimanalysis/HltTree");
  ntmix_b->AddFriend("dataset/mva");
  TTree* ntGen_b = (TTree*)inf_b->Get("Bfinder/ntGen");
  ntGen_b->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen_b->AddFriend("dataset/mva");

  TH2F* hempty = new TH2F("hempty", ";p_{T} (GeV/c);#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 10, MCeff::ptBins[0], MCeff::ptBins[MCeff::nPtBins], 10, 0, 0.2);
  xjjroot::sethempty(hempty, 0, 0.3);
  TH2F* hempty_incl = new TH2F("hempty_incl", ";;#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 5, 0, 5, 10, 0, 0.2);
  xjjroot::sethempty(hempty_incl, 0, 0.3);
  hempty_incl->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hempty_incl->GetXaxis()->SetBinLabel(4, "X(3872)");
  hempty_incl->GetXaxis()->SetLabelSize(hempty_incl->GetXaxis()->GetLabelSize()*1.5);

  TH1F* hmc_a = new TH1F("hmc_a", ";p_{T} (GeV/c);", MCeff::nPtBins, MCeff::ptBins.data()); hmc_a->Sumw2();
  ntmix_a->Project("hmc_a", "Bpt", TCut("pthatweight*Ncoll")*(TCut(cut.c_str())&&TCut("Bgen==23333&&BgencollisionId==0&&hiBin>=0&&hiBin<180")));
  TH1F* hmc_incl_a = new TH1F("hmc_incl_a", "", 5, 0, 5); hmc_incl_a->Sumw2();
  ntmix_a->Project("hmc_incl_a", "1", TCut("pthatweight*Ncoll")*(TCut(cut.c_str())&&TCut("Bgen==23333&&BgencollisionId==0&&hiBin>=0&&hiBin<180")));
  TH1F* hgen_a = new TH1F("hgen_a", ";p_{T} (GeV/c);", MCeff::nPtBins, MCeff::ptBins.data());
  ntGen_a->Project("hgen_a", "Gpt", TCut("pthatweight*Ncoll")*(TCut(cutgen.c_str())&&TCut("hiBin<=200")));
  TH1F* hgen_incl_a = new TH1F("hgen_incl_a", "", 5, 0, 5); hgen_incl_a->Sumw2();
  ntGen_a->Project("hgen_incl_a", "1", TCut("pthatweight*Ncoll")*(TCut(cutgen.c_str())&&TCut("hiBin<=200")));
  TEfficiency* greff_a = new TEfficiency(*hmc_a, *hgen_a); greff_a->SetName("greff_a");
  xjjroot::setthgrstyle(greff_a, 0, 0, 0, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
  TEfficiency* greff_incl_a = new TEfficiency(*hmc_incl_a, *hgen_incl_a); greff_incl_a->SetName("greff_incl_a");
  xjjroot::setthgrstyle(greff_incl_a, fitX::color_a, 20, 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);

  TH1F* hmc_b = new TH1F("hmc_b", ";p_{T} (GeV/c);", MCeff::nPtBins, MCeff::ptBins.data()); hmc_b->Sumw2();
  ntmix_b->Project("hmc_b", "Bpt", TCut("pthatweight*Ncoll")*(TCut(cut.c_str())&&TCut("Bgen==23333&&BgencollisionId==0&&hiBin>=0&&hiBin<180")));
  TH1F* hmc_incl_b = new TH1F("hmc_incl_b", "", 5, 0, 5); hmc_incl_b->Sumw2();
  ntmix_b->Project("hmc_incl_b", "3", TCut("pthatweight*Ncoll")*(TCut(cut.c_str())&&TCut("Bgen==23333&&BgencollisionId==0&&hiBin>=0&&hiBin<180")));
  TH1F* hgen_b = new TH1F("hgen_b", ";p_{T} (GeV/c);", MCeff::nPtBins, MCeff::ptBins.data());
  ntGen_b->Project("hgen_b", "Gpt", TCut("pthatweight*Ncoll")*(TCut(cutgen.c_str())&&TCut("hiBin>=0 && hiBin<180")));
  TH1F* hgen_incl_b = new TH1F("hgen_incl_b", "", 5, 0, 5); hgen_incl_b->Sumw2();
  ntGen_b->Project("hgen_incl_b", "3", TCut("pthatweight*Ncoll")*(TCut(cutgen.c_str())&&TCut("hiBin>=0 && hiBin<180")));
  TEfficiency* greff_b = new TEfficiency(*hmc_b, *hgen_b); greff_b->SetName("greff_b");
  xjjroot::setthgrstyle(greff_b, 0, 0, 0, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
  TEfficiency* greff_incl_b = new TEfficiency(*hmc_incl_b, *hgen_incl_b); greff_incl_b->SetName("greff_incl_b");
  xjjroot::setthgrstyle(greff_incl_b, fitX::color_b, 20, 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
  TH1F* heff_a = (TH1F*)hmc_a->Clone("heff_a");
  heff_a->Divide(hgen_a);
  xjjroot::setthgrstyle(heff_a, 0, 0, 0, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
  TH1F* heff_incl_a = (TH1F*)hmc_incl_a->Clone("heff_incl_a");
  heff_incl_a->Divide(hgen_incl_a);
  xjjroot::setthgrstyle(heff_incl_a, fitX::color_a, 20, 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
  TH1F* heff_b = (TH1F*)hmc_b->Clone("heff_b");
  heff_b->Divide(hgen_b);
  xjjroot::setthgrstyle(heff_b, 0, 0, 0, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
  TH1F* heff_incl_b = (TH1F*)hmc_incl_b->Clone("heff_incl_b");
  heff_incl_b->Divide(hgen_incl_b);
  xjjroot::setthgrstyle(heff_incl_b, fitX::color_b, 20, 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);

  TLegend* leg = new TLegend(0.70, 0.20, 1.10, 0.32);
  xjjroot::setleg(leg, 0.042);
  leg->AddEntry(greff_a, "#psi(2S)", "fl");
  leg->AddEntry(greff_b, "X(3872)", "fl");

  xjjroot::setgstyle();
  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hempty->Draw();
  greff_a->Draw("same3");
  greff_a->Draw("samelX");
  greff_b->Draw("same3");
  greff_b->Draw("samelX");
  leg->Draw();
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  xjjroot::drawtex(0.92, 0.84, "|y| < 1.5", 0.042, 32, 62);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();

  c->cd(2);
  hempty_incl->Draw();
  greff_incl_a->Draw("same ple");
  greff_incl_b->Draw("same ple");
  float effval_a = greff_incl_a->GetEfficiency(2);
  xjjroot::drawtex(0.42, effval_a/0.2 + 0.2, Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-3}", greff_incl_a->GetEfficiency(2)*1000, greff_incl_a->GetEfficiencyErrorUp(2)*1000, greff_incl_a->GetEfficiencyErrorLow(2)*1000), 0.042, 22, 62, fitX::color_a);
  float effval_b = greff_incl_b->GetEfficiency(4);
  xjjroot::drawtex(0.72, effval_b/0.2 + 0.2, Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-3}", greff_incl_b->GetEfficiency(4)*1000, greff_incl_b->GetEfficiencyErrorUp(4)*1000, greff_incl_b->GetEfficiencyErrorLow(4)*1000), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  xjjroot::drawtex(0.92, 0.84, "|y| < 1.5", 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.77, "p_{T} > 15 GeV/c", 0.042, 32, 62);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();

  c->SaveAs(Form("plots/ceff_%s.pdf", output.c_str()));

  TFile* outf = new TFile(Form("rootfiles/root_eff_%s.root", output.c_str()), "recreate");
  outf->cd();
  greff_a->Write();
  greff_b->Write();
  greff_incl_a->Write();
  greff_incl_b->Write();
  heff_a->Write();
  heff_b->Write();
  heff_incl_a->Write();
  heff_incl_b->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==6) { MCefficiency(argv[1], argv[2], argv[3], argv[4], argv[5]); return 0; }
  return 1;
}
