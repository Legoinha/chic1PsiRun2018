#include <TFile.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <string>

#include "xjjcuti.h"
#include "MCefficiency.h"
#include "fitX.h"

void drawkinematic();
void drawcompeff(std::string inputname_a, std::string inputname_b, std::string outputname)
{
  TFile* inf_a = new TFile(Form("%s", inputname_a.c_str()));
  MCeff::MCefficiency mceff_a(inf_a, "");
  MCeff::MCefficiency mceffweight_a(inf_a, "_weight");
  TFile* inf_b = new TFile(Form("%s", inputname_b.c_str()));
  MCeff::MCefficiency mceff_b(inf_b, "", 1);
  MCeff::MCefficiency mceffweight_b(inf_b, "_weight", 1);

  mceff_a.calceff();
  mceffweight_a.calceff();
  mceff_b.calceff();
  mceffweight_b.calceff();

  float ymaxeff = 0.2, ymaxeff_incl = 0.07;
  TH2F* hemptyeff = new TH2F("hemptyeff", ";p_{T} (GeV/c);#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 10, MCeff::ptBins[0], MCeff::ptBins[MCeff::nPtBins], 10, 0, ymaxeff);
  xjjroot::sethempty(hemptyeff, 0, 0.3);
  // TH2F* hemptyeff_incl = new TH2F("hemptyeff_incl", ";;#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 5, 0, 5, 10, 0, ymaxeff);
  // hemptyeff_incl->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  // hemptyeff_incl->GetXaxis()->SetBinLabel(4, "X(3872)");
  TH2F* hemptyeff_incl = new TH2F("hemptyeff_incl", ";;#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 10, MCeff::ptBins_incl[0], MCeff::ptBins_incl[MCeff::nPtBins_incl], 10, 0, ymaxeff_incl);
  xjjroot::sethempty(hemptyeff_incl, 0, 0.3);
  // hemptyeff_incl->GetXaxis()->SetLabelSize(hemptyeff_incl->GetXaxis()->GetLabelSize()*1.5);

  mceff_a.setstyle(xjjroot::mycolor_middle["azure"], 20, 1);
  mceffweight_a.setstyle(xjjroot::mycolor_middle["red"], 47, 1);
  mceff_b.setstyle(xjjroot::mycolor_middle["azure"], 20, 2);
  mceffweight_b.setstyle(xjjroot::mycolor_middle["red"], 47, 2);

  TLegend* legeff = new TLegend(0.35, 0.20, 0.95, 0.32);
  xjjroot::setleg(legeff, 0.038);
  legeff->SetNColumns(2);
  legeff->AddEntry(mceff_a.greff, "#psi(2S)", "fl");
  legeff->AddEntry(mceffweight_a.greff, "#psi(2S) weight p_{T}", "fl");
  legeff->AddEntry(mceff_b.greff, "X(3872)", "fl");
  legeff->AddEntry(mceffweight_b.greff, "X(3872) weight p_{T}", "fl");
  xjjroot::setgstyle();
  TCanvas* ceff = new TCanvas("ceff", "", 1200, 600);
  ceff->Divide(2, 1);
  ceff->cd(1);
  hemptyeff->Draw();
  mceff_a.greff->Draw("same3");
  mceff_a.greff->Draw("samelX");
  mceff_b.greff->Draw("same3");
  mceff_b.greff->Draw("samelX");
  mceffweight_a.greff->Draw("same3");
  mceffweight_a.greff->Draw("samelX");
  mceffweight_b.greff->Draw("same3");
  mceffweight_b.greff->Draw("samelX");
  legeff->Draw();
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  xjjroot::drawCMS("Simulation");
  ceff->cd(2);
  hemptyeff_incl->Draw();
  mceff_a.greff_incl->Draw("same ple");
  mceff_b.greff_incl->Draw("same ple");
  mceffweight_a.greff_incl->Draw("same ple");
  mceffweight_b.greff_incl->Draw("same ple");
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  drawkinematic();
  legeff->Draw();
  xjjroot::drawCMS("Simulation");
  ceff->SaveAs(Form("%s", outputname.c_str()));

}

void drawkinematic()
{
  xjjroot::drawtex(0.92, 0.84, Form("|y| < %s", xjjc::number_remove_zero(fitX::ycut).c_str()), 0.042, 32, 42);
  xjjroot::drawtex(0.92, 0.77, Form("%s < p_{T} < %s GeV/c", xjjc::number_remove_zero(fitX::ptmincut).c_str(), xjjc::number_remove_zero(fitX::ptmaxcut).c_str()), 0.042, 32, 42);
}

int main(int argc, char* argv[])
{
  if(argc==4) { drawcompeff(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}
