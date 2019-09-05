#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitX.h"

void drawkinematics();
void closure_drawhist(std::string input, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = TFile::Open(input.c_str());
  fitX::init(inf);
  TH1F* hyield_binned_a = (TH1F*)inf->Get("hyield_binned_a");
  TH1F* hyield_unbinned_a = (TH1F*)inf->Get("hyield_unbinned_a");
  TEfficiency* greff_a = (TEfficiency*)inf->Get("greff_a");
  TH1F* heffgen_a = (TH1F*)inf->Get("heffgen_a");
  TH1F* hcorr_binned_a = (TH1F*)inf->Get("hcorr_binned_a");
  TH1F* hcorr_unbinned_a = (TH1F*)inf->Get("hcorr_unbinned_a");
  TH1F* hyield_binned_b = (TH1F*)inf->Get("hyield_binned_b");
  TH1F* hyield_unbinned_b = (TH1F*)inf->Get("hyield_unbinned_b");
  TEfficiency* greff_b = (TEfficiency*)inf->Get("greff_b");
  TH1F* heffgen_b = (TH1F*)inf->Get("heffgen_b");
  TH1F* hcorr_binned_b = (TH1F*)inf->Get("hcorr_binned_b");
  TH1F* hcorr_unbinned_b = (TH1F*)inf->Get("hcorr_unbinned_b");

  int nbins = heffgen_a->GetXaxis()->GetNbins();
  std::vector<float> ptbins;
  for(int i=0; i<nbins+1; i++)
    ptbins.push_back((*(heffgen_a->GetXaxis()->GetXbins()))[i]);

  xjjroot::setthgrstyle(hyield_unbinned_a, xjjroot::mycolor_dark["red"], 47, 1.6, xjjroot::mycolor_dark["red"], 1, 2);
  xjjroot::setthgrstyle(hyield_binned_a, xjjroot::mycolor_light["red"], 24, 1.1, xjjroot::mycolor_light["red"], 1, 2);
  xjjroot::setthgrstyle(hyield_unbinned_b, xjjroot::mycolor_dark["red"], 47, 1.6, xjjroot::mycolor_dark["red"], 1, 2);
  xjjroot::setthgrstyle(hyield_binned_b, xjjroot::mycolor_light["red"], 24, 1.1, xjjroot::mycolor_light["red"], 1, 2);
  xjjroot::sethempty(hyield_unbinned_a);
  xjjroot::sethempty(hyield_binned_a);
  xjjroot::sethempty(hyield_unbinned_b);
  xjjroot::sethempty(hyield_binned_b);
  xjjroot::setthgrstyle(hcorr_unbinned_a, xjjroot::mycolor_dark["red"], 47, 1.6, xjjroot::mycolor_dark["red"], 1, 2);
  xjjroot::setthgrstyle(hcorr_binned_a, xjjroot::mycolor_light["red"], 24, 1.1, xjjroot::mycolor_light["red"], 1, 2);
  xjjroot::setthgrstyle(hcorr_unbinned_b, xjjroot::mycolor_dark["red"], 47, 1.6, xjjroot::mycolor_dark["red"], 1, 2);
  xjjroot::setthgrstyle(hcorr_binned_b, xjjroot::mycolor_light["red"], 24, 1.1, xjjroot::mycolor_light["red"], 1, 2);
  xjjroot::sethempty(hcorr_unbinned_a);
  xjjroot::sethempty(hcorr_binned_a);
  xjjroot::sethempty(hcorr_unbinned_b);
  xjjroot::sethempty(hcorr_binned_b);

  xjjroot::setthgrstyle(heffgen_a, xjjroot::mycolor_middle["azure"], 21, 1, xjjroot::mycolor_middle["azure"], 1, 2, xjjroot::mycolor_middle["azure"], 0.1, 1001);
  xjjroot::sethempty(heffgen_a);
  xjjroot::setthgrstyle(heffgen_b, xjjroot::mycolor_middle["azure"], 21, 1, xjjroot::mycolor_middle["azure"], 1, 2, xjjroot::mycolor_middle["azure"], 0.1, 1001);
  xjjroot::sethempty(heffgen_b);

  TLegend* leg = new TLegend(0.61, 0.70-0.045*4, 0.86, 0.70);
  xjjroot::setleg(leg, 0.04);
  leg->AddEntry(heffgen_a, "gen signal", "lf");
  leg->AddEntry((TObject*)0, "corrected signal", NULL);
  leg->AddEntry(hcorr_unbinned_a, "unbinned fits", "pl");
  leg->AddEntry(hcorr_binned_b, "binned fits", "pl");
  TLegend* legy = new TLegend(0.61, 0.70-0.045*2, 0.86, 0.70);
  xjjroot::setleg(legy, 0.04);
  legy->AddEntry(hyield_unbinned_a, "unbinned fits", "pl");
  legy->AddEntry(hyield_binned_b, "binned fits", "pl");

  float ymaxeff = 0.2;
  TH2F* hemptyeff = new TH2F("hemptyeff", ";p_{T} (GeV/c);#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 10, ptbins.front(), ptbins.back(), 10, 0, ymaxeff);
  xjjroot::sethempty(hemptyeff, 0, 0.3);
  xjjroot::setgstyle(1);
  TCanvas* c_a = new TCanvas("c_a", "", 1800, 600);
  c_a->Divide(3, 1);
  c_a->cd(1);
  hyield_binned_a->Draw("pe");
  hyield_unbinned_a->Draw("pe same");
  legy->Draw();
  drawkinematics();
  xjjroot::drawtex(0.24, 0.30, Form("Prompt %s", fitX::title_a.c_str()), 0.038, 12, 62);
  xjjroot::drawCMS("Simulation");
  c_a->cd(2);
  hemptyeff->Draw();
  greff_a->Draw("same3");
  greff_a->Draw("samelX");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.30, Form("Prompt %s", fitX::title_a.c_str()), 0.038, 12, 62);
  xjjroot::drawCMS("Simulation");
  c_a->cd(3);
  gPad->SetLogy();
  heffgen_a->Draw("hist");
  hcorr_binned_a->Draw("pe same");
  hcorr_unbinned_a->Draw("pe same");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.30, Form("Prompt %s", fitX::title_a.c_str()), 0.038, 12, 62);
  xjjroot::drawCMS("Simulation");
  leg->Draw();
  std::string outputname_a = "plots/"+output+"/cclosure_a.pdf";
  xjjroot::mkdir(outputname_a);
  c_a->SaveAs(outputname_a.c_str());

  TCanvas* c_b = new TCanvas("c_b", "", 1800, 600);
  c_b->Divide(3, 1);
  c_b->cd(1);
  hyield_binned_b->Draw("pe");
  hyield_unbinned_b->Draw("pe same");
  legy->Draw();
  drawkinematics();
  xjjroot::drawtex(0.24, 0.30, Form("Prompt %s", fitX::title_b.c_str()), 0.038, 12, 62);
  xjjroot::drawCMS("Simulation");
  c_b->cd(2);
  hemptyeff->Draw();
  greff_b->Draw("same3");
  greff_b->Draw("samelX");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.30, Form("Prompt %s", fitX::title_b.c_str()), 0.038, 12, 62);
  xjjroot::drawCMS("Simulation");
  c_b->cd(3);
  gPad->SetLogy();
  heffgen_b->Draw("hist");
  hcorr_binned_b->Draw("pe same");
  hcorr_unbinned_b->Draw("pe same");
  leg->Draw();
  drawkinematics();
  xjjroot::drawtex(0.24, 0.30, Form("Prompt %s", fitX::title_b.c_str()), 0.038, 12, 62);
  xjjroot::drawCMS("Simulation");

  std::string outputname_b = "plots/"+output+"/cclosure_b.pdf";
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());

}

int main(int argc, char* argv[])
{
  if(argc==3) { closure_drawhist(argv[1], argv[2]); return 0; }
  return 1;
}

void drawkinematics()
{
  xjjroot::drawtex(0.90, 0.84, "PYTHIA + HYDJET", 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 42);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 42);
}

