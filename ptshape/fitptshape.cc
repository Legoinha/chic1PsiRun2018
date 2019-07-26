#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TF1.h>
#include <string>

#include "HEPData-ins1495026-v1-csv/ppATLAS.h"

#include "xjjrootuti.h"
#include "xjjcuti.h"

void fitptshape(std::string inputname, std::string type, std::string inputdir="HEPData-ins1495026-v1-csv")
{
  std::string outputname = xjjc::str_replaceall(xjjc::str_replaceall(inputname, ".root", ".pdf"), "rootfiles/", "plots/");
  gSystem->Exec(Form("mkdir -p %s", xjjc::str_replaceall(outputname, xjjc::str_divide(outputname, "/").back(), "").c_str()));

  TFile* inf = TFile::Open(inputname.c_str());
  TH1F* hmc = (TH1F*)inf->Get("hmc");
  ppRef::ppATLAS getpp(inputdir);

  TGraphAsymmErrors* gg_stat = getpp.gg[type]["stat"];
  TGraphAsymmErrors* gg_syst = getpp.gg[type]["syst"];
  TH1F* hh_stat = getpp.hh[type]["stat"];
  TH1F* hh_syst = getpp.hh[type]["syst"];
  TH2F* hempty = getpp.hempty;

  // hmc->Scale(hh_stat->Integral("width")/hmc->Integral(), "width");

  xjjroot::sethempty(hempty, 0, 0);
  xjjroot::setthgrstyle(gg_stat, xjjroot::mycolor_middle["red"], 21, 1, xjjroot::mycolor_middle["red"], 1, 1);
  xjjroot::setthgrstyle(hh_stat, xjjroot::mycolor_middle["red"], 25, 0.2, xjjroot::mycolor_middle["red"], 1, 1);
  xjjroot::setthgrstyle(gg_syst, xjjroot::mycolor_middle["red"], 21, 1, xjjroot::mycolor_middle["red"], 1, 1, xjjroot::mycolor_middle["red"], 0.4, 1001);
  xjjroot::setthgrstyle(hmc, xjjroot::mycolor_middle["azure"], 20, 0.8, xjjroot::mycolor_middle["azure"], 1, 2);

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  // gPad->SetLogx();
  gPad->SetLogy();
  hempty->Draw();
  gg_syst->Draw("2 same");
  gg_stat->Draw("pe same");
  hh_stat->Draw("pX0 same");

  // TF1* ff = new TF1("ff", "[0]+TMath::Exp([1]+[2]*x+[3]/x+[4]/(x*x))", 10, 70);
  TF1* ff = new TF1("ff", "[0]+TMath::Exp([1]+[2]*x+[3]/x)", 10, 70);
  xjjroot::settfstyle(ff, xjjroot::mycolor_middle["red"], 2, 2);
  // 
  // ff->SetParLimits(0, 0, 1.);
  // ff->SetParLimits(1, -1.e+3, 0);
  // ff->SetParLimits(2, -1., 0);
  // hh_stat->Fit(ff, "IQ", "", 10, 70);
  // hh_stat->Fit(ff, "IQL", "", 10, 70);
  // hh_stat->Fit(ff, "IL", "", 10, 70);
  // gg fit
  ff->SetParLimits(0, 0, 1.);
  ff->SetParLimits(1, -1.e+3, 0);
  ff->SetParLimits(2, -1., 0);
  ff->SetParLimits(3, 0., 1.e+2);
  ff->SetParLimits(4, 0., 1.e+2);
  gg_stat->Fit(ff, "Q", "", 10, 50);
  gg_stat->Fit(ff, "Q", "", 10, 50);
  gg_stat->Fit(ff, "", "", 10, 50);
  // ff->Draw("same");
  float boundmin = xjjc::str_contains(type, "promptX")?15:10, boundmax = 50;
  float normint = hmc->Integral(hmc->GetXaxis()->FindBin(boundmin), hmc->GetXaxis()->FindBin(boundmax));
  // integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin)/
  //   axis->GetBinWidth(bmin);
  // integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/
  // axis->GetBinWidth(bmax);
  hmc->Scale(ff->Integral(boundmin, boundmax)/normint, "width");
  hmc->Draw("hist e same");
  xjjroot::drawCMS("Simulation");
  c->SaveAs(outputname.c_str());
}

int main(int argc, char* argv[])
{
  if(argc==4) { fitptshape(argv[1], argv[2], argv[3]); return 0; }
  if(argc==3) { fitptshape(argv[1], argv[2]); return 0; }
  return 1;
}
