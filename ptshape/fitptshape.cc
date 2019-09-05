#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TF1.h>
#include <string>

#include "HEPData-ins1495026-v1/ppATLAS.h"

#include "xjjrootuti.h"
#include "xjjcuti.h"

void fitptshape(std::string inputname, std::string type, std::string inputdir="HEPData-ins1495026-v1")
{
  std::string outputname = xjjc::str_replaceall(xjjc::str_replaceall(inputname, ".root", ".pdf"), "rootfiles/", "plots/");
  gSystem->Exec(Form("mkdir -p %s", xjjc::str_replaceall(outputname, xjjc::str_divide(outputname, "/").back(), "").c_str()));

  TFile* inf = TFile::Open(inputname.c_str());
  TH1F* hmc = (TH1F*)inf->Get("hmc");
  TH1F* hmcweight = (TH1F*)inf->Get("hmcweight");
  ppref::ppATLAS getpp(inputdir);

  TGraphAsymmErrors* gg_stat = getpp.gg[type]["stat"];
  TGraphAsymmErrors* gg_syst = getpp.gg[type]["syst"];
  TH1F* hh_stat = getpp.hh[type]["stat"];
  TH1F* hh_syst = getpp.hh[type]["syst"];
  TH2F* hempty = getpp.hempty;

  xjjroot::sethempty(hempty, 0, 0);
  xjjroot::setthgrstyle(gg_stat, xjjroot::mycolor_middle["greenblue"], 21, 0.9, xjjroot::mycolor_middle["greenblue"], 1, 1);
  xjjroot::setthgrstyle(hh_stat, xjjroot::mycolor_middle["greenblue"], 25, 0.2, xjjroot::mycolor_middle["greenblue"], 1, 1);
  xjjroot::setthgrstyle(gg_syst, xjjroot::mycolor_middle["greenblue"], 21, 1, xjjroot::mycolor_middle["greenblue"], 1, 1, xjjroot::mycolor_middle["greenblue"], 0.4, 1001);
  xjjroot::setthgrstyle(hmc, xjjroot::mycolor_middle["azure"], 20, 0.8, xjjroot::mycolor_middle["azure"], 1, 2);
  xjjroot::setthgrstyle(hmcweight, xjjroot::mycolor_middle["red"], 47, 1.1, xjjroot::mycolor_middle["red"], 1, 2);

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  gPad->SetLogy();
  hempty->Draw();
  gg_syst->Draw("2 same");
  gg_stat->Draw("pe same");
  // hh_stat->Draw("p same");

  TF1* ff = new TF1("ff", "[0]+TMath::Exp([1]+[2]*x+[3]/x)", 10, 70);
  xjjroot::settfstyle(ff, xjjroot::mycolor_middle["greenblue"], 2, 3);
  ff->SetParLimits(0, 0, 1.);
  ff->SetParLimits(1, -1.e+3, 0);
  ff->SetParLimits(2, -1., 0);
  ff->SetParLimits(3, 0., 1.e+2);
  ff->SetParLimits(4, 0., 1.e+2);
  gg_stat->Fit(ff, "Q", "", 10, 50);
  gg_stat->Fit(ff, "Q", "", 10, 50);
  gg_stat->Fit(ff, "Q", "", 10, 50);
  // TF1* fffix = new TF1("fffix", getpp.formulanum[type].c_str(), 15, 70);
  // xjjroot::settfstyle(fffix, xjjroot::mycolor_middle["greenblue"], 2, 3);
  // fffix->Draw("same");

  //
  float boundmin = 15, boundmax = 50;
  float normint = hmc->Integral(hmc->GetXaxis()->FindBin(boundmin), hmc->GetXaxis()->FindBin(boundmax)-1);
  float norm = ff->Integral(boundmin, boundmax)/normint;
  hmc->Scale(norm, "width");
  hmcweight->Scale(norm, "width");
  //
  hmc->Draw("hist e same");
  TF1* ffmc = new TF1("ffmc", "[0]+TMath::Exp([1]+[2]*x+[3]/x+[4]/(x*x))", 15, 70);
  xjjroot::settfstyle(ffmc, xjjroot::mycolor_middle["blue"], 2, 3);
  ffmc->SetParLimits(0, 0, 1.);
  ffmc->SetParLimits(1, -1.e+1, 0);
  ffmc->SetParLimits(2, -1., 0);
  ffmc->SetParLimits(3, 0., 1.e+2); 
  ffmc->SetParLimits(4, -5.e+2, 0); 
  hmc->Fit(ffmc, "IQ", "", 15, 70);
  hmc->Fit(ffmc, "IQL", "", 15, 70);
  hmc->Fit(ffmc, "IQL", "", 15, 70);
  ffmc->Draw("same");
  hmcweight->Draw("p same");
  // TF1* ffmcfix = new TF1("ffmcfix", getpp.formuladen[type].c_str(), 15, 70);
  // xjjroot::settfstyle(ffmcfix, xjjroot::mycolor_middle["blue"], 2, 3);
  // ffmcfix->Draw("same");
  xjjroot::drawtex(0.88, 0.86, "PYTHIA", 0.04, 33, 62);
  TLegend* leg = new TLegend(0.60, 0.87-0.043-0.045*3, 0.88, 0.86-0.043);
  xjjroot::setleg(leg, 0.04);
  leg->AddEntry(gg_stat, "ATLAS 8TeV pp", "p");
  leg->AddEntry(hmc, "MC", "pl");
  leg->AddEntry(hmcweight, "MC weight p_{T}", "pl");
  leg->Draw();
  xjjroot::drawCMS("Simulation");
  xjjroot::drawtex(0.23, 0.86, getpp.tpar[type].c_str(), 0.04, 13, 62);
  std::cout<<Form("(%f+TMath::Exp(%f+%f*x+%f/x))/(%f+TMath::Exp(%f+%f*x+%f/x+%f/(x*x)))", ff->GetParameter(0), ff->GetParameter(1), ff->GetParameter(2), ff->GetParameter(3),
                  ffmc->GetParameter(0), ffmc->GetParameter(1), ffmc->GetParameter(2), ffmc->GetParameter(3), ffmc->GetParameter(4))<<std::endl;
  c->SaveAs(outputname.c_str());

}

int main(int argc, char* argv[])
{
  if(argc==4) { fitptshape(argv[1], argv[2], argv[3]); return 0; }
  if(argc==3) { fitptshape(argv[1], argv[2]); return 0; }
  return 1;
}
