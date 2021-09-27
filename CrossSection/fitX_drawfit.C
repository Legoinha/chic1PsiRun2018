#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "fit.h"

namespace drawfit
{
  TF1 *f_dump_incl = new TF1(), *bkg_dump_incl = new TF1(), *f_dump_Benr = new TF1(), *bkg_dump_Benr = new TF1();
  TH1F* h_dump = new TH1F("h_dump", "", 1, 0, 1);
  void setlegs(TLegend* leg_both, Size_t tsize);
}

void fitX_drawfit(std::string input, std::string output)
{
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);

  RooWorkspace* ww = (RooWorkspace*)inf->Get("wfit");
  RooDataSet* dsh = (RooDataSet*)ww->data("dsh");
  RooDataSet* dshBenr = (RooDataSet*)ww->data("dshBenr");
  RooAbsPdf* pdf = ww->pdf("pdf_th");
  RooAbsPdf* pdfBenr = ww->pdf("pdf_thBenr");
  RooAbsPdf* bkg = ww->pdf("bkg_th");
  RooAbsPdf* bkgBenr = ww->pdf("bkg_thBenr");
  TF1* fr = (TF1*)inf->Get("fr_th");
  TF1* frBenr = (TF1*)inf->Get("fr_thBenr");
  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);

  xjjroot::setgstyle(2);

  /***************************************************
    Paper figure (1)
  ***************************************************/
  xjjroot::mkdir("plots/"+output+"/paper/x");
  float yincl = 0.49, ybenr = 1-yincl;
  Size_t tsize = 0.07; // 0.05
  float axissize = 1.7;
  RooPlot* frempty_incl = mass->frame(RooFit::Title(""));
  frempty_incl->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  fitX::setmasshist(frempty_incl, 0, -0.65);
  frempty_incl->GetXaxis()->SetTitleSize(frempty_incl->GetXaxis()->GetTitleSize()*axissize);
  frempty_incl->GetYaxis()->SetTitleSize(frempty_incl->GetYaxis()->GetTitleSize()*axissize);
  frempty_incl->GetXaxis()->SetLabelSize(frempty_incl->GetXaxis()->GetLabelSize()*axissize);
  frempty_incl->GetYaxis()->SetLabelSize(frempty_incl->GetYaxis()->GetLabelSize()*axissize); // h->GetYaxis()->SetLabelSize(frempty_incl->GetYaxis()->GetLabelSize());
  RooPlot* frempty_benr = mass->frame(RooFit::Title(""));
  frempty_benr->SetXTitle("m#scale[0.7]{#mu#mu#pi#pi} (GeV/c^{2})");
  frempty_benr->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  fitX::setmasshist(frempty_benr, -0.2, -0.65);
  frempty_benr->GetXaxis()->SetTitleSize(frempty_benr->GetXaxis()->GetTitleSize()*axissize*1.1);
  frempty_benr->GetYaxis()->SetTitleSize(frempty_benr->GetYaxis()->GetTitleSize()*axissize);
  frempty_benr->GetXaxis()->SetLabelSize(frempty_benr->GetXaxis()->GetLabelSize()*axissize);
  frempty_benr->GetYaxis()->SetLabelSize(frempty_benr->GetYaxis()->GetLabelSize()*axissize); // hBenr->GetYaxis()->SetLabelSize(frempty_benr->GetYaxis()->GetLabelSize());
  TLegend* leg_both = new TLegend(0.66, 0.95-0.072*3, 0.91, 0.95);
  drawfit::setlegs(leg_both, tsize);
  TCanvas* cr = new TCanvas("cr", "", 700, 700);
  TPad* pincl = new TPad("pincl", "", 0, 1-yincl, 1, 1);
  pincl->SetMargin(xjjroot::margin_pad_left*0.62/*0.5*/, xjjroot::margin_pad_right*1.8/*1.5*/, 0, xjjroot::margin_pad_top*1.3);
  pincl->Draw();
  pincl->cd();
  dsh->plotOn(frempty_incl, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(drawfit::h_dump->GetMarkerSize()), RooFit::MarkerStyle(drawfit::h_dump->GetMarkerStyle()), RooFit::LineColor(drawfit::h_dump->GetLineColor()), RooFit::LineWidth(drawfit::h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  pdf->plotOn(frempty_incl, RooFit::Name("bkg"), RooFit::Components(*bkg), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(drawfit::bkg_dump_incl->GetLineStyle()), RooFit::LineColor(drawfit::bkg_dump_incl->GetLineColor()), RooFit::LineWidth(drawfit::bkg_dump_incl->GetLineWidth()));
  pdf->plotOn(frempty_incl, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(drawfit::f_dump_incl->GetLineStyle()), RooFit::LineColor(drawfit::f_dump_incl->GetLineColor()), RooFit::LineWidth(drawfit::f_dump_incl->GetLineWidth()));
  dsh->plotOn(frempty_incl, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(drawfit::h_dump->GetMarkerSize()), RooFit::MarkerStyle(drawfit::h_dump->GetMarkerStyle()), RooFit::LineColor(drawfit::h_dump->GetLineColor()), RooFit::LineWidth(drawfit::h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  frempty_incl->SetMinimum(130./fitX::NBIN*38.);
  frempty_incl->SetMaximum(440./fitX::NBIN*38.);
  frempty_incl->Draw();
  TH1F* hforpull_incl = fitX::createhistforpull(frempty_incl, "dshist", fr, "_incl");
  fitX::drawpull(hforpull_incl, fr, fitX::color_data);
  // xjjroot::drawbox(4.02, h->GetMinimum(), 5, h->GetMaximum(), kWhite, 1);
  xjjroot::drawtex(0.99, 0.39, "Pull", frempty_incl->GetYaxis()->GetTitleSize(), 33, 42, fitX::color_data, 270);
  xjjroot::drawtex(0.88, 0.78/*0.79*/, Form("%.0f < p_{T} < %.0f GeV/c", fitX::ptmincut, fitX::ptmaxcut), tsize, 32, 42/*62*/);
  xjjroot::drawtex(0.88, 0.78-0.072 /*0.06*/, Form("%s|y| < %.1f", (fitX::ymincut?Form("%.1f < ", fitX::ymincut):""), fitX::ymaxcut), tsize, 32, 42);
  xjjroot::drawtex(0.88, 0.78-0.072*2, Form("Cent. %.0f-%.0f%s", fitX::centmincut, fitX::centmaxcut, "%"), tsize, 32, 42);
  xjjroot::drawtex(0.22, 0.13, fitX::title_a.c_str(), tsize);
  xjjroot::drawtex(0.59, 0.32, fitX::title_b.c_str(), tsize);
  xjjroot::drawtex(0.15/*0.14*/, 0.79, "#scale[1.25]{#bf{CMS}} #it{Internal}", tsize, 12); // 0.06
  // xjjroot::drawtex(0.15/*0.14*/, 0.79, "#scale[1.25]{#bf{CMS}}", tsize, 12); // 0.06
  xjjroot::drawtex(0.92, 0.92, "1.7 nb^{-1} (PbPb 5.02 TeV)", tsize, 32); // 0.055
  xjjroot::drawtex(0.15, 0.70, "Inclusive", tsize, 12, 52, kBlack);
  cr->cd();

  TPad* pbenr = new TPad("pbenr", "", 0, 0, 1, ybenr);
  pbenr->SetMargin(xjjroot::margin_pad_left*0.62/*0.5*/, xjjroot::margin_pad_right*1.8/*1.5*/, xjjroot::margin_pad_bottom*1.3/*1*/, 0);
  pbenr->Draw();
  pbenr->cd();
  dshBenr->plotOn(frempty_benr, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(drawfit::h_dump->GetMarkerSize()), RooFit::MarkerStyle(drawfit::h_dump->GetMarkerStyle()), RooFit::LineColor(drawfit::h_dump->GetLineColor()), RooFit::LineWidth(drawfit::h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  pdfBenr->plotOn(frempty_benr, RooFit::Name("bkg"), RooFit::Components(*bkgBenr), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(drawfit::bkg_dump_Benr->GetLineStyle()), RooFit::LineColor(drawfit::bkg_dump_Benr->GetLineColor()), RooFit::LineWidth(drawfit::bkg_dump_Benr->GetLineWidth()));
  pdfBenr->plotOn(frempty_benr, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(drawfit::f_dump_Benr->GetLineStyle()), RooFit::LineColor(drawfit::f_dump_Benr->GetLineColor()), RooFit::LineWidth(drawfit::f_dump_Benr->GetLineWidth()));
  dshBenr->plotOn(frempty_benr, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(drawfit::h_dump->GetMarkerSize()), RooFit::MarkerStyle(drawfit::h_dump->GetMarkerStyle()), RooFit::LineColor(drawfit::h_dump->GetLineColor()), RooFit::LineWidth(drawfit::h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  frempty_benr->SetMinimum(0); // hBenr->SetMinimum(0);
  frempty_benr->SetMaximum(75./fitX::NBIN*38.); // hBenr->SetMaximum(75./fitX::NBIN*38.); // maxy[1]
  frempty_benr->Draw();
  leg_both->Draw();
  TH1F* hforpull_benr = fitX::createhistforpull(frempty_benr, "dshist", frBenr, "_benr");
  // hforpull_benr->GetXaxis()->SetTitle("m#scale[0.8]{#mu#mu#pi#pi}");
  fitX::drawpull(hforpull_benr, frBenr, fitX::color_data2);
  // fitX::drawpull(hBenr, frBenr, fitX::color_data2);
  // xjjroot::drawbox(4.02, hBenr->GetMinimum(), 5, hBenr->GetMaximum(), kWhite, 1);
  xjjroot::drawtex(0.99, 0.53, "Pull", frempty_incl->GetYaxis()->GetTitleSize(), 33, 42, fitX::color_data2, 270);
  xjjroot::drawtex(0.15, 0.90, "b-enriched (l#scale[0.8]{xy} > 0.1 mm)", tsize, 12, 52, kBlack);
  cr->SaveAs(Form("plots/%s/paper/chmassr_both_2p.pdf", output.c_str()));

}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX_drawfit(argv[1], argv[2]); return 0; }
  return 1;
}

void drawfit::setlegs(TLegend* leg_both, Size_t tsize)
{
  xjjroot::setthgrstyle(drawfit::h_dump, kBlack, 20, 0.9, kBlack, 1, 1);
  xjjroot::settfstyle(drawfit::f_dump_incl, fitX::color_data, 1, 2);
  xjjroot::settfstyle(drawfit::f_dump_Benr, fitX::color_data2, 1, 2);
  xjjroot::settfstyle(drawfit::bkg_dump_incl, fitX::color_data, 7, 2);
  xjjroot::settfstyle(drawfit::bkg_dump_Benr, fitX::color_data2, 7, 2);

  xjjroot::setleg(leg_both, tsize); // 0.055
  leg_both->SetNColumns(2);
  leg_both->AddEntry((TObject*)0, "", NULL);
  leg_both->AddEntry(h_dump, "data", "pe");
  leg_both->AddEntry(f_dump_incl, "", "l");
  leg_both->AddEntry(f_dump_Benr, "total fit", "l");
  leg_both->AddEntry(bkg_dump_incl, "", "l");
  leg_both->AddEntry(bkg_dump_Benr, "background", "l");
}
