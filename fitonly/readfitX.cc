#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooMinuit.h>
#include <RooAddPdf.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLegend.h>

#include <string>
#include <vector>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "fit.h"
#include "xnll.h"

void drawnll(std::vector<TH1F*> hnlls, std::string name, float ymin, float ymax, std::string outputdir);
TLegend* leg;
void readfitX(std::string outputdir, std::string fitopts="default")
{
  if(xjjc::str_contains(outputdir, "narrow")) { fitX::NBIN = 30; fitX::BIN_MAX = 3.92; }
  std::vector<std::string> options = xjjc::str_divide(fitopts, ",");
  int nopt = options.size();
  std::string outputname = "rootfiles/" + outputdir + "/rootnll.root";
  TFile* inf = TFile::Open(outputname.c_str());
  RooWorkspace* wwnll = (RooWorkspace*)inf->Get("wwnll");
  RooRealVar* par10 = wwnll->var("par10");
  RooAbsData* data = wwnll->data("dsh_ws");
  RooRealVar* mass = wwnll->var("Bmass");

  std::vector<fitX::xnll*> xnlls(nopt, 0);
  for(int i=0; i<nopt; i++)
    {
      std::cout<<std::endl<<"\e[44;1m"<<options[i]<<": read files\e[0m"<<std::endl;
      std::string thisoutput = "rootfiles/" + outputdir + "/" + options[i] +"/rootnll.root";
      TFile* thisinf = TFile::Open(thisoutput.c_str());
      RooWorkspace* ww = (RooWorkspace*)thisinf->Get("wwnll");
      xnlls[i] = new fitX::xnll(ww, options[i]);
    }

  xjjroot::setgstyle();

  // cnllscan
  TCanvas* cnll = new TCanvas("cnll", "", 600, 600);
  RooPlot* framenll = par10->frame(RooFit::Bins(20), RooFit::Range(0., 200.));
  xjjroot::sethempty(framenll);
  TLegend* leg_cnll = new TLegend(0.65, 0.86-nopt*0.04, 0.90, 0.86);
  xjjroot::setleg(leg_cnll, 0.035);
  for(int i=0; i<nopt; i++)
    {
      std::cout<<std::endl<<"\e[44;1m"<<options[i]<<": cnll\e[0m"<<std::endl;
      Color_t cc = xjjroot::colorlist_middle[i];
      Style_t ss = 9;
      if(options[i]=="pol4") { ss = 1; cc = kBlack; }
      xnlls[i]->nll()->plotOn(framenll, RooFit::Name(Form("nll_%s", options[i].c_str())), RooFit::LineColor(cc), RooFit::LineWidth(3), RooFit::LineStyle(ss), RooFit::Normalization(2., RooAbsReal::Raw), RooFit::ShiftToZero());
      leg_cnll->AddEntry(framenll->findObject(Form("nll_%s", options[i].c_str())), options[i].c_str(), "l");
    }
  cnll->cd();
  framenll->SetMinimum(0.);
  framenll->SetMaximum(20.);
  framenll->GetXaxis()->SetTitle("N_{sig}^{(X3872)}");
  framenll->GetYaxis()->SetTitle("Projection of -2log(LL)");
  framenll->Draw();
  leg_cnll->Draw("same");
  xjjroot::drawCMSleft("", 0.05, -0.08);
  xjjroot::drawCMSright("1.7 nb^{-1} (2018 PbPb 5.02 TeV)");
  std::string outputnamenll = "plots/" + outputdir + "/cnllscan_fitopt.pdf";
  xjjroot::mkdir(outputnamenll);
  cnll->SaveAs(outputnamenll.c_str());

  TCanvas* cpll = new TCanvas("cpll", "", 600, 600);
  RooPlot* framepll = par10->frame(RooFit::Bins(20), RooFit::Range(0., 200.));
  xjjroot::sethempty(framepll);
  TLegend* leg_cpll = new TLegend(0.65, 0.86-nopt*0.04, 0.90, 0.86);
  xjjroot::setleg(leg_cpll, 0.035);
  for(int i=0; i<nopt; i++)
    {
      std::cout<<std::endl<<"\e[44;1m"<<options[i]<<": cpll\e[0m"<<std::endl;
      Color_t cc = xjjroot::colorlist_middle[i];
      Style_t ss = 9;
      if(options[i]=="pol4") { ss = 1; cc = kBlack; }
      xnlls[i]->pll()->plotOn(framepll, RooFit::Name(Form("pll_%s", options[i].c_str())), RooFit::LineColor(cc), RooFit::LineStyle(ss), RooFit::LineWidth(3), RooFit::Normalization(2., RooAbsReal::Raw));
      leg_cpll->AddEntry(framepll->findObject(Form("pll_%s", options[i].c_str())), options[i].c_str(), "l");
    }
  cpll->cd();
  framepll->GetXaxis()->SetTitle("N_{sig}^{(X3872)}");
  framepll->GetYaxis()->SetTitle("Projection of -2log(profileLL)");
  framepll->Draw();
  leg_cpll->Draw("same");
  xjjroot::drawCMSleft("", 0.05, -0.08);
  xjjroot::drawCMSright("1.7 nb^{-1} (2018 PbPb 5.02 TeV)");
  std::string outputnamepll = "plots/" + outputdir + "/cpllscan_fitopt.pdf";
  xjjroot::mkdir(outputnamepll);
  cpll->SaveAs(outputnamepll.c_str());

  for(int i=0; i<nopt; i++)
    {
      std::cout<<std::endl<<"\e[44;1m"<<options[i]<<": cnpll\e[0m"<<std::endl;
      TCanvas* cnpll = new TCanvas("cnpll", "", 600, 600);
      RooPlot* framenpll = par10->frame(RooFit::Bins(20), RooFit::Range(0., 200.));
      xjjroot::sethempty(framenpll);
      TLegend* leg_cnpll = new TLegend(0.65, 0.86-2*0.04, 0.90, 0.86);
      xjjroot::setleg(leg_cnpll, 0.035);
      xnlls[i]->nll()->plotOn(framenpll, RooFit::Name(Form("np_nll_%s", options[i].c_str())), RooFit::LineColor(xjjroot::mycolor_satmiddle["azure"]), RooFit::LineStyle(1), RooFit::LineWidth(3), RooFit::Normalization(2., RooAbsReal::Raw), RooFit::ShiftToZero()); // x2
      leg_cnpll->AddEntry(framenpll->findObject(Form("np_nll_%s", options[i].c_str())), "LL", "l");
      xnlls[i]->pll()->plotOn(framenpll, RooFit::Name(Form("np_pll_%s", options[i].c_str())), RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::LineWidth(3), RooFit::Normalization(2., RooAbsReal::Raw)); // x2
      leg_cnpll->AddEntry(framenpll->findObject(Form("np_pll_%s", options[i].c_str())), "profileLL", "l");
      cnpll->cd();
      framenpll->SetMinimum(0);
      framenpll->SetMaximum(20.);
      framenpll->GetXaxis()->SetTitle("N_{sig}^{(X3872)}");
      framenpll->GetYaxis()->SetTitle("Projection of -2log(LL)");
      framenpll->Draw();
      leg_cnpll->Draw("same");
      // xjjroot::drawline(0, 0.50*2, 200, 0.50*2, xjjroot::mycolor_middle["red"], 2, 3);
      // xjjroot::drawline(0, 1.96*2, 200, 1.96*2, xjjroot::mycolor_middle["red"], 2, 3);
      // xjjroot::drawline(xnlls[i]->x68_d(), 0, xnlls[i]->x68_d(), 0.50*2, xjjroot::mycolor_middle["red"], 2, 3);
      // xjjroot::drawline(xnlls[i]->x68_u(), 0, xnlls[i]->x68_u(), 0.50*2, xjjroot::mycolor_middle["red"], 2, 3);
      // xjjroot::drawline(xnlls[i]->x95_d(), 0, xnlls[i]->x95_d(), 1.96*2, xjjroot::mycolor_middle["red"], 2, 3);
      // xjjroot::drawline(xnlls[i]->x95_u(), 0, xnlls[i]->x95_u(), 1.96*2, xjjroot::mycolor_middle["red"], 2, 3);
      xjjroot::drawCMSleft("", 0.05, -0.08);
      xjjroot::drawCMSright("1.7 nb^{-1} (2018 PbPb 5.02 TeV)");
      std::string outputnamenpll = "plots/" + outputdir + "/" + options[i] + "/cnpllscan.pdf";
      xjjroot::mkdir(outputnamenpll);
      cnpll->SaveAs(outputnamenpll.c_str());
      delete cnpll;
    }

  leg = new TLegend(0.65, 0.86-nopt*0.04, 0.90, 0.86);
  xjjroot::setleg(leg, 0.035);
  std::vector<TH1F*> hnlls(nopt, 0), hnllscorr1(nopt, 0), hnllscorr2(nopt, 0), hplls(nopt, 0);
  float ymin = 1.e+10, ymincorr1 = 1.e+10, ymincorr2 = 1.e+10;
  float ymax = -1.e+10, ymaxcorr1 = -1.e+10, ymaxcorr2 = -1.e+10;
  std::vector<float> aymin(nopt, 1.e+10), aymax(nopt, -1.e+10);
  for(int i=0; i<nopt; i++)
    {
      std::cout<<std::endl<<"\e[44;1m"<<options[i]<<": createHistogram\e[0m"<<std::endl;
      xnlls[i]->createHist();

      xjjroot::setthgrstyle(xnlls[i]->hnll(), xjjroot::colorlist_middle[i], 20, 0.5, xjjroot::colorlist_middle[i], 9, 3);
      if(options[i]=="pol4") { xnlls[i]->hnll()->SetLineStyle(1); xnlls[i]->hnll()->SetLineColor(kBlack); }
      leg->AddEntry(xnlls[i]->hnll(), options[i].c_str(), "l");
      xjjroot::setthgrstyle(xnlls[i]->hnllcorr1(), xnlls[i]->hnll()->GetMarkerColor(), 20, 0.5, xnlls[i]->hnll()->GetLineColor(), xnlls[i]->hnll()->GetLineStyle(), 3);
      xjjroot::setthgrstyle(xnlls[i]->hnllcorr2(), xnlls[i]->hnll()->GetMarkerColor(), 20, 0.5, xnlls[i]->hnll()->GetLineColor(), xnlls[i]->hnll()->GetLineStyle(), 3);
      aymin[i] = xnlls[i]->min_nll();
      aymax[i] = xnlls[i]->max_nll();

      hnlls[i] = xjjroot::copyobject(xnlls[i]->hnll(), Form("%s_copy", xnlls[i]->hnll()->GetName()));
      hnllscorr1[i] = xjjroot::copyobject(xnlls[i]->hnll(), Form("%s_copy", xnlls[i]->hnllcorr1()->GetName()));
      hnllscorr2[i] = xjjroot::copyobject(xnlls[i]->hnll(), Form("%s_copy", xnlls[i]->hnllcorr2()->GetName()));

      if(fabs(aymax[i]-ymax) > 1.e+4 && i!=0) continue;
      if(xnlls[i]->min_nll() < ymin) { ymin = xnlls[i]->min_nll(); }
      if(xnlls[i]->max_nll() > ymax) { ymax = xnlls[i]->max_nll(); }
      if(xnlls[i]->min_nllcorr1() < ymincorr1) { ymincorr1 = xnlls[i]->min_nllcorr1(); }
      if(xnlls[i]->max_nllcorr1() > ymaxcorr1) { ymaxcorr1 = xnlls[i]->max_nllcorr1(); }
      if(xnlls[i]->min_nllcorr2() < ymincorr2) { ymincorr2 = xnlls[i]->min_nllcorr2(); }
      if(xnlls[i]->max_nllcorr2() > ymaxcorr2) { ymaxcorr2 = xnlls[i]->max_nllcorr2(); }
    }
  // std::cout<<"\e[44;1m"<<ymin<<" "<<ymax<<"\e[0m"<<std::endl;
  // std::cout<<"\e[44;1m"<<ymincorr1<<" "<<ymaxcorr1<<"\e[0m"<<std::endl;
  // std::cout<<"\e[44;1m"<<ymincorr2<<" "<<ymaxcorr2<<"\e[0m"<<std::endl;
  drawnll(hnlls, "", ymin, ymax, outputdir);
  drawnll(hnllscorr1, "corr1", ymincorr1, ymaxcorr1, outputdir);
  drawnll(hnllscorr2, "corr2", ymincorr2, ymaxcorr2, outputdir);

  std::cout<<"\e[36;1m"<<std::endl;
  UInt_t stdwid = 15;
  std::cout<<std::string(stdwid*4+1, '-')<<std::endl;
  std::cout<<std::left<<std::setw(stdwid)<<"| Title"<<std::setw(stdwid)<<"| Sig (PLL)"<<std::setw(stdwid)<<"| Sig (Null)"<<std::setw(stdwid)<<"| Yield"<<"|"<<std::endl;
  std::cout<<std::string(stdwid*4+1, '-')<<std::endl;
  for(int i=0; i<nopt; i++)
    {
      std::cout<<std::left<<std::setw(stdwid)<<Form("| %s", options[i].c_str())<<std::setw(stdwid)<<Form("| %.4f", xnlls[i]->sig())<<std::setw(stdwid)<<Form("| %.4f", xnlls[i]->nullsig())<<std::setw(stdwid)<<Form("| %.1f", xnlls[i]->xmean_pll())<<"|"<<std::endl;
      std::cout<<std::string(stdwid*4+1, '-')<<std::endl;      
    }
  std::cout<<"\e[0m"<<std::endl;

  // fit
  // RooPlot* frempty = mass->frame(RooFit::Title(""));
  // xjjroot::sethempty(frempty);
  // frempty->SetXTitle("m_{#mu#mu#pi#pi} (GeV/c^{2})");
  // frempty->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  // frempty->SetMinimum(50);
  // frempty->SetMaximum(450);
  // TCanvas* cmass = new TCanvas("cmass", "", 700, 600);
  // data->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  // for(int i=0; i<nopt; i++)
  //   {
  //     pdfs[i]->plotOn(frempty, RooFit::Name("bkg"), RooFit::Components(*(bkgs[i])), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(xjjroot::colorlist_middle[i]), RooFit::LineWidth(2));
  //     pdfs[i]->plotOn(frempty, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(xjjroot::colorlist_middle[i]), RooFit::LineWidth(2));
  //   }
  // data->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  // frempty->Draw();
  // std::string outputnamemass = "plots/" + outputdir + "/cmassscan_fitopt.pdf";
  // xjjroot::mkdir(outputnamemass);
  // cmass->SaveAs(outputnamemass.c_str());
}

int main(int argc, char* argv[])
{
  if(argc==3) { readfitX(argv[1], argv[2]); return 0; }
}

void drawnll(std::vector<TH1F*> hnlls, std::string name, float ymin, float ymax, std::string outputdir)
{
  TCanvas* cnllNpar = new TCanvas("cnllNpar", "", 600, 600);
  std::string ytitle = "";
  if(name=="corr1") ytitle = " + N_{par}";
  if(name=="corr2") ytitle = " + 2N_{par}";
  TH2F* hframenll = new TH2F("hframenll", Form(";N_{sig}^{(X3872)};Projection of -2log(LL)%s", ytitle.c_str()), 10, 0, 200, 10, 0, 1.05*(ymax-ymin));
  // std::cout<<ymax-ymin<<std::endl;
  xjjroot::sethempty(hframenll);
  hframenll->Draw();
  for(int i=0; i<hnlls.size(); i++)
    {
      for(int j=0; j<hnlls[i]->GetXaxis()->GetNbins(); j++) { hnlls[i]->SetBinContent(j+1, hnlls[i]->GetBinContent(j+1)-ymin); }
      hnlls[i]->Draw("c same");
    }
  xjjroot::drawCMSleft("", 0.05, -0.08);
  xjjroot::drawCMSright("1.7 nb^{-1} (2018 PbPb 5.02 TeV)");
  leg->Draw();
  std::string outputnamenll = "plots/" + outputdir + "/cnllscanNpar_fitopt" + name + ".pdf";
  xjjroot::mkdir(outputnamenll);
  cnllNpar->SaveAs(outputnamenll.c_str());

  delete hframenll;
  delete cnllNpar;
}
