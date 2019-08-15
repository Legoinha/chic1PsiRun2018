#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "fitX.h"

void drawkinematic();

namespace pdfvar
{
  std::vector<std::string> fitopt = {"poly4", "poly3", "poly2", "cheb4", "cheb3", "cheb2"};
  int nfit = fitopt.size();
}

void fitX_pdfvar(std::string input, std::string output)
{
  TFile* inf = new TFile(Form("%s.root", input.c_str()));
  TH1F* h = (TH1F*)inf->Get("h");
  TH1F* hBenr = (TH1F*)inf->Get("hBenr");
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");

  std::vector<TH1F*> hyield(pdfvar::nfit), hBenryield(pdfvar::nfit);
  for(int ff=0; ff<pdfvar::nfit; ff++)
    {
      hyield[ff] = new TH1F(Form("hyield_%s", pdfvar::fitopt[ff].c_str()), "", 5, 0, 5);
      hBenryield[ff] = new TH1F(Form("hBenryield_%s", pdfvar::fitopt[ff].c_str()), "", 5, 0, 5);
    }
  TH1F* hyields_a = new TH1F("hyields_a", ";;Raw Yield", pdfvar::nfit, 0, pdfvar::nfit);
  xjjroot::sethempty(hyields_a, 0, 0.3);
  hyields_a->GetXaxis()->SetLabelSize(hyields_a->GetXaxis()->GetLabelSize());
  TH1F* hBenryields_a = new TH1F("hBenryields_a", ";;Raw Yield", pdfvar::nfit, 0, pdfvar::nfit);
  xjjroot::sethempty(hBenryields_a, 0, 0.3);
  hBenryields_a->GetXaxis()->SetLabelSize(hBenryields_a->GetXaxis()->GetLabelSize());
  TH1F* hyields_b = new TH1F("hyields_b", ";;Raw Yield", pdfvar::nfit, 0, pdfvar::nfit);
  xjjroot::sethempty(hyields_b, 0, 0.3);
  hyields_b->GetXaxis()->SetLabelSize(hyields_b->GetXaxis()->GetLabelSize());
  TH1F* hBenryields_b = new TH1F("hBenryields_b", ";;Raw Yield", pdfvar::nfit, 0, pdfvar::nfit);
  xjjroot::sethempty(hBenryields_b, 0, 0.3);
  hBenryields_b->GetXaxis()->SetLabelSize(hBenryields_b->GetXaxis()->GetLabelSize());
  for(int ff=0; ff<pdfvar::nfit; ff++)
    {
      hyields_a->GetXaxis()->SetBinLabel(ff+1, pdfvar::fitopt[ff].c_str());
      hBenryields_a->GetXaxis()->SetBinLabel(ff+1, pdfvar::fitopt[ff].c_str());
      hyields_b->GetXaxis()->SetBinLabel(ff+1, pdfvar::fitopt[ff].c_str());
      hBenryields_b->GetXaxis()->SetBinLabel(ff+1, pdfvar::fitopt[ff].c_str());
    }
  // fit + yield
  std::vector<TH1F*> vh = {h    , hBenr};
  std::vector<TH1F*> vhyields_a = {hyields_a    , hBenryields_a};
  std::vector<TH1F*> vhyields_b = {hyields_b    , hBenryields_b};
  std::vector<std::vector<TH1F*>> vhyield = {hyield    , hBenryield};
  std::vector<std::string> vname           = {"th"       , "thBenr"};
  std::vector<std::string> vtitle          = {"Inclusive", "B-enriched (l_{xy} > 0.1 mm)"};
  std::vector<Style_t> mstyle              = {20         , 24};

  for(int l=0; l<2; l++)
    {
      for(int ff=0; ff<pdfvar::nfit; ff++)
        {
          // ====>
          // std::vector<TF1*> funsft = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, "plots/fltm", false, false); // fix mean = false
          std::vector<TF1*> funs   = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, Form("plots/%s/pdfvar/idx", output.c_str()), true, 
                                               true, "_"+vname[l]+"_"+pdfvar::fitopt[ff], pdfvar::fitopt[ff]); // fix mean = true
          // <====
          float ysig_a = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
          float ysigerr_a = funs[0]->GetParError(5)*ysig_a/funs[0]->GetParameter(5);
          float ysig_b = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
          float ysigerr_b = funs[0]->GetParError(10)*ysig_b/funs[0]->GetParameter(10);
          // float ybkg1 = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH;
          // float ybkg2 = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH;
          // yield
          vhyield[l].at(ff)->SetBinContent(fitX::ibin_a, ysig_a);
          vhyield[l].at(ff)->SetBinError(fitX::ibin_a, ysigerr_a);
          vhyield[l].at(ff)->SetBinContent(fitX::ibin_b, ysig_b);
          vhyield[l].at(ff)->SetBinError(fitX::ibin_b, ysigerr_b);
          vhyields_a[l]->SetBinContent(ff+1, ysig_a);
          vhyields_a[l]->SetBinError(ff+1, ysigerr_a);
          vhyields_b[l]->SetBinContent(ff+1, ysig_b);
          vhyields_b[l]->SetBinError(ff+1, ysigerr_b);
          xjjroot::setthgrstyle(vhyield[l].at(ff), xjjroot::colorlist_middle[ff], mstyle[l], 1.2, xjjroot::colorlist_middle[ff], 2, 2);
        }
      xjjroot::setthgrstyle(vhyields_a[l], fitX::color_a, mstyle[l], 1.2, fitX::color_a, 2, 2);
      xjjroot::setthgrstyle(vhyields_b[l], fitX::color_b, mstyle[l], 1.2, fitX::color_b, 2, 2);
    }

  xjjroot::setgstyle();
  TCanvas* cy = new TCanvas("cy", "", 1200, 600);
  cy->Divide(2, 1);
  for(int l=0; l<2; l++)
    {
      cy->cd(l+1);
      float ymax = 0;
      for(int ff=0; ff<pdfvar::nfit; ff++) { if(vhyield[l].at(ff)->GetMaximum() > ymax) { ymax = vhyield[l].at(ff)->GetMaximum(); } }
      TH2F* hempty = new TH2F(Form("hempty%d",l), ";;Raw Yield", 5, 0, 5, 10, 0, ymax*1.4);
      xjjroot::sethempty(hempty, 0, 0.3);
      hempty->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
      hempty->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
      hempty->GetXaxis()->SetLabelSize(hempty->GetXaxis()->GetLabelSize()*1.5);
      hempty->Draw();
      for(auto& hh : vhyield[l]) { hh->Draw("same ple"); }
      xjjroot::drawtexgroup(0.22, 0.24, pdfvar::fitopt, 3, 0.13, 0.04, 11, 42, xjjroot::colorlist_middle);
      drawkinematic();
      xjjroot::drawtex(0.22, 0.84, vtitle[l].c_str(), 0.042, 12, 62);
      xjjroot::drawCMS();
    }
  cy->SaveAs(Form("plots/%s/pdfvar/cpdfvar1.pdf", output.c_str()));

  TCanvas* cyper = new TCanvas("cyper", "", 1200, 600);
  cyper->Divide(2, 1);
  for(int l=0; l<2; l++)
    {
      cyper->cd(l+1);
      vhyields_a[l]->Scale(1./vhyields_a[l]->GetBinContent(1));
      vhyields_b[l]->Scale(1./vhyields_b[l]->GetBinContent(1));
      float ymax = vhyields_a[l]->GetMaximum();
      TH2F* hempty = new TH2F(Form("hemptyper%d",l), ";;PDF variation / nominal", pdfvar::nfit, 0, pdfvar::nfit, 10, 0, ymax*1.8);
      xjjroot::sethempty(hempty, 0, 0.3);
      for(int ff=0; ff<pdfvar::nfit; ff++)
        { hempty->GetXaxis()->SetBinLabel(ff+1, pdfvar::fitopt[ff].c_str()); }
      hempty->GetXaxis()->SetLabelSize(hempty->GetXaxis()->GetLabelSize()*1.5);
      hempty->Draw();
      xjjroot::drawline(0, 1, pdfvar::nfit, 1, kGray+1, 10, 2);
      vhyields_a[l]->Draw("same ple");
      vhyields_b[l]->Draw("same ple");
      drawkinematic();
      xjjroot::drawtex(0.22, 0.84, vtitle[l].c_str(), 0.042, 12, 62);
      xjjroot::drawCMS();
      xjjroot::drawtex(0.24, 0.23+0.045, fitX::title_a.c_str(), 0.042, 11, 62, fitX::color_a);
      xjjroot::drawtex(0.24, 0.23, fitX::title_b.c_str(), 0.042, 11, 62, fitX::color_b);
    }
  cyper->SaveAs(Form("plots/%s/pdfvar/cpdfvar2.pdf", output.c_str()));

  // print
  float per_a = 0;
  for(int i=0; i<hyields_a->GetNbinsX(); i++)
    {
      if(TMath::Abs(hyields_a->GetBinContent(i+1)-1) > per_a) per_a = TMath::Abs(hyields_a->GetBinContent(i+1)-1);
    }
  float per_b = 0;
  for(int i=0; i<hyields_b->GetNbinsX(); i++)
    {
      if(TMath::Abs(hyields_b->GetBinContent(i+1)-1) > per_b) per_b = TMath::Abs(hyields_b->GetBinContent(i+1)-1);
    }
  float per_ab = TMath::Sqrt(per_a*per_a + per_b*per_b);
  std::cout<<"Yield Extraction & "<<Form("%.1f", per_a*1.e+2)<<"\\% & "<<Form("%.1f", per_b*1.e+2)<<"\\% & "<<Form("%.1f", per_ab*1.e+2)<<"\\% \\\\"<<std::endl;

}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX_pdfvar(argv[1], argv[2]); return 0; }
  return 1;
}

void drawkinematic()
{
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.042, 32, 42);
  xjjroot::drawtex(0.92, 0.77, fitX::pttag().c_str(), 0.042, 32, 42);
}
