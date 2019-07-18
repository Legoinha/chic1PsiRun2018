#include "xjjrootuti.h"

#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>

#include <string>

#include "ntuple.h"

namespace fitX
{
  void setmasshist(TH1* h, float xoffset=0, float yoffset=0, Color_t pcolor=kBlack);
  std::vector<TF1*> fit(TH1F* h, TH1F* h_ss, TH1F* hmc_a, TH1F* hmc_b, std::string name, bool fixmean, bool saveplot);

  const int NBIN = 38, NBIN_L = 50, NBIN_H = 50;
  // const int NBIN = 57, NBIN_L = 50, NBIN_H = 50;
  const float BIN_MIN = 3.62, BIN_MAX = 4.0, BIN_MIN_L = 3.64, BIN_MAX_L = 3.74, BIN_MIN_H = 3.82, BIN_MAX_H = 3.92;
  float BIN_WIDTH = (BIN_MAX-BIN_MIN)/NBIN*1.0, BIN_WIDTH_L = (BIN_MAX_L-BIN_MIN_L)/NBIN_L*1.0, BIN_WIDTH_H = (BIN_MAX_H-BIN_MIN_H)/NBIN_H*1.0;
  void drawpull(TH1* hmc, TF1* f, Color_t color);

  Color_t color_data = kRed-3, color_a = kAzure+4, color_b = kGreen-1, color_ss = kGray+1, color_bkg = color_data;
}


// --->
std::vector<TF1*> fitX::fit(TH1F* h, TH1F* h_ss, TH1F* hmc_a, TH1F* hmc_b, std::string name, bool fixmean, bool saveplot)
{
  std::vector<TF1*> funs;
  xjjroot::setgstyle(2);
  gStyle->SetPadLeftMargin(gStyle->GetPadLeftMargin()*0.7);

  fitX::setmasshist(h, 0, -0.2);
  h->SetMinimum(0);
  if(h_ss) fitX::setmasshist(h_ss, 0, -0.2, color_ss);
  fitX::setmasshist(hmc_a, 0, 0.05);
  hmc_a->SetMaximum(hmc_a->GetMaximum()*1.2);
  hmc_a->SetMinimum(0 - hmc_a->GetMaximum()*0.1);
  hmc_a->SetNdivisions(-505);
  fitX::setmasshist(hmc_b, 0, 0.05);
  hmc_b->SetMaximum(hmc_b->GetMaximum()*1.2);
  hmc_b->SetMinimum(0 - hmc_b->GetMaximum()*0.1);
  hmc_b->SetNdivisions(-505);

  TCanvas* c = new TCanvas("c", "", 800, 600);
  TCanvas* cmc = new TCanvas("cmc", "", 1200, 600);
  cmc->Divide(2, 1);
  
  TString str_bkg = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x";
  TString str_sig1 = "[5]*([9]*TMath::Gaus(x,[6],[7])/(TMath::Sqrt(2*3.14159)*[7]) + (1-[9])*TMath::Gaus(x,[6],[8])/(TMath::Sqrt(2*3.14159)*[8]))";
  TString str_sig2 = "[10]*([14]*TMath::Gaus(x,[11],[12])/(TMath::Sqrt(2*3.14159)*[7]) + (1-[14])*TMath::Gaus(x,[11],[13])/(TMath::Sqrt(2*3.14159)*[13]))";
  TF1 *f = new TF1("f", str_bkg+"+"+str_sig1+"+"+str_sig2, BIN_MIN, BIN_MAX);
  const double parinit[15] = {0., 0.     , 0.     , 0.  , 0. , // 0-4
                              1., 3.686  , 0.00357, 0.03, 0.5, // 5-9
                              1., 3.87169, 0.001  , 0.01, 0.5, // 10-14
  };
  f->SetParameters(parinit);
  f->SetParLimits(9, 0, 1);
  f->SetParLimits(14, 0, 1);
  f->FixParameter(0, 0);
  f->FixParameter(1, 0);
  f->FixParameter(2, 0);
  f->FixParameter(3, 0);
  f->FixParameter(4, 0);
  //
  cmc->cd(1);
  hmc_a->Draw("pe");
  xjjroot::settfstyle(f, color_a, 1, 2);
  f->FixParameter(10, 0);
  f->FixParameter(6, 3.68610);
  hmc_a->Fit("f", "Nq");
  hmc_a->Fit("f", "NLLq");
  f->ReleaseParameter(6);
  hmc_a->Fit("f", "NLL");
  f->FixParameter(6, f->GetParameter(6));
  f->FixParameter(7, f->GetParameter(7));
  f->FixParameter(8, f->GetParameter(8));
  f->FixParameter(9, f->GetParameter(9));
  f->ReleaseParameter(10);
  fitX::drawpull(hmc_a, f, color_a);
  xjjroot::copyobject(f, "fmc_a")->Draw("same");
  hmc_a->Draw("pesame");
  xjjroot::drawtex(0.18, 0.85, "Gen-matched #psi(2S)", 0.04, 12, 62);
  xjjroot::drawtex(0.65, 0.86, Form("#bar{m} = %.4f GeV", f->GetParameter(6)), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05, Form("#sigma_{1} = %.4f GeV", std::min(f->GetParameter(7),f->GetParameter(8))), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05*2, Form("#sigma_{2} = %.4f GeV", std::max(f->GetParameter(7),f->GetParameter(8))), 0.038);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  // 
  cmc->cd(2);
  hmc_b->Draw("pe");
  xjjroot::settfstyle(f, color_b, 1, 2);
  f->SetLineColor(color_b);
  f->FixParameter(5, 0);
  f->FixParameter(11, 3.87169);
  hmc_b->Fit("f", "Nq");
  hmc_b->Fit("f", "NLLq");
  f->ReleaseParameter(11);
  hmc_b->Fit("f", "NLL");
  f->FixParameter(11, f->GetParameter(11));
  f->FixParameter(12, f->GetParameter(12));
  f->FixParameter(13, f->GetParameter(13));
  f->FixParameter(14, f->GetParameter(14));
  f->ReleaseParameter(5);
  fitX::drawpull(hmc_b, f, color_b);
  xjjroot::copyobject(f, "fmc_b")->Draw("same");
  hmc_b->Draw("pesame");
  xjjroot::drawtex(0.18, 0.85, "Gen-matched X(3872)", 0.042, 12, 62);
  xjjroot::drawtex(0.65, 0.86, Form("#bar{m} = %.4f GeV", f->GetParameter(11)), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05, Form("#sigma_{1} = %.4f GeV", std::min(f->GetParameter(12),f->GetParameter(13))), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05*2, Form("#sigma_{2} = %.4f GeV", std::max(f->GetParameter(12),f->GetParameter(13))), 0.038);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  cmc->RedrawAxis();
  if(saveplot) cmc->SaveAs(Form("plots/chmassmc_%s.pdf", name.c_str()));

  //
  xjjroot::settfstyle(f, color_data, 1, 3);
  f->ReleaseParameter(0);
  f->ReleaseParameter(1);
  f->ReleaseParameter(2);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  // >>>
  if(!fixmean)
    {
      f->ReleaseParameter(6);
      f->SetParLimits(6, MASS_PSI2S-mytmva::sigwindowL, MASS_PSI2S+mytmva::sigwindowL);
      f->ReleaseParameter(11);
      f->SetParLimits(11, MASS_X-mytmva::sigwindowH, MASS_X+mytmva::sigwindowH);
    }
  // <<<
  f->SetParLimits(5, 0, 1.e+5);
  f->SetParLimits(10, 0, 1.e+5);

  c->cd();
  f->SetLineColor(color_data);
  h->SetMaximum(h->GetMaximum()*1.25);
  h->Draw("pe");
  h->Fit("f","Nq");
  h->Fit("f","NLLq");
  // TFitResultPtr r = histo->Fit(func, "S");
  h->Fit("f","NLL");

  fitX::drawpull(h, f, color_data);
  if(h_ss) h_ss->Draw("pe same");

  TF1 *fsig1 = xjjroot::copyobject(f, "fsig1");
  fsig1->SetRange(BIN_MIN_L, BIN_MAX_L);
  fsig1->FixParameter(10, 0);
  fsig1->FixParameter(0, 0);
  fsig1->FixParameter(1, 0);
  fsig1->FixParameter(2, 0);
  fsig1->FixParameter(3, 0);
  fsig1->FixParameter(4, 0);
  fsig1->SetParError(5, f->GetParError(5));
  xjjroot::settfstyle(fsig1, color_a, 2, 3, color_a, 0.3, 1001);
  fsig1->Draw("same");

  TF1 *fsig2 = xjjroot::copyobject(f, "fsig2");
  fsig2->SetRange(BIN_MIN_H, BIN_MAX_H);
  fsig2->FixParameter(5, 0);
  fsig2->FixParameter(0, 0);
  fsig2->FixParameter(1, 0);
  fsig2->FixParameter(2, 0);
  fsig2->FixParameter(3, 0);
  fsig2->FixParameter(4, 0);
  fsig2->SetParError(10, f->GetParError(10));
  xjjroot::settfstyle(fsig2, color_b, 2, 3, color_b, 0.3, 1001);
  fsig2->Draw("same");

  TF1 *fbkg = xjjroot::copyobject(f, "fbkg");
  fbkg->SetRange(BIN_MIN, BIN_MAX);
  fbkg->FixParameter(5, 0);
  fbkg->FixParameter(10, 0);
  xjjroot::settfstyle(fbkg, color_bkg, 2, 3);
  fbkg->Draw("same");

  f->Draw("same");
  xjjroot::copyobject(h, "hclone")->Draw("pesame");  
  c->RedrawAxis();
  float ysig1 = fsig1->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;
  float ysig2 = fsig2->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;

  xjjroot::drawtex(0.92, 0.84, "p_{T} > 15 GeV/c", 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.79, "|y| < 1.5", 0.042, 32, 62);
  xjjroot::drawtex(0.17, 0.84, Form("#chi^{2} Prob = %.1f%s", TMath::Prob(f->GetChisquare(), f->GetNDF())*100., "%"), 0.042, 12, 62);

  xjjroot::drawtex(0.32, 0.36, Form("#bar{m}_{#psi(2S)} = %.4f GeV", f->GetParameter(6)), 0.04, 12, 62, color_a);
  xjjroot::drawtex(0.32, 0.31, Form("N_{#psi(2S)} = %.0f #pm %.0f", ysig1, f->GetParError(5)*ysig1/f->GetParameter(5)), 0.04, 12, 62, color_a);
  xjjroot::drawtex(0.68, 0.36, Form("#bar{m}_{X(3872)} = %.4f GeV", f->GetParameter(11)), 0.04, 12, 62, color_b);
  xjjroot::drawtex(0.68, 0.31, Form("N_{X(3872)} = %.0f #pm %.0f", ysig2, f->GetParError(10)*ysig2/f->GetParameter(10)), 0.04, 12, 62, color_b);
  xjjroot::drawCMS();

  if(saveplot) c->SaveAs(Form("plots/chmass_%s.pdf", name.c_str()));

  funs.push_back(f);
  funs.push_back(fsig1);
  funs.push_back(fsig2);
  funs.push_back(fbkg);
  return funs;
}

void fitX::setmasshist(TH1* h, float xoffset/*=0*/, float yoffset/*=0*/, Color_t pcolor/*=kBlack*/)
{
  xjjroot::sethempty(h, xoffset, yoffset);
  h->Sumw2();
  h->SetLineColor(pcolor);
  h->SetMarkerColor(pcolor);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.0);
}

void fitX::drawpull(TH1* hmc, TF1* f, Color_t color)
{
  int nbin = hmc->GetXaxis()->GetNbins();
  float binmin = hmc->GetBinCenter(1) - hmc->GetBinWidth(1)/2.;
  float binmax = hmc->GetBinCenter(nbin) + hmc->GetBinWidth(nbin)/2.;
  for(int bb=0; bb<nbin; bb++) 
    { 
      float realval = hmc->GetBinError(bb+1)==0?0:(hmc->GetBinContent(bb+1)-f->Eval(hmc->GetBinCenter(bb+1)))/hmc->GetBinError(bb+1);
      float fillval = ((realval+4)/(4*2))*hmc->GetMaximum();
      xjjroot::drawbox(hmc->GetBinCenter(bb+1)-hmc->GetBinWidth(bb+1)/2., hmc->GetMaximum()/2., hmc->GetBinCenter(bb+1)+hmc->GetBinWidth(bb+1)/2., fillval, color, 0.1, 1001);
    }
  xjjroot::drawline(binmin, hmc->GetMaximum()/2., binmax, hmc->GetMaximum()/2., kGray, 2, 2);
  xjjroot::drawaxis(binmax, 0, binmax, hmc->GetMaximum(), -4, 4, color, 1, 2, "+L");
  xjjroot::drawtex(0.93, 0.55, "Pull", 0.04, 33, 62, color);
}
