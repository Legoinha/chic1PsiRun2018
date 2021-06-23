#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <string>
#include <vector>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitX.h"

float correrr(TH1* h, float ref);
int nbins = 9; float minbin = 0.04, divbin = 0.01, maxbin = minbin+nbins*divbin;
float valref = 0.06;
void fitX_drawvary(std::string inputname, std::string output)
{
  TH1F* hyieldpromptCorr_a = new TH1F("hyieldpromptCorr_a", ";BDT >;", nbins, minbin, maxbin);
  xjjroot::sethempty(hyieldpromptCorr_a, 0, 0);
  TH1F* hyieldpromptCorr_b = new TH1F("hyieldpromptCorr_b", ";BDT >;", nbins, minbin, maxbin);
  xjjroot::sethempty(hyieldpromptCorr_b, 0, 0);
  TH1F* hratio_ab = new TH1F("hratio_ab", ";BDT >;", nbins, minbin, maxbin);
  xjjroot::sethempty(hratio_ab, 0, 0);
  for(int i=0; i<nbins; i++)
    {
      float mvaval = minbin + divbin*i;
      TFile* inf = TFile::Open(xjjc::str_replaceall(inputname, "BDTQvalue", "BDTQvalue"+xjjc::number_to_string(mvaval)).c_str());
      TH1F* hyieldpromptCorr = (TH1F*)inf->Get("hyieldpromptCorr");
      TH1F* hratio = (TH1F*)inf->Get("hratio");
      hyieldpromptCorr_a->SetBinContent(i+1, hyieldpromptCorr->GetBinContent(fitX::ibin_a));
      hyieldpromptCorr_a->SetBinError(i+1, hyieldpromptCorr->GetBinError(fitX::ibin_a));
      hyieldpromptCorr_b->SetBinContent(i+1, hyieldpromptCorr->GetBinContent(fitX::ibin_b));
      hyieldpromptCorr_b->SetBinError(i+1, hyieldpromptCorr->GetBinError(fitX::ibin_b));
      hratio_ab->SetBinContent(i+1, hratio->GetBinContent(1));
      hratio_ab->SetBinError(i+1, hratio->GetBinError(1));
      delete hyieldpromptCorr;
      delete hratio;
    }

  int ibinref = hratio_ab->FindBin(valref);
  float errref_hyieldpromptCorr_a = correrr(hyieldpromptCorr_a, valref), ref_hyieldpromptCorr_a = hyieldpromptCorr_a->GetBinContent(ibinref);
  TGraphErrors* gyieldpromptCorr_a = xjjroot::shifthistcenter(hyieldpromptCorr_a, "gyieldpromptCorr_a");
  xjjroot::setthgrstyle(gyieldpromptCorr_a, fitX::color_a, 47, 1.4, fitX::color_a, 1, 2);
  float errref_hyieldpromptCorr_b = correrr(hyieldpromptCorr_b, valref), ref_hyieldpromptCorr_b = hyieldpromptCorr_b->GetBinContent(ibinref);
  TGraphErrors* gyieldpromptCorr_b = xjjroot::shifthistcenter(hyieldpromptCorr_b, "gyieldpromptCorr_b");
  xjjroot::setthgrstyle(gyieldpromptCorr_b, fitX::color_b, 47, 1.4, fitX::color_b, 1, 2);
  float errref_hratio_ab = correrr(hratio_ab, valref), ref_hratio_ab = hratio_ab->GetBinContent(ibinref);
  TGraphErrors* gratio_ab = xjjroot::shifthistcenter(hratio_ab, "gratio_ab");
  xjjroot::setthgrstyle(gratio_ab, kBlack, 47, 1.4, kBlack, 1, 2);

  TH2F* hempty = new TH2F("hempty", ";BDT >;N_{signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", 
                          nbins+1, minbin-divbin, maxbin, 
                          10, 0, std::max(hyieldpromptCorr_a->GetMaximum(), hyieldpromptCorr_b->GetMaximum())*1.5);
  hempty->GetXaxis()->SetNdivisions(505);
  xjjroot::sethempty(hempty, 0, 0.3);
  TH2F* hempty_r = new TH2F("hempty_r", ";BDT >;R", 
                            nbins+1, minbin-divbin, maxbin, 
                            10, 0, 2);
  hempty_r->GetXaxis()->SetNdivisions(505);
  xjjroot::sethempty(hempty_r, 0, 0.3);

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hempty->Draw();
  xjjroot::drawbox(hempty->GetXaxis()->GetXmin(), ref_hyieldpromptCorr_a-errref_hyieldpromptCorr_a, hempty->GetXaxis()->GetXmax(), ref_hyieldpromptCorr_a+errref_hyieldpromptCorr_a, 
                   fitX::color_a, 0.1, 1001);
  xjjroot::drawline(hempty->GetXaxis()->GetXmin(), ref_hyieldpromptCorr_a, hempty->GetXaxis()->GetXmax(), ref_hyieldpromptCorr_a,
                    fitX::color_a, 2, 3, 0.5);
  xjjroot::drawbox(hempty->GetXaxis()->GetXmin(), ref_hyieldpromptCorr_b-errref_hyieldpromptCorr_b, hempty->GetXaxis()->GetXmax(), ref_hyieldpromptCorr_b+errref_hyieldpromptCorr_b, 
                   fitX::color_b, 0.1, 1001);
  xjjroot::drawline(hempty->GetXaxis()->GetXmin(), ref_hyieldpromptCorr_b, hempty->GetXaxis()->GetXmax(), ref_hyieldpromptCorr_b,
                    fitX::color_b, 2, 3, 0.5);
  gyieldpromptCorr_a->Draw("pe same");
  gyieldpromptCorr_b->Draw("pe same");
  xjjroot::drawCMS("Internal");
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.04, fitX::title_b.c_str(), 0.038, 12, 62, fitX::color_b);
  fitX::drawkinematics();
  c->cd(2);
  hempty_r->Draw();
  xjjroot::drawbox(hempty_r->GetXaxis()->GetXmin(), ref_hratio_ab-errref_hratio_ab, hempty_r->GetXaxis()->GetXmax(), ref_hratio_ab+errref_hratio_ab, 
                   kGray+2, 0.1, 1001);
  xjjroot::drawline(hempty_r->GetXaxis()->GetXmin(), ref_hratio_ab, hempty_r->GetXaxis()->GetXmax(), ref_hratio_ab,
                    kGray+2, 2, 3, 0.5);
  gratio_ab->Draw("pe same");
  xjjroot::drawCMS("Internal");
  fitX::drawkinematics();

  std::string outputname = "plots/"+output+"/cmvavary_BDT.pdf";
  xjjroot::mkdir(outputname);
  c->SaveAs(outputname.c_str());

  float dmax_a = 0, dmax_b = 0, dmax_r = 0;
  for(int i=0; i<nbins; i++)
    {
      if(fabs(hyieldpromptCorr_a->GetBinContent(i+1)-ref_hyieldpromptCorr_a) > dmax_a) { dmax_a = fabs(hyieldpromptCorr_a->GetBinContent(i+1)-ref_hyieldpromptCorr_a); }
      if(fabs(hyieldpromptCorr_b->GetBinContent(i+1)-ref_hyieldpromptCorr_b) > dmax_b) { dmax_b = fabs(hyieldpromptCorr_b->GetBinContent(i+1)-ref_hyieldpromptCorr_b); }
      if(fabs(hratio_ab->GetBinContent(i+1)-ref_hratio_ab) > dmax_r) { dmax_r = fabs(hratio_ab->GetBinContent(i+1)-ref_hratio_ab); }
    }
  dmax_a/=ref_hyieldpromptCorr_a;
  dmax_b/=ref_hyieldpromptCorr_b;
  dmax_r/=ref_hratio_ab;
  // std::cout<<dmax_a<<" "<<dmax_b<<" "<<dmax_r<<std::endl;
  std::cout<<"      Efficiency       & "<<Form("%.1f", dmax_a*100.)<<"\\%  & "<<Form("%.1f", dmax_b*100.)<<"\\%  & "<<Form("%.1f", dmax_r*100.)<<"\\% \\\\"<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX_drawvary(argv[1], argv[2]); return 0; }
  return 1;
}

float correrr(TH1* h, float ref)
{
  int ibin = h->FindBin(ref);
  float errref = h->GetBinError(ibin);
  for(int i=0; i<h->GetXaxis()->GetNbins(); i++)
    h->SetBinError(i+1, sqrt(fabs(h->GetBinError(i+1)*h->GetBinError(i+1)-errref*errref)));
  return errref;
}
