#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TColor.h>
#include "xjjrootuti.h"

namespace ppref
{
  class CMS_2013_Rprompt
  {
  public:
    CMS_2013_Rprompt();
    TGraphAsymmErrors* grae_stat() { return fgrae_stat; }
    TGraphAsymmErrors* grae_syst() { return fgrae_syst; }
    void Draw();
  private:
    TGraphAsymmErrors* fgrae_stat;
    TGraphAsymmErrors* fgrae_syst;
  };
}

ppref::CMS_2013_Rprompt::CMS_2013_Rprompt()
{
  //=========Macro generated from canvas: c1/c1
  //=========  (Wed Jul 17 15:11:47 2019) by ROOT version 6.16/00
  // TCanvas *c1 = new TCanvas("c1", "c1",358,142,700,500);
  // gStyle->SetOptFit(1);
  // gStyle->SetOptStat(0);
  // c1->Range(-4.329115,0.04199192,56.43038,0.09815367);
  // c1->SetFillColor(0);
  // c1->SetBorderMode(0);
  // c1->SetBorderSize(0);
  // c1->SetTickx(1);
  // c1->SetTicky(1);
  // c1->SetLeftMargin(0.17);
  // c1->SetRightMargin(0.04);
  // c1->SetTopMargin(0.05);
  // c1->SetBottomMargin(0.15);
  // c1->SetFrameLineColor(0);
  // c1->SetFrameBorderMode(0);
  // c1->SetFrameLineColor(0);
  // c1->SetFrameBorderMode(0);
   
  Double_t Graph1D_y2_fx3001[5] = {
    11.75,
    14.25,
    16.5,
    24};
  Double_t Graph1D_y2_fy3001[5] = {
    0.0878,
    0.1009,
    0.1010,
    0.0987};
  Double_t Graph1D_y2_felx3001[5] = {
    1.75,
    0.75,
    1.5,
    6};
  // Double_t Graph1D_y2_fely3001[5] = {
  //   0.01251,
  //   0.008438009,
  //   0.007500667,
  //   0.005939697,
  //   0.01360147};
  Double_t Graph1D_y2_fehx3001[5] = {
    1.75,
    0.75,
    1.5,
    6};
  // Double_t Graph1D_y2_fehy3001[5] = {
  //   0.01251,
  //   0.008438009,
  //   0.007500667,
  //   0.005939697,
  //   0.01360147};
  Double_t Graph1D_y2_festaty3001[5] = {
    0.0122,
    0.0134,
    0.0110,
    0.0108};
  Double_t Graph1D_y2_fesysty3001[5] = {
    0.0123,
    0.0086,
    0.0085,
    0.0099};
  Double_t Graph1D_y2_fexsmall3001[5] = {
    0.5,
    0.5,
    0.5,
    0.5};
  fgrae_stat = new TGraphAsymmErrors(5,Graph1D_y2_fx3001,Graph1D_y2_fy3001,Graph1D_y2_felx3001,Graph1D_y2_fehx3001,Graph1D_y2_festaty3001,Graph1D_y2_festaty3001);
  fgrae_stat->SetName("grae_stat_prompt");
  fgrae_stat->SetTitle("doi:10.17182/hepdata.60421.v1/t1");
  fgrae_syst = new TGraphAsymmErrors(5,Graph1D_y2_fx3001,Graph1D_y2_fy3001,Graph1D_y2_fexsmall3001,Graph1D_y2_fexsmall3001,Graph1D_y2_fesysty3001,Graph1D_y2_fesysty3001);
  fgrae_syst->SetName("grae_syst_prompt");
  fgrae_syst->SetTitle("doi:10.17182/hepdata.60421.v1/t1");

  xjjroot::setthgrstyle(fgrae_syst, kGray+3, 21, 1.2, 0, 1, 1, kGray+3, 0.3, 1001);
  xjjroot::setthgrstyle(fgrae_stat, kGray+3, 21, 1.2, kGray+3, 1, 1, kGray+3, 0.3, 1001);
}

void ppref::CMS_2013_Rprompt::Draw()
{
  fgrae_syst->Draw("2 same");
  fgrae_stat->Draw("pe same");
}
