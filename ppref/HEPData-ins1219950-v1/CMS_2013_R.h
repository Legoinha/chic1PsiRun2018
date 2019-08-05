#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TColor.h>
#include "xjjrootuti.h"

namespace ppref
{
  void CMS_2013_R();
  TGraphAsymmErrors* grae_stat;
  TGraphAsymmErrors* grae_syst;
}

void ppref::CMS_2013_R()
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
    24,
    40};
  Double_t Graph1D_y2_fy3001[5] = {
    0.0727,
    0.0671,
    0.0687,
    0.0601,
    0.078};
  Double_t Graph1D_y2_felx3001[5] = {
    1.75,
    0.75,
    1.5,
    6,
    10};
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
    6,
    10};
  // Double_t Graph1D_y2_fehy3001[5] = {
  //   0.01251,
  //   0.008438009,
  //   0.007500667,
  //   0.005939697,
  //   0.01360147};
  Double_t Graph1D_y2_festaty3001[5] = {
    0.0079,
    0.0072,
    0.0055,
    0.0042,
    0.013};
  Double_t Graph1D_y2_fesysty3001[5] = {
    0.0097,
    0.0044,
    0.0051,
    0.0042,
    0.004};
  Double_t Graph1D_y2_fexsmall3001[5] = {
    0.5,
    0.5,
    0.5,
    0.5,
    0.5};
  grae_stat = new TGraphAsymmErrors(5,Graph1D_y2_fx3001,Graph1D_y2_fy3001,Graph1D_y2_felx3001,Graph1D_y2_fehx3001,Graph1D_y2_festaty3001,Graph1D_y2_festaty3001);
  grae_stat->SetName("grae_stat");
  grae_stat->SetTitle("doi:10.17182/hepdata.60421.v1/t1");
  grae_syst = new TGraphAsymmErrors(5,Graph1D_y2_fx3001,Graph1D_y2_fy3001,Graph1D_y2_fexsmall3001,Graph1D_y2_fexsmall3001,Graph1D_y2_fesysty3001,Graph1D_y2_fesysty3001);
  grae_syst->SetName("grae_syst");
  grae_syst->SetTitle("doi:10.17182/hepdata.60421.v1/t1");

  xjjroot::setthgrstyle(grae_syst, kOrange+8, 20, 1.2, 0, 1, 1, kOrange+8, 0.3, 1001);
  xjjroot::setthgrstyle(grae_stat, kOrange+8, 20, 1.2, kOrange+8, 1, 1, kOrange+8, 0.3, 1001);
  
  grae_syst->Draw("2 same");
  grae_stat->Draw("pe same");
  // TH1F *Graph_Graph1D_y23001 = new TH1F("Graph_Graph1D_y23001","doi:10.17182/hepdata.60421.v1/t1",100,6,54);
  // Graph_Graph1D_y23001->SetMinimum(0.05041619);
  // Graph_Graph1D_y23001->SetMaximum(0.09534559);
  // Graph_Graph1D_y23001->SetDirectory(0);
  // Graph_Graph1D_y23001->SetStats(0);

  // Int_t ci;      // for color index setting
  // TColor *color; // for color definition with alpha
  // ci = TColor::GetColor("#000099");
  // Graph_Graph1D_y23001->SetLineColor(ci);
  // Graph_Graph1D_y23001->GetXaxis()->SetTitle("PT [GEV]");
  // Graph_Graph1D_y23001->GetXaxis()->SetLabelFont(42);
  // Graph_Graph1D_y23001->GetXaxis()->SetLabelSize(0.035);
  // Graph_Graph1D_y23001->GetXaxis()->SetTitleSize(0.035);
  // Graph_Graph1D_y23001->GetXaxis()->SetTitleOffset(1);
  // Graph_Graph1D_y23001->GetXaxis()->SetTitleFont(42);
  // Graph_Graph1D_y23001->GetYaxis()->SetTitle("R");
  // Graph_Graph1D_y23001->GetYaxis()->SetLabelFont(42);
  // Graph_Graph1D_y23001->GetYaxis()->SetLabelSize(0.035);
  // Graph_Graph1D_y23001->GetYaxis()->SetTitleSize(0.035);
  // Graph_Graph1D_y23001->GetYaxis()->SetTitleFont(42);
  // Graph_Graph1D_y23001->GetZaxis()->SetLabelFont(42);
  // Graph_Graph1D_y23001->GetZaxis()->SetLabelSize(0.035);
  // Graph_Graph1D_y23001->GetZaxis()->SetTitleSize(0.035);
  // Graph_Graph1D_y23001->GetZaxis()->SetTitleOffset(1);
  // Graph_Graph1D_y23001->GetZaxis()->SetTitleFont(42);
  // grae->SetHistogram(Graph_Graph1D_y23001);
   
  // grae->Draw("same pe");

  // TPaveText *pt = new TPaveText(0.2220917,0.9338535,0.7779083,0.995,"blNDC");
  // pt->SetName("title");
  // pt->SetBorderSize(0);
  // pt->SetFillColor(0);
  // pt->SetFillStyle(0);
  // pt->SetTextFont(42);
  // TText *pt_LaTex = pt->AddText("doi:10.17182/hepdata.60421.v1/t1");
  // pt->Draw();
  // c1->Modified();
  // c1->cd();
  // c1->SetSelected(c1);
}
