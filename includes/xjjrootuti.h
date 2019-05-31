#ifndef _XJJROOT_UTI_H_
#define _XJJROOT_UTI_H_

/* xjjrootuti.h */

#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TH1.h>
#include <TGaxis.h>

namespace xjjroot
{
  const float margin_pad_left = 0.18;
  const float margin_pad_right = 0.043;
  const float margin_pad_bottom = 0.145;
  const float margin_pad_top = 0.1;

  const Color_t colorlist_light[]  = {kGreen-8, kOrange-4, kRed-9, kAzure-9, kMagenta-8, kCyan-10, kYellow-10, kBlue-8, kPink+1, kViolet-9};
  const Color_t colorlist_middle[] = {kGreen+2, kOrange-3, kRed-3, kAzure-3, kMagenta-5, kCyan-6,  kYellow+2,  kBlue-5, kPink+2, kViolet+7};
  const Color_t colorlist_dark[]   = {kGreen+3, kOrange+5, kRed+2, kAzure-6, kMagenta-1, kCyan+3,  kYellow+3,  kBlue-1, kPink+3, kViolet+4};
  const int ncolor = sizeof(colorlist_light)/sizeof(colorlist_light[0]);

  void setgstyle(Int_t padtick=0, Width_t lwidth=2);
  template <class T> void sethempty(T* hempty, Float_t xoffset=0, Float_t yoffset=0, Float_t xsize=0.05, Float_t ysize=0.05);
  template <class T> void setthgr(T* hempty, Float_t xoffset=0, Float_t yoffset=0);
  template <class T> void setthgrstyle(T* h, Color_t mcolor=-1, Style_t mstyle=-1, Size_t msize=-1, Color_t lcolor=-1, Style_t lstyle=-1, Width_t lwidth=-1, Color_t fcolor=-1, Float_t falpha=-1, Style_t fstyle=-1);
  template <class T> void setlinestyle(T* h, Color_t lcolor=-1, Style_t lstyle=-1, Width_t lwidth=-1);
  template <class T> void settfstyle(T* h, Color_t lcolor=-1, Style_t lstyle=-1, Width_t lwidth=-1, Color_t fcolor=-1, Float_t falpha=-1, Style_t fstyle=-1);
  template <class T> void setmarkerstyle(T* h, Color_t mcolor=-1, Style_t mstyle=-1, Size_t msize=-1);
  void drawCMSleft(TString content="#scale[1.25]{#bf{CMS}} #it{Preliminary}", Float_t xpos=0, Float_t ypos=0);
  void drawCMSright(TString content="PbPb #sqrt{s_{NN}} = 5.02 TeV", Float_t xpos=0, Float_t ypos=0);
  void drawCMS(TString collsyst="PbPb", TString snn="5.02");
  void settex(TLatex* tex, Float_t tsize=0.04, Short_t align=12, Style_t font=42, Color_t color=kBlack);
  void drawtex(Double_t x, Double_t y, const char *text, Float_t tsize=0.04, Short_t align=12, Style_t font=42, Color_t color=kBlack);
  void setleg(TLegend* leg, Float_t tsize=0.04);
  void setlegndraw(TLegend* leg, Float_t tsize=0.04);
  void setline(TLine* l, Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=2);
  void drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=2);
  void drawbox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t fcolor=kGray, Float_t falpha=0.4, Style_t fstyle=1001, Color_t lcolor=0, Style_t lstyle=1, Width_t lwidth=0);

  void drawaxis(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, 
                Double_t wmin, Double_t wmax, 
                Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=1,
                Option_t *chopt="", Int_t ndiv=510, Double_t gridlength=0);

  void dividebinwid(TH1* h);

  void setbranchaddress(TTree* nt, const char* bname, void* addr);
  template <class T> T* copyobject(const T* obj, TString objname);
}

/* ---------- */

void xjjroot::setgstyle(Int_t padtick/*=0*/, Width_t lwidth/*=2*/)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(xjjroot::margin_pad_right);
  gStyle->SetPadLeftMargin(xjjroot::margin_pad_left);
  gStyle->SetPadTopMargin(xjjroot::margin_pad_top);
  gStyle->SetPadBottomMargin(xjjroot::margin_pad_bottom);
  gStyle->SetTitleX(.0f);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetLineWidth(1);
  gStyle->SetFrameLineWidth(1);
  if(padtick==1 || padtick==3)
    {
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
    }
  if(padtick==2 || padtick==3) { gStyle->SetLineWidth(lwidth); gStyle->SetFrameLineWidth(lwidth);}
}

template <class T>
void xjjroot::sethempty(T* hempty, Float_t xoffset/*=0*/, Float_t yoffset/*=0*/, Float_t xsize/*=0.05*/, Float_t ysize/*=0.05*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(xsize);
  hempty->GetYaxis()->SetTitleSize(ysize);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
}

template <class T>
void xjjroot::setthgr(T* hempty, Float_t xoffset/*=0*/, Float_t yoffset/*=0*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
  hempty->SetMarkerStyle(20);
  hempty->SetMarkerSize(1.1);
  hempty->SetMarkerColor(kBlack);
  hempty->SetLineColor(kBlack);
}

template <class T>
void xjjroot::setthgrstyle(T* h, Color_t mcolor/*=-1*/, Style_t mstyle/*=-1*/, Size_t msize/*=-1*/, Color_t lcolor/*=-1*/, Style_t lstyle/*=-1*/, Width_t lwidth/*=-1*/, Color_t fcolor/*=-1*/, Float_t falpha/*=-1*/, Style_t fstyle/*=-1*/)
{
  if(mcolor>=0) h->SetMarkerColor(mcolor);
  if(mstyle>=0) h->SetMarkerStyle(mstyle);
  if(msize>=0)  h->SetMarkerSize(msize);
  if(lcolor>=0) h->SetLineColor(lcolor);
  if(lstyle>=0) h->SetLineStyle(lstyle);
  if(lwidth>=0) h->SetLineWidth(lwidth);
  if(fcolor>=0) h->SetFillColor(fcolor);
  if(falpha>=0) h->SetFillColorAlpha(fcolor, falpha);
  if(fstyle>=0) h->SetFillStyle(fstyle);
}

template <class T>
void xjjroot::setlinestyle(T* h, Color_t lcolor/*=-1*/, Style_t lstyle/*=-1*/, Width_t lwidth/*=-1*/)
{
  if(lcolor>=0) h->SetLineColor(lcolor);
  if(lstyle>=0) h->SetLineStyle(lstyle);
  if(lwidth>=0) h->SetLineWidth(lwidth);
}

template <class T>
void xjjroot::settfstyle(T* h, Color_t lcolor/*=-1*/, Style_t lstyle/*=-1*/, Width_t lwidth/*=-1*/, Color_t fcolor/*=-1*/, Float_t falpha/*=-1*/, Style_t fstyle/*=-1*/)
{
  h->SetNpx(1000);
  if(lcolor>=0) h->SetLineColor(lcolor);
  if(lstyle>=0) h->SetLineStyle(lstyle);
  if(lwidth>=0) h->SetLineWidth(lwidth);
  if(fcolor>=0) h->SetFillColor(fcolor);
  if(falpha>=0) h->SetFillColorAlpha(fcolor, falpha);
  if(fstyle>=0) h->SetFillStyle(fstyle);
}

template <class T>
void xjjroot::setmarkerstyle(T* h, Color_t mcolor/*=-1*/, Style_t mstyle/*=-1*/, Size_t msize/*=-1*/)
{
  if(mcolor>=0) h->SetMarkerColor(mcolor);
  if(mstyle>=0) h->SetMarkerStyle(mstyle);
  if(msize>=0)  h->SetMarkerSize(msize);
}

void xjjroot::drawCMS(TString collsyst/*="PbPb"*/, TString snn/*="5.02"*/)
{
  drawCMSleft();
  drawCMSright(collsyst+" #sqrt{s_{NN}} = "+snn+" TeV");
}

void xjjroot::drawCMSleft(TString content/*=""*/, Float_t xpos/*=0*/, Float_t ypos/*=0*/)
{
  TLatex* texCms = new TLatex(gStyle->GetPadLeftMargin()+xpos,(1-gStyle->GetPadTopMargin())*1.02+ypos, content.Data());
  texCms->SetNDC();
  texCms->SetTextAlign(11);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();
}

void xjjroot::drawCMSright(TString content/*="PbPb #sqrt{s_{NN}} = 5.02 TeV"*/, Float_t xpos/*=0*/, Float_t ypos/*=0*/)
{
  TLatex* texCol = new TLatex((1-gStyle->GetPadRightMargin())+xpos,(1-gStyle->GetPadTopMargin())*1.02+ypos, content.Data());
  texCol->SetNDC();
  texCol->SetTextAlign(31);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();
}

void xjjroot::settex(TLatex* tex, Float_t tsize/*=0.04*/, Short_t align/*=12*/, Style_t font/*=42*/, Color_t color/*=kBlack*/)
{
  tex->SetNDC();
  tex->SetTextFont(font);
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
  tex->SetTextColor(color);
}

void xjjroot::drawtex(Double_t x, Double_t y, const char* text, Float_t tsize/*=0.04*/, Short_t align/*=12*/, Style_t font/*=42*/, Color_t color/*=kBlack*/)
{
  TLatex* tex = new TLatex(x, y, text);
  xjjroot::settex(tex, tsize, align, font, color);
  tex->Draw();
}

void xjjroot::setleg(TLegend* leg, Float_t tsize/*=0.04*/)
{
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(tsize);
}

void xjjroot::setlegndraw(TLegend* leg, Float_t tsize/*=0.04*/)
{
  xjjroot::setleg(leg, tsize);
  leg->Draw();
}

void xjjroot::setline(TLine* l, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/)
{
  l->SetLineColor(lcolor);
  l->SetLineStyle(lstyle);
  l->SetLineWidth(lwidth);
}

void xjjroot::drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/)
{
  TLine* l = new TLine(x1, y1, x2, y2);
  xjjroot::setline(l, lcolor, lstyle, lwidth);
  l->Draw();
}

void xjjroot::drawbox(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t fcolor/*=kGray*/, Float_t falpha/*=0.4*/, Style_t fstyle/*=1001*/, Color_t lcolor/*=0*/, Style_t lstyle/*=1*/, Width_t lwidth/*=0*/)
{
  TBox* b = new TBox(x1, y1, x2, y2);
  b->SetFillColor(fcolor);
  b->SetFillColorAlpha(fcolor, falpha);
  b->SetFillStyle(fstyle);
  b->SetLineColor(lcolor);
  b->SetLineStyle(lstyle);
  b->SetLineWidth(lwidth);
  b->Draw();
}

void xjjroot::drawaxis(Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, 
                       Double_t wmin, Double_t wmax, 
                       Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=1*/,
                       Option_t *chopt/*=""*/, Int_t ndiv/*=510*/, Double_t gridlength/*=0*/)
{
  TGaxis* g = new TGaxis(xmin, ymin, xmax, ymax,
                         wmin, wmax, ndiv, chopt, gridlength);
  g->SetLabelColor(lcolor);
  g->SetTitleColor(lcolor);
  g->SetLineColor(lcolor);
  g->SetLineStyle(lstyle);
  g->SetLineWidth(lwidth);
  g->Draw();
}

/* ----- */

void xjjroot::setbranchaddress(TTree* nt, const char* bname, void* addr)
{
  nt->SetBranchStatus(bname, 1);
  nt->SetBranchAddress(bname, addr);
}

template <class T>
T* xjjroot::copyobject(const T* obj, TString objname)
{
  T* newobj = new T(*obj);
  newobj->SetName(objname);
  return newobj;
}

void xjjroot::dividebinwid(TH1* h)
{
  for(int i=0;i<h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i+1)/h->GetBinWidth(i+1);
      Float_t valErr = h->GetBinError(i+1)/h->GetBinWidth(i+1);
      h->SetBinContent(i+1,val);
      h->SetBinError(i+1,valErr);
    }
}

#endif
