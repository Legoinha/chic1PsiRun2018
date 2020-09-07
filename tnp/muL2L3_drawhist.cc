#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>

#include "xjjrootuti.h"
#include "xjjcuti.h"
#include <string>
#include "fitX.h"

void draw(std::string inputname, std::string dirname, std::string name)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::string tparticle("untitled");
  if(xjjc::str_contains(inputname, "_a.root")) tparticle = fitX::title_a;
  if(xjjc::str_contains(inputname, "_b.root")) tparticle = fitX::title_b;
  
  TFile* inf = TFile::Open(inputname.c_str());
  fitX::init(inf);
  TH2D* hL2L3pair = (TH2D*)inf->Get("hL2L3pair");
  hL2L3pair->Scale(1./hL2L3pair->Integral());
  
  xjjroot::sethempty(hL2L3pair, 0, 0);
  hL2L3pair->GetXaxis()->SetLabelSize(0.07);
  hL2L3pair->GetYaxis()->SetLabelSize(0.07);
  hL2L3pair->GetXaxis()->SetNdivisions(-303);
  hL2L3pair->GetYaxis()->SetNdivisions(-303);
  hL2L3pair->GetXaxis()->SetBinLabel(1, "Other");
  hL2L3pair->GetXaxis()->SetBinLabel(2, "L2");
  hL2L3pair->GetXaxis()->SetBinLabel(3, "L3");
  hL2L3pair->GetYaxis()->SetBinLabel(1, "Other");
  hL2L3pair->GetYaxis()->SetBinLabel(2, "L2");
  hL2L3pair->GetYaxis()->SetBinLabel(3, "L3");

  xjjroot::setgstyle(1);
  gStyle->SetPaintTextFormat("1.3f");
  hL2L3pair->SetMarkerSize(2);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  c->SetGrid();
  hL2L3pair->Draw("col text");
  xjjroot::drawCMS("Simulation");
  xjjroot::drawtex(0.25, 0.25, tparticle.c_str(), 0.04, 11, 62);

  std::string outputname = dirname+fitX::tagname()+"/muL2L3"+name;
  xjjroot::mkdir("plots/"+outputname+".pdf");
  c->SaveAs(std::string("plots/"+outputname+".pdf").c_str());

  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==4) { draw(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}

