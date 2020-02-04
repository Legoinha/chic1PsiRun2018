#ifndef __FITX_H_
#define __FITX_H_

#include "xjjrootuti.h"
#include "xjjcuti.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <string>

namespace fitX
{
  float ptmincut = 15.;
  float ptmaxcut = 50.;
  float centmincut = 0.;
  float centmaxcut = 90;
  float ymincut = 0;
  float ymaxcut = 1.6;

  void init(float ptmin, float ptmax, float centmin, float centmax, float ymin, float ymax);
  void init(TFile* inf);
  void write();
  std::string tagname() { return std::string(Form("_pt%.0f-%.0f",fitX::ptmincut,fitX::ptmaxcut)) + std::string(Form("_cent%.0f%.0f",fitX::centmincut,fitX::centmaxcut)) + std::string(Form("_y%s-%s",xjjc::number_to_string(fitX::ymincut).c_str(),xjjc::number_to_string(fitX::ymaxcut).c_str())); }
  std::string ytag() { return std::string(Form("%s|y| < %s", (fitX::ymincut?Form("%s < ",xjjc::number_remove_zero(fitX::ymincut).c_str()):""), xjjc::number_remove_zero(fitX::ymaxcut).c_str())); }
  std::string pttag() { return std::string(Form("%s < p_{T} < %s GeV/c", xjjc::number_remove_zero(fitX::ptmincut).c_str(), xjjc::number_remove_zero(fitX::ptmaxcut).c_str())); }
  std::string centtag() { return std::string(Form("Cent. %.0f-%.0f%s", fitX::centmincut, fitX::centmaxcut, "%")); }
  void drawkinematics();
  template<class T> void printhist(T* hh) { std::cout<<"\e[2m"<<hh->GetName()<<"\e[0m\e[36;1m ("<<hh->GetEntries()<<")\e[0m"<<std::endl; }
  template<class T> void setaxis(T* hh);

  Color_t color_data = kRed-3, color_a = kAzure+4, color_b = kGreen-1, color_ss = kGray+1, color_bkg = color_data;
  int ibin_a = 2, ibin_b = 4, nbin = 5;
  std::string title_a = "#psi(2S)", title_b = "X(3872)";
}

void fitX::init(float ptmin, float ptmax, float centmin, float centmax, float ymin, float ymax)
{
  if(ptmin >= 0) ptmincut = ptmin;
  if(ptmax >= 0) ptmaxcut = ptmax;
  if(centmin >= 0) centmincut = centmin;
  if(centmax >= 0) centmaxcut = centmax;
  if(ymin >= 0) ymincut = ymin;
  if(ymax >= 0) ymaxcut = ymax;
}

void fitX::init(TFile* inf)
{
  TTree* kinfo = (TTree*)inf->Get("kinfo");
  float ptmin, ptmax, centmin, centmax, ymin, ymax;
  kinfo->SetBranchAddress("ptmincut", &ptmin);
  kinfo->SetBranchAddress("ptmaxcut", &ptmax);
  kinfo->SetBranchAddress("centmincut", &centmin);
  kinfo->SetBranchAddress("centmaxcut", &centmax);
  kinfo->SetBranchAddress("ymincut", &ymin);
  kinfo->SetBranchAddress("ymaxcut", &ymax);
  kinfo->GetEntry(0);
  init(ptmin, ptmax, centmin, centmax, ymin, ymax);
}

void fitX::write()
{
  TTree* kinfo = new TTree("kinfo", "kinematics");
  kinfo->Branch("ptmincut", &ptmincut);
  kinfo->Branch("ptmaxcut", &ptmaxcut);
  kinfo->Branch("centmincut", &centmincut);
  kinfo->Branch("centmaxcut", &centmaxcut);
  kinfo->Branch("ymincut", &ymincut);
  kinfo->Branch("ymaxcut", &ymaxcut);
  kinfo->Fill();
  kinfo->Write();
}

void fitX::drawkinematics()
{
  xjjroot::drawtex(0.90, 0.84, fitX::pttag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 62);
}

template<class T>
void fitX::setaxis(T* hh)
{
  xjjroot::sethempty(hh, 0, 0.3);
  hh->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hh->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hh->GetXaxis()->SetLabelSize(hh->GetXaxis()->GetLabelSize()*1.5);
}

#endif
