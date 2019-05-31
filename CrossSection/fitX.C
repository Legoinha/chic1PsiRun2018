#include "xjjrootuti.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>

#include <string>

void fitX(std::string input, std::string cut)
{
  xjjroot::setgstyle(1);

  TFile* inf = TFile::Open(input.c_str());
  TTree* ntmix = (TTree*)inf->Get("Bfinder/ntmix");
  ntmix->AddFriend("hiEvtAnalyzer/HiTree");
  ntmix->AddFriend("dataset/mva");

  TH1F* h = new TH1F("h", "", 40, 3.6, 4);
  ntmix->Project("h", "Bmass", TCut(cut.c_str())&&"hiBin<=200");
  // TH1F* hmc_a = new TH1F("hmc_a", "", 40, 3.6, 4);
  // ntmix->Project("h", "Bmass", TCut(cut.c_str())&&"hiBin<=200");
  TCanvas* c = new TCanvas("c", "", 800, 600);
  xjjroot::sethempty(h);
  h->SetMinimum(0);
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
  h->Draw("pe");
  TF1 *f = new TF1("f","[0]+[1]*x+[2]*x*x+[8]*x*x*x+[9]*x*x*x*x+[3]*TMath::Gaus(x,[4],[5])/(TMath::Sqrt(2*3.14159)*[5])+[6]*TMath::Gaus(x,[7],[5])/(TMath::Sqrt(2*3.14159)*[5])");
  f->SetLineWidth(3);
  f->SetLineColor(kRed);
  f->SetParameters(0,0,0,0,3.686,0.00357,1,3.8725,0.0054);
  f->FixParameter(4,3.68610);
  f->FixParameter(5,0.00357);
  f->FixParameter(7,3.87169);
  f->SetNpx(2000);
  h->Sumw2();
  h->Fit("f","");
  h->Fit("f","LL");
  h->Fit("f","LL");
  h->Fit("f","LL");

  TF1 *fsig1 = new TF1("fsig1", "[0]+[1]*x+[2]*x*x+[8]*x*x*x+[9]*x*x*x*x+[3]*TMath::Gaus(x,[4],[5])/(TMath::Sqrt(2*3.14159)*[5])+[6]*TMath::Gaus(x,[7],[5])/(TMath::Sqrt(2*3.14159)*[5])", 3.6, 4.0);
  fsig1->SetParameters(0, 0, 0, f->GetParameter(3), f->GetParameter(4), f->GetParameter(5), 0, 0, 0, 0);
  fsig1->SetParError(3, f->GetParError(3));
  fsig1->SetLineColor(kOrange-3);
  fsig1->SetFillColor(kOrange-3);
  fsig1->SetFillStyle(3002);
  fsig1->SetLineWidth(3);
  fsig1->SetNpx(2000);
  fsig1->Draw("same");
  std::cout<<fsig1->GetParameter(3)<<" "<<fsig1->GetParError(3)<<std::endl;

  TF1 *fsig2 = new TF1("fsig2", "[0]+[1]*x+[2]*x*x+[8]*x*x*x+[9]*x*x*x*x+[3]*TMath::Gaus(x,[4],[5])/(TMath::Sqrt(2*3.14159)*[5])+[6]*TMath::Gaus(x,[7],[5])/(TMath::Sqrt(2*3.14159)*[5])", 3.6, 4.0);
  fsig2->SetParameters(0, 0, 0, 0, 0, f->GetParameter(5), f->GetParameter(6), f->GetParameter(7), 0, 0);
  fsig2->SetParError(6, f->GetParError(6));
  fsig2->SetLineColor(kGreen+3);
  fsig2->SetFillColor(kGreen+3);
  fsig2->SetFillStyle(3002);
  fsig2->SetLineWidth(3);
  fsig2->SetNpx(2000);
  fsig2->Draw("same");
  std::cout<<fsig2->GetParameter(6)<<" "<<fsig2->GetParError(6)<<std::endl;
  c->RedrawAxis();

  float ysig1 = fsig1->Integral(3.6, 4.0)/0.01;
  float ysig2 = fsig2->Integral(3.6, 4.0)/0.01;

  xjjroot::drawtex(0.23, 0.83, Form("N_{#Psi'} = %.1f #pm %.1f", ysig1, fsig1->GetParError(3)*ysig1*0.01/fsig1->GetParameter(3)), 0.042, 12, 62, kOrange-3);
  xjjroot::drawtex(0.23, 0.78, Form("N_{X} = %.1f #pm %.1f", ysig2, fsig2->GetParError(6)*ysig2*0.01/fsig2->GetParameter(6)), 0.042, 12, 62, kGreen+3);

  c->SaveAs("plots/chmass.pdf");
}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX(argv[1], argv[2]); return 0; }
  return 1;
}
