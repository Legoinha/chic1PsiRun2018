#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include <vector>
#include <string>
#include <iostream>

#include "fitX.h"
#include "xjjrootuti.h"

const int n = 50000;
void drawkinematics();
void acc_fit(std::string input, std::string output, std::string type)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  // output += (fitX::tagname()+"/"+type);
  TH1F* hratiodis_a = (TH1F*)inf->Get("hratiodis_a");
  xjjroot::sethempty(hratiodis_a, 0, 0);
  TH1F* hratiodis_b = (TH1F*)inf->Get("hratiodis_b");
  xjjroot::sethempty(hratiodis_b, 0, 0);

  int nbin = hratiodis_a->GetXaxis()->GetNbins();
  TH1F* h_a = (TH1F*)hratiodis_a->Clone("h_a");
  TH1F* h_b = (TH1F*)hratiodis_b->Clone("h_b");
  std::vector<float> means_a(nbin), means_b(nbin), errs_a(nbin), errs_b(nbin);
  for(int k=0; k<nbin; k++) 
    { 
      means_a[k] = hratiodis_a->GetBinContent(k+1);
      means_b[k] = hratiodis_b->GetBinContent(k+1);
      errs_a[k] = hratiodis_a->GetBinError(k+1);
      errs_b[k] = hratiodis_b->GetBinError(k+1);
    }
  float xmin_a = h_a->GetBinCenter(1) - h_a->GetBinWidth(1)/2.;
  float xmax_a = h_a->GetBinCenter(nbin) + h_a->GetBinWidth(nbin)/2.;
  TF1* f_a = new TF1("f_a", "exp([0]+[1]*x)", xmin_a, xmax_a);
  // f_a->SetParameters(1, 0);
  f_a->SetParameters(0, 0);
  xjjroot::settfstyle(f_a, xjjroot::mycolor_middle["red"], 3, 4);
  float xmin_b = h_b->GetBinCenter(1) - h_b->GetBinWidth(1)/2.;
  float xmax_b = h_b->GetBinCenter(nbin) + h_b->GetBinWidth(nbin)/2.;
  TF1* f_b = new TF1("f_b", "exp([0]+[1]*x)", xmin_b, xmax_b);
  // f_b->SetParameters(1, 0);
  f_b->SetParameters(0, 0);
  xjjroot::settfstyle(f_b, xjjroot::mycolor_middle["red"], 3, 4);

  TRandom3* rmd = new TRandom3();

  std::string outputname = "rootfiles/"+output+"/acc_fit.root";
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  TTree* par = new TTree("par", "");
  float a0; par->Branch("a0", &a0);
  float a1; par->Branch("a1", &a1);
  float b0; par->Branch("b0", &b0);
  float b1; par->Branch("b1", &b1);
  
  xjjroot::setgstyle(1);
  std::string outputnamec = "plots/"+output+"/accidx/cacc";
  xjjroot::mkdir(outputnamec);
  for(int i=0, ic=0; i<n; i++)
    {
      f_a->SetParameters(0, 0);
      f_b->SetParameters(0, 0);
      for(int k=0; k<nbin; k++)
        {
          h_a->SetBinContent(k+1, rmd->Gaus(means_a[k], errs_a[k]));
          h_b->SetBinContent(k+1, rmd->Gaus(means_b[k], errs_b[k]));
        }
      h_a->Fit("f_a", "Nq", "", xmin_a, xmax_a);
      h_a->Fit("f_a", "Nq", "", xmin_a, xmax_a);
      a0 = f_a->GetParameter(0);
      a1 = f_a->GetParameter(1);
      h_b->Fit("f_b", "Nq", "", xmin_b, xmax_b);
      h_b->Fit("f_b", Form("N%s", i%(n/5)!=0?"q":""), "", xmin_b, xmax_b);
      b0 = f_b->GetParameter(0);
      b1 = f_b->GetParameter(1);
      if(i%(n/5)==0) 
        {
          TCanvas* c = new TCanvas("c", "", 1200, 600);
          c->Divide(2, 1);
          c->cd(1);
          h_a->Draw("pe");
          f_a->Draw("same");
          xjjroot::drawline(xmin_a, 1, xmax_a, 1, fitX::color_a, 7, 2, 0.3);
          xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.043, 12, 62, fitX::color_a);
          drawkinematics();
          xjjroot::drawCMS();
          xjjroot::drawcomment(Form("Toy MC --> %d", i));
          c->cd(2);
          h_b->Draw("pe");
          f_b->Draw("same");
          xjjroot::drawline(xmin_b, 1, xmax_b, 1, fitX::color_b, 7, 2, 0.3);
          xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.043, 12, 62, fitX::color_b);
          drawkinematics();
          xjjroot::drawCMS();
          c->SaveAs(Form("%s-%d.pdf", outputnamec.c_str(), ic)); 
          delete c;
          ic++; 
        }
      par->Fill();
    }
  par->Write();
  fitX::write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==4) { acc_fit(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}


void drawkinematics()
{
  xjjroot::drawtex(0.24, 0.84-0.05, "Pseudo-data", 0.043, 12, 42, kBlack);
  xjjroot::drawtex(0.90, 0.84, fitX::pttag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 62);
}
