#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include <vector>
#include <string>
#include <iostream>

#include "fitX.h"
#include "xjjrootuti.h"

void drawkinematics();
void eff_fit(std::string input, std::string output, std::string type)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  // output += (fitX::tagname()+"/"+type);
  TH1F* hratiodis_a = (TH1F*)inf->Get("hratiodis_a");
  xjjroot::sethempty(hratiodis_a, 0, 0);
  TH1F* hratiodis_b = (TH1F*)inf->Get("hratiodis_b");
  xjjroot::sethempty(hratiodis_b, 0, 0);
  TGraphAsymmErrors* gratiodis_a = (TGraphAsymmErrors*)inf->Get("gratiodis_a");
  TGraphAsymmErrors* gratiodis_b = (TGraphAsymmErrors*)inf->Get("gratiodis_b");

  int nbin = hratiodis_a->GetXaxis()->GetNbins();
  TH1F* h_a = (TH1F*)hratiodis_a->Clone("h_a");
  TH1F* h_b = (TH1F*)hratiodis_b->Clone("h_b");
  TGraphAsymmErrors* g_a = (TGraphAsymmErrors*)gratiodis_a->Clone("g_a");
  TGraphAsymmErrors* g_b = (TGraphAsymmErrors*)gratiodis_b->Clone("g_b");
  // std::vector<double> gx(nbin), gy(nbin);
  // for(int i=0; i<nbin; i++) { g_a->GetPoint(i, gx[i], gy[i]); }
  // std::vector<float> means_a(nbin), means_b(nbin), errs_a(nbin), errs_b(nbin);
  // for(int k=0; k<nbin; k++) 
  //   { 
  //     means_a[k] = hratiodis_a->GetBinContent(k+1);
  //     means_b[k] = hratiodis_b->GetBinContent(k+1);
  //     errs_a[k] = hratiodis_a->GetBinError(k+1);
  //     errs_b[k] = hratiodis_b->GetBinError(k+1);
  //   }
  float xmin_a = h_a ->GetBinCenter(1)    - h_a ->GetBinWidth(1)/2.;
  float xmax_a = h_a ->GetBinCenter(nbin) + h_a ->GetBinWidth(nbin)/2.;
  float xmin_b = h_b ->GetBinCenter(1)    - h_b ->GetBinWidth(1)/2.;
  float xmax_b = h_b ->GetBinCenter(nbin) + h_b ->GetBinWidth(nbin)/2.;

  std::map<std::string, std::vector<double>> parms = {
    // std::pair<std::string, std::vector<double>>("[0]+[1]*x"            , std::vector<double>({1, 0})),
    // std::pair<std::string, std::vector<double>>("[0]+[1]*x+[2]*x*x"    , std::vector<double>({1, 0, 0})),
    std::pair<std::string, std::vector<double>>("exp([0]+[1]*x)"       , std::vector<double>({0, 0})),
    // std::pair<std::string, std::vector<double>>("[0]+[1]*x+[2]/x"      , std::vector<double>({1, 0, 0})),
    // std::pair<std::string, std::vector<double>>("[0]+[1]/sqrt(x)+[2]/x", std::vector<double>({1, 0, 0})),    
    // std::pair<std::string, std::vector<double>>("[0]+[1]/x"            , std::vector<double>({1, 0})),
  };

  xjjroot::setgstyle(1);
  std::map<std::string, TF1*> f_a, f_b;
  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  h_a->Draw("AXIS");
  g_a->Draw("pe same");
  xjjroot::drawline(xmin_a, 1, xmax_a, 1, kGray+2, 7, 2, 0.3);
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.043, 12, 62);
  drawkinematics();
  xjjroot::drawCMS();
  c->cd(2);
  h_b->Draw("AXIS");
  g_b->Draw("pe same");
  xjjroot::drawline(xmin_b, 1, xmax_b, 1, kGray+2, 7, 2, 0.3);
  xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.043, 12, 62);
  drawkinematics();
  xjjroot::drawCMS();
  int n=0;
  for(auto& p : parms)
    {
      std::string skey(p.first);
      f_a[skey] = new TF1(Form("f_a(%s)", skey.c_str()), skey.c_str(), xmin_a, xmax_a);
      f_b[skey] = new TF1(Form("f_b(%s)", skey.c_str()), skey.c_str(), xmin_b, xmax_b);

      for(int k=0; k<p.second.size(); k++) { f_a[skey]->SetParameter(k+1, (p.second)[k]); }
      for(int k=0; k<p.second.size(); k++) { f_b[skey]->SetParameter(k+1, (p.second)[k]); }

      xjjroot::settfstyle(f_a[skey], xjjroot::mycolor_middle[xjjroot::cc[n]], 3, 4);
      xjjroot::settfstyle(f_b[skey], xjjroot::mycolor_middle[xjjroot::cc[n]], 3, 4);

      std::cout<<std::endl<<"---> "<<skey<<std::endl;

      c->cd(1);
      g_a->Fit(Form("f_a(%s)", skey.c_str()), "Nq", "", xmin_a, xmax_a);
      g_a->Fit(Form("f_a(%s)", skey.c_str()), "N", "", xmin_a, xmax_a);
      f_a[skey]->Draw("same");
      xjjroot::drawtex(0.35, 0.75-n*0.04, skey.c_str(), 0.036, 13, 62, xjjroot::mycolor_middle[xjjroot::cc[n]]);
      xjjroot::drawtex(0.33, 0.75-n*0.04, Form("%.0f%s", TMath::Prob(f_a[skey]->GetChisquare(), f_a[skey]->GetNDF())*1.e+2, "%"), 0.036, 33, 62, xjjroot::mycolor_middle[xjjroot::cc[n]]);
      c->cd(2);
      g_b->Fit(Form("f_b(%s)", skey.c_str()), "Nq", "", xmin_b, xmax_b);
      g_b->Fit(Form("f_b(%s)", skey.c_str()), "N", "", xmin_b, xmax_b);
      f_b[skey]->Draw("same");
      xjjroot::drawtex(0.35, 0.75-n*0.04, skey.c_str(), 0.036, 13, 62, xjjroot::mycolor_middle[xjjroot::cc[n]]);
      xjjroot::drawtex(0.33, 0.75-n*0.04, Form("%.0f%s", TMath::Prob(f_b[skey]->GetChisquare(), f_b[skey]->GetNDF())*1.e+2, "%"), 0.036, 33, 62, xjjroot::mycolor_middle[xjjroot::cc[n]]);

      n++;      
    }
  std::string outputnamec = "plots/"+output+"/efficiency/cefffit.pdf";
  xjjroot::mkdir(outputnamec);
  c->SaveAs(outputnamec.c_str());

  std::vector<TString> par_a(n), par_b(n), fun_a(n), fun_b(n);
  std::string outputname = "rootfiles/"+output+"/eff_fit.root";
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  fitX::write();
  outf->cd();
  TTree* ntf = new TTree("ntf", "functions");
  ntf->Branch("n", &n);
  do {
    int i=0;
    for(auto& ff : f_a)
      {
        par_a[i] = ff.second->GetExpFormula("CLING P");
        fun_a[i] = ff.second->GetExpFormula();
        fun_a[i].ReplaceAll("[p", "[");
        ntf->Branch(Form("par_a-%d", i+1), &(par_a[i]));
        ntf->Branch(Form("fun_a-%d", i+1), &(fun_a[i]));
        i++;
      }
  } while(false);
  do {
    int i=0;
    for(auto& ff : f_b)
      {
        par_b[i] = ff.second->GetExpFormula("CLING P");
        fun_b[i] = ff.second->GetExpFormula();
        fun_b[i].ReplaceAll("[p", "[");
        ntf->Branch(Form("par_b-%d", i+1), &(par_b[i]));
        ntf->Branch(Form("fun_b-%d", i+1), &(fun_b[i]));
        i++;
      }
  } while(false);
  ntf->Fill();
  ntf->Write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==4) { eff_fit(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}

void drawkinematics()
{
  xjjroot::drawtex(0.24, 0.84-0.05, "Pseudo-data", 0.043, 12, 42, kBlack);
  xjjroot::drawtex(0.90, 0.84, fitX::pttag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 62);
}

