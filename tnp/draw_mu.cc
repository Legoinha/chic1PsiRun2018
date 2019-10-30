#include "tnpcc_tmp.h"

#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "xjjrootuti.h"
#include "xjjcuti.h"

void draw_mu(std::string inputname, std::string dirname, std::string name)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::string tparticle("untitled");
  if(xjjc::str_contains(inputname, "_b.root")) tparticle = "X(3872)";
  if(xjjc::str_contains(inputname, "_a.root")) tparticle = "#psi(2S)";

  TFile* inf = TFile::Open(inputname.c_str());
  fitX::init(inf);
  std::vector<TH1D*> hmupt(tnpcc::nmueta, 0);
  std::map<std::string, std::vector<TH2D*>> hmuscale;
  for(int m=0; m<tnpcc::nmueta; m++)
    {
      hmupt[m] = (TH1D*)inf->Get(Form("hmupt-%d", m));
      xjjroot::setthgrstyle(hmupt[m], kBlack, 20, 1, kBlack, 1, 2);
      xjjroot::sethempty(hmupt[m]);
      for(auto& tt : tnpcc::types)
        {
          hmuscale[tt].push_back((TH2D*)inf->Get(Form("hmuscale-%s-%d", tt.c_str(), m)));
          xjjroot::sethempty(hmuscale[tt][m]);
        }
    }

  xjjroot::setgstyle();
  TCanvas* c = new TCanvas("c", "", 1800, 2400);
  c->Divide(3, 4);
  for(int m=0; m<tnpcc::nmueta; m++)
    {
      c->cd(m+1);
      hmupt[m]->Draw("pe");
      fitX::drawkinematics();
      xjjroot::drawtex(0.90, 0.84-0.04*3, Form("%s < #eta^{#mu} < %s", 
                                               xjjc::number_remove_zero(tnpcc::muetabins[m]).c_str(),
                                               xjjc::number_remove_zero(tnpcc::muetabins[m+1]).c_str()), 0.038, 32, 42);
      xjjroot::drawtex(0.30, 0.25, tparticle.c_str(), 0.038, 12, 62);
      xjjroot::drawCMS("Simulation");
      for(int t=0, n=1; t<tnpcc::types.size(); t++)
        {
          if(tnpcc::types[t]=="nominal" || tnpcc::types[t]=="total") continue;
          c->cd(m+1 + tnpcc::nmueta*n);
          hmuscale[tnpcc::types[t]][m]->Draw("colz");
          fitX::drawkinematics();
          xjjroot::drawtex(0.90, 0.84-0.04*3, Form("%s < #eta^{#mu} < %s", 
                                                   xjjc::number_remove_zero(tnpcc::muetabins[m]).c_str(),
                                                   xjjc::number_remove_zero(tnpcc::muetabins[m+1]).c_str()), 0.038, 32, 42);
          xjjroot::drawtex(0.30, 0.25, tparticle.c_str(), 0.038, 12, 62);
          xjjroot::drawCMS("Simulation");
          n++;
        }
    }
  std::string outputname = "plots/"+dirname+fitX::tagname()+"/drawmu"+name+".pdf";
  xjjroot::mkdir(outputname);
  c->SaveAs(outputname.c_str());
}

int main(int argc, char* argv[])
{
  if(argc==4) { draw_mu(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}
