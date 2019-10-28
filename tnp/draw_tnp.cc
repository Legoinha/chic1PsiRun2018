#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TLegend.h>
#include "xjjrootuti.h"
#include "xjjcuti.h"
#include <string>
#include "tnpcc_tmp.h"
#include "fitX.h"

void draw_tnp(std::vector<std::string> inputname, std::string dirname, std::string name)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if(inputname.size() < 1) return;
  int n = inputname.size();
  std::string tparticle("untitled");
  if(xjjc::str_contains(inputname[0], "_b.root")) tparticle = "X(3872)";
  if(xjjc::str_contains(inputname[0], "_a.root")) tparticle = "#psi(2S)";
  
  std::vector<tnpcc::drawtnp*> dd(n);
  std::vector<std::string> tname(n);
  for(int i=0; i<n; i++) { 
    dd[i] = new tnpcc::drawtnp(inputname[i].c_str(), Form("_%d", i+1), (n==1?0:xjjroot::mycolor_middle[xjjroot::cc[i]])); 
    tname[i] = n==1?"":((TH2F*)(dd[i]->inf()->Get("hptweight")))->GetYaxis()->GetTitle();
  }
  fitX::init(dd[0]->inf());

  // draw
  xjjroot::setgstyle(1);
  if(n==1)
    {
      float ty = 0.85, tx = 0.89;
      dd[0]->setleg(ty);
      TCanvas* c = new TCanvas("c", "", 1200, 600);
      c->Divide(2, 1);
      c->cd(1);
      // gPad->SetLogy();
      dd[0]->draw();
      xjjroot::drawCMS("Simulation");
      gPad->RedrawAxis();
      c->cd(2);
      dd[0]->drawscale();
      xjjroot::drawCMS("Simulation");
      gPad->RedrawAxis();

      c->cd(1);
      xjjroot::drawtex(0.22, ty, tparticle.c_str(), 0.04, 13, 62);
      dd[0]->drawkinematics(ty);
      dd[0]->drawleg();
      c->cd(2);
      xjjroot::drawtex(0.22, ty, tparticle.c_str(), 0.04, 13, 62);
      dd[0]->drawkinematics(ty);
      dd[0]->drawleg();
      c->cd();

      std::string outputname = dirname+fitX::tagname()+"/drawtnp"+name;
      xjjroot::mkdir("plots/"+outputname+".pdf");
      c->SaveAs(std::string("plots/"+outputname+".pdf").c_str());

      xjjroot::mkdir("rootfiles/"+outputname+".root");
      TFile* outf = new TFile(std::string("rootfiles/"+outputname+".root").c_str(), "recreate");
      for(auto& hp : dd[0]->hhp) (hp.second)["nominal"]->Write();
      for(auto& gp : dd[0]->ggp) 
        {
          for(auto& gr : gp.second)
            gr.second->Write();
        }
      fitX::write();
      outf->Close();

      float stat_u = dd[0]->ggp["total"]["stat"]->GetErrorYhigh(0)/dd[0]->hhp["total"]["nominal"]->GetBinContent(1);
      float stat_d = dd[0]->ggp["total"]["stat"]->GetErrorYlow(0)/dd[0]->hhp["total"]["nominal"]->GetBinContent(1);
      float syst_u = dd[0]->ggp["total"]["syst"]->GetErrorYhigh(0)/dd[0]->hhp["total"]["nominal"]->GetBinContent(1);
      float syst_d = dd[0]->ggp["total"]["syst"]->GetErrorYlow(0)/dd[0]->hhp["total"]["nominal"]->GetBinContent(1);
      float err_u = sqrt(stat_u*stat_u + syst_u*syst_u);
      float err_d = sqrt(stat_d*stat_d + syst_d*syst_d);
      std::cout<<"err_u "<<err_u<<std::endl;
      std::cout<<"err_d "<<err_d<<std::endl;
    }
  else
    {
      float ty = 0.85, tx = 0.89;
      TLegend* leg = new TLegend(0.59, ty-0.04-0.05*n, 0.79, ty-0.04, "p_{T} weight");
      xjjroot::setleg(leg, 0.04);

      TCanvas* c = new TCanvas("c", "", 600, 600);
      dd[0]->hhp["nominal"]["nominal"]->SetLineColor(kBlack);
      dd[0]->hhp["nominal"]["nominal"]->SetLineWidth(3);
      dd[0]->hhp["nominal"]["nominal"]->Draw("0");
      for(int i=0; i<n; i++)
        {
          // dd[i]->ggp["total"]["stat"]->Draw("2 same");
          dd[i]->ggp["total"]["syst"]->Draw("5p same");
          leg->AddEntry(dd[i]->ggp["total"]["syst"], tname[i].c_str(), "f");
        }
      leg->Draw();
      xjjroot::drawtex(0.22, ty, tparticle.c_str(), 0.04, 13, 62);
      dd[0]->drawkinematics(ty);
      xjjroot::drawCMS("Simulation");
      gPad->RedrawAxis();

      std::string outputname = dirname+fitX::tagname()+"/funs/drawtnp"+name;
      xjjroot::mkdir("plots/"+outputname+".pdf");
      c->SaveAs(std::string("plots/"+outputname+".pdf").c_str());

      xjjroot::mkdir("rootfiles/"+outputname+".root");
      TFile* outf = new TFile(std::string("rootfiles/"+outputname+".root").c_str(), "recreate");
      for(auto& d : dd)
        {
          for(auto& hp : d->hhp) (hp.second)["nominal"]->Write();
          for(auto& gp : d->ggp) 
            {
              for(auto& gr : gp.second)
                gr.second->Write();
            }
        }
      fitX::write();
      outf->Close();
    }

  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc>=4) { 
    std::vector<std::string> inputs(argc-3);
    for(int i=0; i<argc-3; i++) { inputs[i] = argv[i+1]; }
    draw_tnp(inputs, argv[argc-2], argv[argc-1]); 
    return 0; }
  return 1;
}

