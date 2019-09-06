#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TLegend.h>
#include "xjjrootuti.h"
#include "xjjcuti.h"
#include <string>
#include <algorithm>
#include "tnpcc_tmp.h"
#include "fitX.h"

void draw_tnp_ratio(std::string inputname_a, std::string inputname_b, std::string output)
{
  if(xjjc::str_contains(inputname_b, "_a") && xjjc::str_contains(inputname_a, "_b")) std::swap(inputname_a, inputname_b);
  TFile* inf_a = TFile::Open(inputname_a.c_str());
  fitX::init(inf_a);
  TFile* inf_b = TFile::Open(inputname_b.c_str());
  std::map<std::string, std::map<std::string, TH1D*>> hhp_a, hhp_b;
  std::map<std::string, std::map<std::string, TH1D*>> hhp_ratio;
  std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> ggp_a, ggp_b;
  std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> ggp_ratio;
  for(auto& tt : tnpcc::types)
    {
      hhp_a[tt]["nominal"] = (TH1D*)inf_a->Get(Form("scale_htnp_%s_%s", tt.c_str(), "nominal"));
      hhp_a[tt]["nominal"]->SetName(Form("%s_a", hhp_a[tt]["nominal"]->GetName()));
      xjjroot::sethempty(hhp_a[tt]["nominal"], 0, 0.3);
      for(auto& ee : tnpcc::err)
        {
          ggp_a[tt][ee] = (TGraphAsymmErrors*)inf_a->Get(Form("ggptnp_%s_%s", tt.c_str(), ee.c_str()));
          ggp_a[tt][ee]->SetName(Form("%s_a", ggp_a[tt][ee]->GetName()));
        }
      hhp_b[tt]["nominal"] = (TH1D*)inf_b->Get(Form("scale_htnp_%s_%s", tt.c_str(), "nominal"));
      hhp_b[tt]["nominal"]->SetName(Form("%s_b", hhp_b[tt]["nominal"]->GetName()));
      xjjroot::sethempty(hhp_b[tt]["nominal"], 0, 0.3);
      for(auto& ee : tnpcc::err)
        {
          ggp_b[tt][ee] = (TGraphAsymmErrors*)inf_b->Get(Form("ggptnp_%s_%s", tt.c_str(), ee.c_str()));
          ggp_b[tt][ee]->SetName(Form("%s_b", ggp_b[tt][ee]->GetName()));
        }

      hhp_ratio[tt]["nominal"] = (TH1D*)hhp_b[tt]["nominal"]->Clone(Form("ratio_%s_a", hhp_b[tt]["nominal"]->GetName()));
      hhp_ratio[tt]["nominal"]->Divide(hhp_a[tt]["nominal"]);
      hhp_ratio[tt]["nominal"]->GetYaxis()->SetTitle("#beta^{TnP}_{X(3872)} / #beta^{TnP}_{#psi(2S)}");
      xjjroot::sethempty(hhp_ratio[tt]["nominal"], 0, 0.3);
      xjjroot::setthgrstyle(hhp_ratio[tt]["nominal"], tnpcc::typecolor[tt], 21, 0.8, tnpcc::typecolor[tt], 1, 2, 0, 0, 1001);
    }

  std::vector<double> yy, yye_d, yye_u, xx, xxe;
  for(auto& tt : tnpcc::types)
    {
      for(auto& ee : tnpcc::err)
        {
          xx.clear(); xxe.clear(); yy.clear(); yye_d.clear(); yye_u.clear();
          for(int i=0; i<tnpcc::nptbins; i++)
            {
              double xe = ggp_a[tt][ee]->GetErrorXlow(i);
              double x_a,y_a,y_d_a,y_u_a;
              double x_b,y_b,y_d_b,y_u_b;
              ggp_a[tt][ee]->GetPoint(i, x_a, y_a);
              y_d_a = hhp_a[tt]["nominal"]->GetBinContent(i+1) - ggp_a[tt][ee]->GetErrorYlow(i);
              y_u_a = ggp_a[tt][ee]->GetErrorYhigh(i) - hhp_a[tt]["nominal"]->GetBinContent(i+1);
              ggp_b[tt][ee]->GetPoint(i, x_b, y_b);
              y_d_b = hhp_b[tt]["nominal"]->GetBinContent(i+1) - ggp_b[tt][ee]->GetErrorYlow(i);
              y_u_b = ggp_b[tt][ee]->GetErrorYhigh(i) - hhp_b[tt]["nominal"]->GetBinContent(i+1);

              xx.push_back(x_a);
              xxe.push_back(xe);
              yy.push_back(y_b/y_a);
              yye_d.push_back(TMath::Abs(y_b/y_a - y_d_b/y_d_a));
              yye_u.push_back(TMath::Abs(y_b/y_a - y_u_b/y_u_a));
            }
          ggp_ratio[tt][ee] = new TGraphAsymmErrors(tnpcc::nptbins, xx.data(), yy.data(), xxe.data(), xxe.data(), yye_d.data(), yye_u.data());
          ggp_ratio[tt][ee]->SetName(Form("ggp_ratio_%s_%s", tt.c_str(), ee.c_str()));
          if(ee=="stat") xjjroot::setthgrstyle(ggp_ratio[tt][ee], tnpcc::typecolor[tt], 21, 0.8, 0, 0, 0, tnpcc::typecolor[tt], 0.3, 1001);
          if(ee=="syst") xjjroot::setthgrstyle(ggp_ratio[tt][ee], tnpcc::typecolor[tt], 21, 0.8, tnpcc::typecolor[tt], 1, 1, 0, 0, 1001);
        }
    }

  float ty = 0.85, tx = 0.89;
  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  hhp_ratio["nominal"]["nominal"]->SetMaximum(1.05);
  hhp_ratio["nominal"]["nominal"]->SetMinimum(0.982);
  hhp_ratio["nominal"]["nominal"]->Draw("0");
  for(auto& ee : tnpcc::err)
    {
      std::string drawopt = (ee=="stat"?"2p same":"p5 same");
      ggp_ratio["total"][ee]->Draw(drawopt.c_str());
      ggp_ratio["trg"][ee]->Draw(drawopt.c_str());
      ggp_ratio["muid"][ee]->Draw(drawopt.c_str());
      ggp_ratio["trk"][ee]->Draw(drawopt.c_str());
    }
  hhp_ratio["total"]["nominal"]->Draw("hist same");
  hhp_ratio["nominal"]["nominal"]->Draw("hist same");
  xjjroot::drawCMS("Simulation");
  xjjroot::drawtex(0.89, ty, "PYTHIA + HYDJET", 0.04, 33, 62);
  xjjroot::drawtex(0.22, ty-0.043, "TnP Correction", 0.04, 13, 62);
  TLegend* leg = new TLegend(0.60, ty-0.04-0.045*6, 0.89, ty-0.04, "Stat.          Syst.");
  xjjroot::setleg(leg, 0.04);
  leg->SetNColumns(2);
  for(auto& tt : tnpcc::types)
    {
      if(tt == "nominal") continue;
      for(auto& ee : tnpcc::err)
        {
          leg->AddEntry(ggp_ratio[tt][ee], tt.c_str(), "f");
        }
    }
  leg->AddEntry(hhp_ratio["nominal"]["nominal"], "Nominal", "l");
  leg->Draw();
  gPad->RedrawAxis();
  std::string outputname = "plots/"+output+fitX::tagname()+"/drawtnp_ratio.pdf";
  xjjroot::mkdir(outputname);
  c->SaveAs(outputname.c_str());  

  // print
  float per_a_u, per_a_d, per_b_u, per_b_d, per_ab_u, per_ab_d;
  for(int i=0; i<tnpcc::nptbins; i++)
    {
      double ix,iy;
      ggp_a["total"]["stat"]->GetPoint(i, ix, iy);
      float per_a_u_stat = ggp_a["total"]["stat"]->GetErrorYhigh(i)/iy;
      float per_a_d_stat = ggp_a["total"]["stat"]->GetErrorYlow(i)/iy;
      ggp_a["total"]["syst"]->GetPoint(i, ix, iy);
      float per_a_u_syst = ggp_a["total"]["syst"]->GetErrorYhigh(i)/iy;
      float per_a_d_syst = ggp_a["total"]["syst"]->GetErrorYlow(i)/iy;
      per_a_u = TMath::Sqrt(per_a_u_stat*per_a_u_stat+per_a_u_syst*per_a_u_syst);
      per_a_d = TMath::Sqrt(per_a_d_stat*per_a_d_stat+per_a_d_syst*per_a_d_syst);

      ggp_b["total"]["stat"]->GetPoint(i, ix, iy);
      float per_b_u_stat = ggp_b["total"]["stat"]->GetErrorYhigh(i)/iy;
      float per_b_d_stat = ggp_b["total"]["stat"]->GetErrorYlow(i)/iy;
      ggp_b["total"]["syst"]->GetPoint(i, ix, iy);
      float per_b_u_syst = ggp_b["total"]["syst"]->GetErrorYhigh(i)/iy;
      float per_b_d_syst = ggp_b["total"]["syst"]->GetErrorYlow(i)/iy;
      per_b_u = TMath::Sqrt(per_b_u_stat*per_b_u_stat+per_b_u_syst*per_b_u_syst);
      per_b_d = TMath::Sqrt(per_b_d_stat*per_b_d_stat+per_b_d_syst*per_b_d_syst);

      ggp_ratio["total"]["stat"]->GetPoint(i, ix, iy);
      float per_ab_u_stat = ggp_ratio["total"]["stat"]->GetErrorYhigh(i)/iy;
      float per_ab_d_stat = ggp_ratio["total"]["stat"]->GetErrorYlow(i)/iy;
      ggp_ratio["total"]["syst"]->GetPoint(i, ix, iy);
      float per_ab_u_syst = ggp_ratio["total"]["syst"]->GetErrorYhigh(i)/iy;
      float per_ab_d_syst = ggp_ratio["total"]["syst"]->GetErrorYlow(i)/iy;
      per_ab_u = TMath::Sqrt(per_ab_u_stat*per_ab_u_stat+per_ab_u_syst*per_ab_u_syst);
      per_ab_d = TMath::Sqrt(per_ab_d_stat*per_ab_d_stat+per_ab_d_syst*per_ab_d_syst);
    }
  std::cout<<"TnP & "<<Form("+%.1f", per_a_u*1.e+2)<<"\\% & "<<Form("+%.1f", per_b_u*1.e+2)<<"\\% & "<<Form("+%.1f", per_ab_u*1.e+2)<<"\\% \\\\"<<std::endl;
  std::cout<<" & "<<Form("-%.1f", per_a_d*1.e+2)<<"\\% & "<<Form("-%.1f", per_b_d*1.e+2)<<"\\% & "<<Form("-%.1f", per_ab_d*1.e+2)<<"\\% \\\\"<<std::endl;

}

int main(int argc, char* argv[])
{
  if(argc==4) { draw_tnp_ratio(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}
