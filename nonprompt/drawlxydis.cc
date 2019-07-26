#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "lxydis.h"
#include "fitX.h"

void drawkinematic();

void drawlxydis(std::string input, std::string output)
{
  std::string plotdirname(xjjc::str_replaceall(xjjc::str_getdir(output), "rootfiles/", "plots/"));

  TFile* inf = new TFile(Form("%s.root", input.c_str()));
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
  TH1F* hlxymcnp_a = (TH1F*)inf->Get("hlxymcnp_a");
  TH1F* hlxymcnp_b = (TH1F*)inf->Get("hlxymcnp_b");
  TH1F* h = (TH1F*)inf->Get("h");
  std::vector<TH1F*> hBenr(lxydis::lxycut.size());
  for(int k=0; k<lxydis::lxycut.size(); k++)
    { hBenr[k] = (TH1F*)inf->Get(Form("hBenr_%d", k)); }
  // fit + yield
  TH1F* hyield_a = new TH1F("hyield_a", "", 5, 0, 5);
  TH1F* hyield_b = new TH1F("hyield_b", "", 5, 0, 5);
  std::vector<TH1F*> hBenryield_a(lxydis::lxycut.size()), hBenryield_b(lxydis::lxycut.size());
  for(int k=0; k<lxydis::lxycut.size(); k++)
    {
      hBenryield_a[k] = new TH1F(Form("hBenryield_a_%d", k), "", 5, 0, 5);
      hBenryield_b[k] = new TH1F(Form("hBenryield_b_%d", k), "", 5, 0, 5);
    }
  std::vector<TH1F*> vh({h}); for(auto& hh : hBenr) { vh.push_back(hh); } ;
  std::vector<TH1F*> vhyield_a({hyield_a}); for(auto& hh : hBenryield_a) { vhyield_a.push_back(hh); } ;
  std::vector<TH1F*> vhyield_b({hyield_b}); for(auto& hh : hBenryield_b) { vhyield_b.push_back(hh); } ;
  std::vector<Style_t> mstyle({20}); for(auto& cc: lxydis::lxycut) { mstyle.push_back(24); } ;
  std::vector<Color_t> mcolor({kGray+3, (Color_t)xjjroot::mycolor_middle["azure"], (Color_t)xjjroot::mycolor_middle["red"], (Color_t)xjjroot::mycolor_middle["green"], 
        (Color_t)xjjroot::mycolor_middle["orange"], (Color_t)xjjroot::mycolor_middle["cyan"], (Color_t)xjjroot::mycolor_middle["magenta"]});
  float ymax = 0;
  for(int l=0; l<vh.size(); l++)
    {
      // ====>
      // std::vector<TF1*> funsft = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, "plots/fltm", false, false); // fix mean = false
      std::vector<TF1*> funs   = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, plotdirname+"/idx/", true, true, Form("_%d", l)); // fix mean = true
      xjjroot::setgstyle();
      // <====
      float ysig1 = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
      float ysig1err = funs[0]->GetParError(5)*ysig1/funs[0]->GetParameter(5);
      float ysig2 = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
      float ysig2err = funs[0]->GetParError(10)*ysig2/funs[0]->GetParameter(10);

      // yield
      xjjroot::setthgrstyle(vhyield_a[l], mcolor[l], mstyle[l], 1.2, mcolor[l], 1, 2, mcolor[l], 0.1, 1001);
      vhyield_a[l]->SetBinContent(fitX::ibin_a, ysig1);
      vhyield_a[l]->SetBinError(fitX::ibin_a, ysig1err);
      xjjroot::setthgrstyle(vhyield_b[l], mcolor[l], mstyle[l], 1.2, mcolor[l], 1, 2, mcolor[l], 0.1, 1001);
      vhyield_b[l]->SetBinContent(fitX::ibin_b, ysig2);
      vhyield_b[l]->SetBinError(fitX::ibin_b, ysig2err);

      if(l && (vhyield_a[l]->GetMaximum()>ymax)) { ymax = vhyield_a[l]->GetMaximum(); }
    } 
  ymax *= 2.;
  TH2F* hempty = new TH2F("hempty", ";;Raw Yield", 5, 0, 5, 10, 0, hyield_a->GetMaximum()*2.);
  xjjroot::sethempty(hempty, 0, 0.3);
  hempty->GetXaxis()->SetBinLabel(fitX::ibin_a, "#psi(2S)");
  hempty->GetXaxis()->SetBinLabel(fitX::ibin_b, "X(3872)");
  hempty->GetXaxis()->SetLabelSize(hempty->GetXaxis()->GetLabelSize()*1.5);
  TH2F* hemptyBenr = new TH2F("hemptyBenr", ";;Raw Yield", 5, 0, 5, 10, 0, ymax);
  xjjroot::sethempty(hemptyBenr, 0, 0.3);
  hemptyBenr->GetXaxis()->SetBinLabel(fitX::ibin_a, "#psi(2S)");
  hemptyBenr->GetXaxis()->SetBinLabel(fitX::ibin_b, "X(3872)");
  hemptyBenr->GetXaxis()->SetLabelSize(hemptyBenr->GetXaxis()->GetLabelSize()*1.5);

  TLegend* leg = new TLegend(0.60, 0.60-0.045*lxydis::lxycut.size(), 0.90, 0.60);
  xjjroot::setleg(leg, 0.04);
  for(int k=0; k<lxydis::lxycut.size(); k++) { leg->AddEntry(hBenryield_a[k], Form("l_{xy} > %s mm", xjjc::number_remove_zero(lxydis::lxycut[k]).c_str()), "lp"); }
  xjjroot::setgstyle();
  TCanvas* cy = new TCanvas("cy", "", 1200, 600);
  cy->Divide(2, 1);
  cy->cd(1);
  hempty->Draw();
  hyield_a->Draw("ple same");
  hyield_b->Draw("ple same");
  drawkinematic();
  xjjroot::drawtex(0.22, 0.84, "Inclusive", 0.042, 12, 62);  
  xjjroot::drawCMS();
  cy->cd(2);
  hemptyBenr->Draw();
  for(auto& hh : hBenryield_a) { hh->Draw("ple same"); }
  for(auto& hh : hBenryield_b) { hh->Draw("ple same"); }
  leg->Draw();
  drawkinematic();
  xjjroot::drawtex(0.22, 0.84, "B-enriched", 0.042, 12, 62);  
  xjjroot::drawCMS();
  cy->SaveAs(Form("%s/cfpromptyield.pdf", plotdirname.c_str()));

  // fprompt
  hlxymcnp_a->Scale(1./hlxymcnp_a->Integral(), "width");
  hlxymcnp_b->Scale(1./hlxymcnp_b->Integral(), "width");
  std::vector<TH1F*> hlxyfrac_a(lxydis::lxycut.size()), hlxyfrac_b(lxydis::lxycut.size()), hyieldprompt_a(lxydis::lxycut.size()), hyieldprompt_b(lxydis::lxycut.size());
  std::vector<TEfficiency*> grfprompt_a(lxydis::lxycut.size()), grfprompt_b(lxydis::lxycut.size());
  for(int k=0; k<lxydis::lxycut.size(); k++)
    {
      std::vector<double> vlxyfrac = lxydis::nplxyfrac(hlxymcnp_a, hlxymcnp_b, lxydis::lxycut[k]);
      hlxyfrac_a[k] = new TH1F(Form("hlxyfrac_a_%d", k), "", 5, 0, 5);
      hlxyfrac_a[k]->SetBinContent(fitX::ibin_a, vlxyfrac[0]);
      hlxyfrac_a[k]->SetBinError(fitX::ibin_a, vlxyfrac[1]);
      hlxyfrac_b[k] = new TH1F(Form("hlxyfrac_b_%d", k), "", 5, 0, 5);
      hlxyfrac_b[k]->SetBinContent(fitX::ibin_b, vlxyfrac[2]);
      hlxyfrac_b[k]->SetBinError(fitX::ibin_b, vlxyfrac[3]);
      grfprompt_a[k] = lxydis::calclxyfprompt(hyield_a, hBenryield_a[k], hlxyfrac_a[k], Form("grfprompt_a_%d", k), &(hyieldprompt_a[k]));
      grfprompt_b[k] = lxydis::calclxyfprompt(hyield_b, hBenryield_b[k], hlxyfrac_b[k], Form("grfprompt_b_%d", k), &(hyieldprompt_b[k]));
      xjjroot::setthgrstyle(grfprompt_a[k], mcolor[k+1], 24, 1.5, mcolor[k+1], 1, 3);
      xjjroot::setthgrstyle(grfprompt_b[k], mcolor[k+1], 24, 1.5, mcolor[k+1], 1, 3);
    }
  xjjroot::setgstyle(2);
  TCanvas* cfprompt = new TCanvas("cfprompt", "", 600, 600);
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";;f_{prompt}", 5, 0, 5, 10, 0, 1.3);
  xjjroot::sethempty(hemptyfprompt, 0, 0.3);
  hemptyfprompt->GetXaxis()->SetBinLabel(fitX::ibin_a, "#psi(2S)");
  hemptyfprompt->GetXaxis()->SetBinLabel(fitX::ibin_b, "X(3872)");
  hemptyfprompt->GetXaxis()->SetLabelSize(hemptyfprompt->GetXaxis()->GetLabelSize()*1.5);
  hemptyfprompt->Draw();
  xjjroot::drawline(0, 1, 5, 1, kGray+1, 9, 2);
  for(auto& gg : grfprompt_a) gg->Draw("pe same");
  for(auto& gg : grfprompt_b) gg->Draw("pe same");
  leg->Draw();
  drawkinematic();
  xjjroot::drawCMS();
  cfprompt->SaveAs(Form("%s/cfprompt.pdf", plotdirname.c_str()));

  xjjroot::mkdir(output);
  TFile* outf = new TFile(output.c_str(), "recreate");
  outf->cd();
  for(auto& gg : grfprompt_a) gg->Write();
  for(auto& gg : grfprompt_b) gg->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==3) { drawlxydis(argv[1], argv[2]); return 0; }
  return 1;
}

void drawkinematic()
{
  xjjroot::drawtex(0.92, 0.84, Form("|y| < %s", xjjc::number_remove_zero(fitX::ycut).c_str()), 0.042, 32, 42);
  xjjroot::drawtex(0.92, 0.77, Form("%s < p_{T} < %s GeV/c", xjjc::number_remove_zero(fitX::ptmincut).c_str(), xjjc::number_remove_zero(fitX::ptmaxcut).c_str()), 0.042, 32, 42);
}
