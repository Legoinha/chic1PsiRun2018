#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "lxydis.h"
#include "fit.h"

void fprompt_fithist(std::string input, std::string output, std::string lxyvar="lxy")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  output += ("/"+lxyvar);
  int ncut = lxydis::lxycut[lxyvar].size();
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  RooDataSet* dsh = (RooDataSet*)ww->data("dsh");
  std::vector<RooDataSet*> dshBenr(ncut);
  for(int i=0; i<ncut; i++)
    { dshBenr[i] = (RooDataSet*)ww->data(Form("dshBenr-%d",i)); }
  RooDataSet* dshmcp_a = (RooDataSet*)ww->data("dshmcp_a");
  RooDataSet* dshmcp_b = (RooDataSet*)ww->data("dshmcp_b");
  TH1F* h = (TH1F*)inf->Get("h");
  std::vector<TH1F*> hBenr(ncut);
  for(int i=0; i<ncut; i++)
    { hBenr[i] = (TH1F*)inf->Get(Form("hBenr-%d",i)); }
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
  TH1F* hlxymcnp_a = (TH1F*)inf->Get("hlxymcnp_a");
  TH1F* hlxymcnp_b = (TH1F*)inf->Get("hlxymcnp_b");
  TH1F* hlxymcp_a = (TH1F*)inf->Get("hlxymcp_a");
  TH1F* hlxymcp_b = (TH1F*)inf->Get("hlxymcp_b");

  // fit + yield
  TH1F* hyield_a = new TH1F("hyield_a", "", 5, 0, 5); //hyield_a->Sumw2();
  TH1F* hyield_b = new TH1F("hyield_b", "", 5, 0, 5); //hyield_b->Sumw2();
  std::vector<TH1F*> hBenryield_a(ncut), hBenryield_b(ncut);
  for(int i=0; i<ncut; i++)
    {
      hBenryield_a[i] = new TH1F(Form("hBenryield_a-%d", i), "", 5, 0, 5); //hBenryield_a->Sumw2();
      hBenryield_b[i] = new TH1F(Form("hBenryield_b-%d", i), "", 5, 0, 5); //hBenryield_b->Sumw2();
    }

  std::map<std::string, fitX::fitXresult*> resulth = fitX::fit(h, 0, hmcp_a, hmcp_b, 
                                                               dsh, dshmcp_a, dshmcp_b,
                                                               Form("plots/%s/didx", output.c_str()), 0, true, "th", 
                                                               "Inclusive");
  hyield_a->SetBinContent(fitX::ibin_a, resulth["unbinned"]->ysig_a());
  hyield_a->SetBinError(fitX::ibin_a, resulth["unbinned"]->ysigerr_a());
  hyield_b->SetBinContent(fitX::ibin_b, resulth["unbinned"]->ysig_b());
  hyield_b->SetBinError(fitX::ibin_b, resulth["unbinned"]->ysigerr_b());
  xjjroot::setthgrstyle(hyield_a, kGray+3, 20, 1.2, kGray+3, 1, 2, kGray+3, 0.1, 1001);
  xjjroot::setthgrstyle(hyield_b, kGray+3, 20, 1.2, kGray+3, 1, 2, kGray+3, 0.1, 1001);
  TLegend* leg = new TLegend(0.60, 0.65-0.045*ncut, 0.90, 0.65);
  xjjroot::setleg(leg, 0.04);
  float ymaxh = std::max(hyield_a->GetMaximum(), hyield_b->GetMaximum()), ymaxhBenr = 0;
  for(int i=0; i<ncut; i++)
    {
      std::map<std::string, fitX::fitXresult*> resulthBenr = fitX::fit(hBenr[i], 0, hmcp_a, hmcp_b, 
                                                                       dshBenr[i], dshmcp_a, dshmcp_b,
                                                                       Form("plots/%s/didx", output.c_str()), resulth["unbinned"]->msig_b(), true, Form("thBenr-%d",i), 
                                                                       Form("B-enr (%s > %s)", lxydis::vars[lxyvar].c_str(), xjjc::number_remove_zero(lxydis::lxycut[lxyvar][i]).c_str()));
      hBenryield_a[i]->SetBinContent(fitX::ibin_a, resulthBenr["unbinned"]->ysig_a());
      hBenryield_a[i]->SetBinError(fitX::ibin_a, resulthBenr["unbinned"]->ysigerr_a());
      hBenryield_b[i]->SetBinContent(fitX::ibin_b, resulthBenr["unbinned"]->ysig_b());
      hBenryield_b[i]->SetBinError(fitX::ibin_b, resulthBenr["unbinned"]->ysigerr_b());
      float thisymaxhBenr = std::max(resulthBenr["unbinned"]->ysig_a(), resulthBenr["unbinned"]->ysig_b());
      xjjroot::setthgrstyle(hBenryield_a[i], xjjroot::mycolor_middle[xjjroot::cc[i]], 24, 1.2, xjjroot::mycolor_middle[xjjroot::cc[i]], 1, 2, xjjroot::mycolor_middle[xjjroot::cc[i]], 0.1, 1001);
      xjjroot::setthgrstyle(hBenryield_b[i], xjjroot::mycolor_middle[xjjroot::cc[i]], 24, 1.2, xjjroot::mycolor_middle[xjjroot::cc[i]], 1, 2, xjjroot::mycolor_middle[xjjroot::cc[i]], 0.1, 1001);
      ymaxhBenr = std::max(ymaxhBenr, thisymaxhBenr);
      leg->AddEntry(hBenryield_a[i], Form("%s > %s", lxydis::vars[lxyvar].c_str(), xjjc::number_remove_zero(lxydis::lxycut[lxyvar][i]).c_str()), "pl");
    }

  TH2F* hemptyh = new TH2F("hemptyh", ";;Raw Yield", 5, 0, 5, 10, 0, ymaxh*2.5);
  xjjroot::sethempty(hemptyh, 0, 0.3);
  hemptyh->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hemptyh->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hemptyh->GetXaxis()->SetLabelSize(hemptyh->GetXaxis()->GetLabelSize()*1.5);
  TH2F* hemptyhBenr = new TH2F("hemptyhBenr", ";;Raw Yield", 5, 0, 5, 10, 0, ymaxhBenr*2.5);
  xjjroot::sethempty(hemptyhBenr, 0, 0.3);
  hemptyhBenr->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hemptyhBenr->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hemptyhBenr->GetXaxis()->SetLabelSize(hemptyhBenr->GetXaxis()->GetLabelSize()*1.5);

  // draw raw yield
  xjjroot::setgstyle(1);
  TCanvas* cy = new TCanvas("cy", "", 1200, 600);
  cy->Divide(2, 1);
  cy->cd(1);
  hemptyh->Draw();
  hyield_a->Draw("pe same");
  hyield_b->Draw("pe same");
  xjjroot::drawtex(0.22, 0.84, "Inclusive", 0.042, 12, 62);
  xjjroot::drawCMS();
  fitX::drawkinematics();
  cy->cd(2);
  hemptyhBenr->Draw();
  for(int i=0; i<ncut; i++) 
    {
      hBenryield_a[i]->Draw("pe same");
      hBenryield_b[i]->Draw("pe same");
    }
  leg->Draw();
  xjjroot::drawtex(0.22, 0.84, "B-enriched", 0.042, 12, 62);
  xjjroot::drawCMS();
  fitX::drawkinematics();
  xjjroot::mkdir(Form("plots/%s/cyield.pdf", output.c_str()));
  cy->SaveAs(Form("plots/%s/cyield.pdf", output.c_str()));

  // fprompt
  hlxymcnp_a->Scale(1./hlxymcnp_a->Integral(), "width");
  hlxymcnp_b->Scale(1./hlxymcnp_b->Integral(), "width");
  hlxymcp_a->Scale(1./hlxymcp_a->Integral(), "width");
  hlxymcp_b->Scale(1./hlxymcp_b->Integral(), "width");
  xjjroot::sethempty(hlxymcnp_a, 0, 0);
  xjjroot::sethempty(hlxymcp_a, 0, 0);
  xjjroot::setthgrstyle(hlxymcp_a, fitX::color_a, 21, 1, fitX::color_a, 1, 2);
  xjjroot::setthgrstyle(hlxymcp_b, fitX::color_b, 21, 1, fitX::color_b, 1, 2);
  xjjroot::setthgrstyle(hlxymcnp_a, fitX::color_a, 21, 1, fitX::color_a, 1, 2);
  xjjroot::setthgrstyle(hlxymcnp_b, fitX::color_b, 21, 1, fitX::color_b, 1, 2);
  hlxymcp_a->SetMaximum(hlxymcp_a->GetMaximum()*10.);
  hlxymcnp_a->SetMaximum(hlxymcnp_a->GetMaximum()*10.);
  TCanvas* clxy = new TCanvas("clxy", "", 1200, 600);
  clxy->Divide(2, 1);
  clxy->cd(1);
  gPad->SetLogy();
  hlxymcp_a->Draw("histe");
  hlxymcp_b->Draw("histe same");
  for(int i=0; i<ncut; i++) 
    xjjroot::drawline(lxydis::lxycut[lxyvar][i], 0, lxydis::lxycut[lxyvar][i], hlxymcp_a->GetMaximum(), xjjroot::mycolor_middle[xjjroot::cc[i]], 2, 2, 0.6);
  fitX::drawkinematics();
  xjjroot::drawtex(0.23, 0.85, "#psi(2S)", 0.042, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.23, 0.78, "X(3872)", 0.042, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.57, 0.25, "Prompt", 0.042, 22, 62);
  xjjroot::drawCMS("Simulation");
  clxy->cd(2);
  gPad->SetLogy();
  hlxymcnp_a->Draw("histe");
  hlxymcnp_b->Draw("histe same");
  for(int i=0; i<ncut; i++) 
    xjjroot::drawline(lxydis::lxycut[lxyvar][i], 0, lxydis::lxycut[lxyvar][i], hlxymcnp_a->GetMaximum(), xjjroot::mycolor_middle[xjjroot::cc[i]], 2, 2, 0.6);
  fitX::drawkinematics();
  xjjroot::drawtex(0.23, 0.85, "#psi(2S)", 0.042, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.23, 0.78, "X(3872)", 0.042, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.57, 0.25, "Nonprompt", 0.042, 22, 62);
  xjjroot::drawCMS("Simulation");
  xjjroot::mkdir(Form("plots/%s/clxy.pdf", output.c_str()));
  clxy->SaveAs(Form("plots/%s/clxy.pdf", output.c_str()));

  TLegend* legfprompt = new TLegend(0.60, 0.50-0.045*ncut, 0.90, 0.50);
  xjjroot::setleg(legfprompt, 0.04);
  std::vector<TEfficiency*> grfprompt_a(ncut), grfprompt_b(ncut); 
  for(int i=0; i<ncut; i++)
    {
      std::vector<double> vlxyfrac = lxydis::nplxyfrac(hlxymcnp_a, hlxymcnp_b, lxydis::lxycut[lxyvar][i]);
      TH1F* hlxyfrac_a = new TH1F("hlxyfrac_a", "", 5, 0, 5); //hlxyfrac_a->Sumw2();
      hlxyfrac_a->SetBinContent(fitX::ibin_a, vlxyfrac[0]);
      hlxyfrac_a->SetBinError(fitX::ibin_a, vlxyfrac[1]);
      TH1F* hlxyfrac_b = new TH1F("hlxyfrac_b", "", 5, 0, 5); //hlxyfrac_b->Sumw2();
      hlxyfrac_b->SetBinContent(fitX::ibin_b, vlxyfrac[2]);
      hlxyfrac_b->SetBinError(fitX::ibin_b, vlxyfrac[3]);
      TH1F *hyieldprompt_a, *hyieldprompt_b;
      grfprompt_a[i] = lxydis::calclxyfprompt(hyield_a, hBenryield_a[i], hlxyfrac_a, Form("grfprompt_a-%d", i), &hyieldprompt_a);
      grfprompt_b[i] = lxydis::calclxyfprompt(hyield_b, hBenryield_b[i], hlxyfrac_b, Form("grfprompt_b-%d", i), &hyieldprompt_b);
      xjjroot::setthgrstyle(grfprompt_a[i], xjjroot::mycolor_middle[xjjroot::cc[i]], 25, 1.5, xjjroot::mycolor_middle[xjjroot::cc[i]], 1, 3);
      xjjroot::setthgrstyle(grfprompt_b[i], xjjroot::mycolor_middle[xjjroot::cc[i]], 25, 1.5, xjjroot::mycolor_middle[xjjroot::cc[i]], 1, 3);
      legfprompt->AddEntry(grfprompt_a[i], Form("%s > %s", lxydis::vars[lxyvar].c_str(), xjjc::number_remove_zero(lxydis::lxycut[lxyvar][i]).c_str()), "pl");
      delete hlxyfrac_a;
      delete hlxyfrac_b;
      delete hyieldprompt_a;
      delete hyieldprompt_b;
    }
  xjjroot::setgstyle(2);
  TCanvas* cfprompt = new TCanvas("cfprompt", "", 600, 600);
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";;f_{prompt}", 5, 0, 5, 10, 0, 1.4);
  xjjroot::sethempty(hemptyfprompt, 0, 0.3);
  hemptyfprompt->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hemptyfprompt->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hemptyfprompt->GetXaxis()->SetLabelSize(hemptyfprompt->GetXaxis()->GetLabelSize()*1.5);
  hemptyfprompt->Draw();
  xjjroot::drawline(0, 1, 5, 1, kGray+1, 9, 2);
  for(int i=0; i<ncut; i++)
    {
      grfprompt_a[i]->Draw("pe same");
      grfprompt_b[i]->Draw("pe same");
    }
  legfprompt->Draw();
  fitX::drawkinematics();
  xjjroot::drawcomment(output.c_str(), "r");
  xjjroot::drawCMS();
  xjjroot::mkdir(Form("plots/%s/cfprompt.pdf", output.c_str()));
  cfprompt->SaveAs(Form("plots/%s/cfprompt.pdf", output.c_str()));

  // write
  std::string outputname(Form("rootfiles/%s/fprompt_fithist.root", output.c_str()));
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hyield_a->Write();
  hyield_b->Write();
  for(auto& hh : hBenryield_a) hh->Write();
  for(auto& hh : hBenryield_b) hh->Write();
  for(auto& gg : grfprompt_a) gg->Write();
  for(auto& gg : grfprompt_b) gg->Write();
  fitX::write();
  outf->Close();
  std::cout<<"output: "<<outputname<<std::endl;
  std::cout<<std::endl;

  float ix, iy;
  iy = grfprompt_a[0]->GetEfficiency(fitX::ibin_a);
  float nominal_a = iy; std::cout<<iy<<std::endl;
  float per_a = 0;
  for(auto& gg : grfprompt_a)
    {
      iy = gg->GetEfficiency(fitX::ibin_a);
      if(TMath::Abs(iy - nominal_a)/nominal_a > per_a) { per_a = TMath::Abs(iy - nominal_a)/nominal_a; }
    }
  iy = grfprompt_b[0]->GetEfficiency(fitX::ibin_b);
  float nominal_b = iy;
  float per_b = 0;
  for(auto& gg : grfprompt_b)
    {
      iy = gg->GetEfficiency(fitX::ibin_b);
      if(TMath::Abs(iy - nominal_b)/nominal_b > per_b) { per_b = TMath::Abs(iy - nominal_b)/nominal_b; }
    }
  float per_ab = TMath::Sqrt(per_a*per_a + per_b*per_b);
  std::cout<<"Prompt Fraction & "<<Form("%.1f", per_a*1.e+2)<<"\\% & "<<Form("%.1f", per_b*1.e+2)<<"\\% & "<<Form("%.1f", per_ab*1.e+2)<<"\\% \\\\"<<std::endl;

}

int main(int argc, char* argv[])
{
  fitX::init(TFile::Open(argv[1]));
  std::string outputname = std::string(argv[2])+fitX::tagname();
  if(argc==4) { fprompt_fithist(argv[1], outputname, argv[3]); return 0; }
  return 1;
}
