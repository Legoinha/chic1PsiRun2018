#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "fit.h"
#include "MCefficiency.h"

void closure_fithist(std::string input, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  MCeff::MCefficiency mceff_a(inf, "_a");
  MCeff::MCefficiency mceff_b(inf, "_b");
  int nbins = mceff_a.nbins();
  std::vector<float> ptbins = mceff_a.ptbins();
  std::vector<TH1F*> h_a(nbins, 0), h_b(nbins, 0), hmcp_a(nbins, 0), hmcp_b(nbins, 0);
  std::vector<RooDataSet*> dsh_a(nbins, 0), dsh_b(nbins, 0), dshmcp_a(nbins, 0), dshmcp_b(nbins, 0);
  for(int i=0; i<mceff_a.nbins(); i++)
    {
      dsh_a[i]    = (RooDataSet*)ww->data(Form("dsh_a_%d", i));
      dsh_b[i]    = (RooDataSet*)ww->data(Form("dsh_b_%d", i));
      dshmcp_a[i] = (RooDataSet*)ww->data(Form("dshmcp_a_%d", i));
      dshmcp_b[i] = (RooDataSet*)ww->data(Form("dshmcp_b_%d", i));
      h_a[i]      = (TH1F*)inf->Get(Form("h_a_%d", i));
      h_b[i]      = (TH1F*)inf->Get(Form("h_b_%d", i));
      hmcp_a[i]   = (TH1F*)inf->Get(Form("hmcp_a_%d", i));
      hmcp_b[i]   = (TH1F*)inf->Get(Form("hmcp_b_%d", i));
      // hmcp_a[i]->Scale(hmcp_a[i]->GetEntries()/hmcp_a[i]->Integral());
      // hmcp_b[i]->Scale(hmcp_b[i]->GetEntries()/hmcp_b[i]->Integral());
    }

  // fit + yield
  TH1F* hyield_binned_a = (TH1F*)mceff_a.heffgen()->Clone("hyield_binned_a");
  TH1F* hyield_binned_b = (TH1F*)mceff_b.heffgen()->Clone("hyield_binned_b");
  TH1F* hyield_unbinned_a = (TH1F*)mceff_a.heffgen()->Clone("hyield_unbinned_a");
  TH1F* hyield_unbinned_b = (TH1F*)mceff_b.heffgen()->Clone("hyield_unbinned_b");

  for(int i=0; i<nbins; i++)
    {
      std::map<std::string, fitX::fitXresult*> result_a = fitX::fit(h_a[i], 0, hmcp_a[i], hmcp_b[i],
                                                                    dsh_a[i], dshmcp_a[i], dshmcp_b[i],
                                                                    Form("plots/%s/idx_a", output.c_str()), false, true, Form("_%d", i), 
                                                                    Form("%s < p_{T} < %s GeV/c", xjjc::number_remove_zero(ptbins[i]).c_str(), xjjc::number_remove_zero(ptbins[i+1]).c_str()));
      hyield_unbinned_a->SetBinContent(i+1, result_a["unbinned"]->ysig_a());
      hyield_unbinned_a->SetBinError(i+1, result_a["unbinned"]->ysigerr_a());
      hyield_binned_a->SetBinContent(i+1, result_a["binned"]->ysig_a());
      hyield_binned_a->SetBinError(i+1, result_a["binned"]->ysigerr_a());
      
      std::map<std::string, fitX::fitXresult*> result_b = fitX::fit(h_b[i], 0, hmcp_a[i], hmcp_b[i],
                                                                    dsh_b[i], dshmcp_a[i], dshmcp_b[i],
                                                                    Form("plots/%s/idx_b", output.c_str()), false, true, Form("_%d", i), 
                                                                    Form("%s < p_{T} < %s GeV/c", xjjc::number_remove_zero(ptbins[i]).c_str(), xjjc::number_remove_zero(ptbins[i+1]).c_str()));
      hyield_unbinned_b->SetBinContent(i+1, result_b["unbinned"]->ysig_b());
      hyield_unbinned_b->SetBinError(i+1, result_b["unbinned"]->ysigerr_b());
      hyield_binned_b->SetBinContent(i+1, result_b["binned"]->ysig_b());
      hyield_binned_b->SetBinError(i+1, result_b["binned"]->ysigerr_b());
    }

  // raw yields
  xjjroot::setthgrstyle(hyield_unbinned_a, xjjroot::mycolor_dark["red"], 20, 1.1, xjjroot::mycolor_dark["red"], 1, 1);
  xjjroot::setthgrstyle(hyield_binned_a, xjjroot::mycolor_light["red"], 24, 1.1, xjjroot::mycolor_light["red"], 1, 1);
  xjjroot::setthgrstyle(hyield_unbinned_b, xjjroot::mycolor_dark["red"], 20, 1.1, xjjroot::mycolor_dark["red"], 1, 1);
  xjjroot::setthgrstyle(hyield_binned_b, xjjroot::mycolor_light["red"], 24, 1.1, xjjroot::mycolor_light["red"], 1, 1);
  xjjroot::sethempty(hyield_unbinned_a);
  xjjroot::sethempty(hyield_binned_a);
  xjjroot::sethempty(hyield_unbinned_b);
  xjjroot::sethempty(hyield_binned_b);
  // efficiency
  mceff_a.calceff();
  mceff_a.setstyle(xjjroot::mycolor_middle["azure"]);
  mceff_b.calceff();
  mceff_b.setstyle(xjjroot::mycolor_middle["azure"]);
  // correct
  TH1F* hcorr_unbinned_a = (TH1F*)hyield_unbinned_a->Clone("hcorr_unbinned_a");
  hcorr_unbinned_a->Divide(mceff_a.heff());
  TH1F* hcorr_binned_a = (TH1F*)hyield_binned_a->Clone("hcorr_binned_a");
  hcorr_binned_a->Divide(mceff_a.heff());
  xjjroot::setthgrstyle(mceff_a.heffgen(), xjjroot::mycolor_middle["azure"], 21, 1, xjjroot::mycolor_middle["azure"], 1, 2, xjjroot::mycolor_middle["azure"], 0.2, 1001);
  xjjroot::sethempty(mceff_a.heffgen());
  TH1F* hcorr_unbinned_b = (TH1F*)hyield_unbinned_b->Clone("hcorr_unbinned_b");
  hcorr_unbinned_b->Divide(mceff_b.heff());
  TH1F* hcorr_binned_b = (TH1F*)hyield_binned_b->Clone("hcorr_binned_b");
  hcorr_binned_b->Divide(mceff_b.heff());
  xjjroot::setthgrstyle(mceff_b.heffgen(), xjjroot::mycolor_middle["azure"], 21, 1, xjjroot::mycolor_middle["azure"], 1, 2, xjjroot::mycolor_middle["azure"], 0.2, 1001);
  xjjroot::sethempty(mceff_b.heffgen());

  // draw
  float ymaxeff = 0.2;
  TH2F* hemptyeff = new TH2F("hemptyeff", ";p_{T} (GeV/c);#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 10, ptbins.front(), ptbins.back(), 10, 0, ymaxeff);
  xjjroot::sethempty(hemptyeff, 0, 0.3);
  xjjroot::setgstyle(1);
  TCanvas* c_a = new TCanvas("c_a", "", 1800, 600);
  c_a->SetLogy();
  c_a->Divide(3, 1);
  c_a->cd(1);
  hyield_binned_a->Draw("pe");
  hyield_unbinned_a->Draw("pe same");
  xjjroot::drawCMS("Simulation");
  c_a->cd(2);
  hemptyeff->Draw();
  mceff_a.greff()->Draw("same3");
  mceff_a.greff()->Draw("samelX");
  xjjroot::drawCMS("Simulation");
  c_a->cd(3);
  mceff_a.heffgen()->Draw("hist");
  hcorr_binned_a->Draw("pe same");
  hcorr_unbinned_a->Draw("pe same");
  xjjroot::drawCMS("Simulation");
  std::string outputname_a = "plots/"+output+"/cclosure_a.pdf";
  xjjroot::mkdir(outputname_a);
  c_a->SaveAs(outputname_a.c_str());
  TCanvas* c_b = new TCanvas("c_b", "", 1800, 600);
  c_b->SetLogy();
  c_b->Divide(3, 1);
  c_b->cd(1);
  hyield_binned_b->Draw("pe");
  hyield_unbinned_b->Draw("pe same");
  xjjroot::drawCMS("Simulation");
  c_b->cd(2);
  hemptyeff->Draw();
  mceff_b.greff()->Draw("same3");
  mceff_b.greff()->Draw("samelX");
  xjjroot::drawCMS("Simulation");
  c_b->cd(3);
  mceff_b.heffgen()->Draw("hist");
  hcorr_binned_b->Draw("pe same");
  hcorr_unbinned_b->Draw("pe same");
  xjjroot::drawCMS("Simulation");
  std::string outputname_b = "plots/"+output+"/cclosure_b.pdf";
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());

  // save
  std::string outputname = "rootfiles/"+output+"/closure_fit.root";
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  hemptyeff->Write();
  hyield_binned_a->Write();
  hyield_unbinned_a->Write();
  mceff_a.greff()->Write();
  mceff_a.heffgen()->Write();
  hcorr_binned_a->Write();
  hcorr_unbinned_a->Write();
  hyield_binned_b->Write();
  hyield_unbinned_b->Write();
  mceff_b.greff()->Write();
  mceff_b.heffgen()->Write();
  hcorr_binned_b->Write();
  hcorr_unbinned_b->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==3) { closure_fithist(argv[1], argv[2]); return 0; }
  return 1;
}
