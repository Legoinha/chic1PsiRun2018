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
  TH1F* hcorr_unbinned_b = (TH1F*)hyield_unbinned_b->Clone("hcorr_unbinned_b");
  hcorr_unbinned_b->Divide(mceff_b.heff());
  TH1F* hcorr_binned_b = (TH1F*)hyield_binned_b->Clone("hcorr_binned_b");
  hcorr_binned_b->Divide(mceff_b.heff());

  // save
  std::string outputname = "rootfiles/"+output+"/closure_fithist.root";
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  hyield_binned_a->Write();
  hyield_unbinned_a->Write();
  mceff_a.heff()->Write();
  mceff_a.greff()->Write();
  mceff_a.heffgen()->Write();
  hcorr_binned_a->Write();
  hcorr_unbinned_a->Write();
  hyield_binned_b->Write();
  hyield_unbinned_b->Write();
  mceff_b.heff()->Write();
  mceff_b.greff()->Write();
  mceff_b.heffgen()->Write();
  hcorr_binned_b->Write();
  hcorr_unbinned_b->Write();
  fitX::write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==3) { closure_fithist(argv[1], argv[2]); return 0; }
  return 1;
}

