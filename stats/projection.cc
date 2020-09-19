#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooMinuit.h>

#include <string>

#include "xjjrootuti.h"
#include "fit.h"
// #include "smear.h"

void projection(std::string input, std::string output, float lumiscale)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  RooDataSet* dsh = (RooDataSet*)ww->data("dsh");
  RooDataSet* dshBenr = (RooDataSet*)ww->data("dshBenr");
  RooDataSet* dshmcp_a = (RooDataSet*)ww->data("dshmcp_a");
  RooDataSet* dshmcp_b = (RooDataSet*)ww->data("dshmcp_b");
  TH1F* h = (TH1F*)inf->Get("h");
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  hmcp_a->Scale(hmcp_a->GetEntries()/hmcp_a->Integral());
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
  hmcp_b->Scale(hmcp_b->GetEntries()/hmcp_b->Integral());

  std::map<std::string, fitX::fitXresult*> result = fitX::fit(h, 0, hmcp_a, hmcp_b,
                                                              dsh, dshmcp_a, dshmcp_b,
                                                              Form("plots/%s/idx", output.c_str()), 0, true, "nominal", "2018 data", "default", true); // fix mean = false
  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooWorkspace* wo = new RooWorkspace("wo");
  TF1* f = result["unbinned"]->f();
  int nentries = h->GetEntries();
  TH1F* hh = new TH1F("hh", ";;", fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
  RooDataSet* dshh = new RooDataSet("dshh", "", RooArgSet(*mass));
  for(int j=0; j<nentries*lumiscale; j++)
    {
      float mm = f->GetRandom(fitX::BIN_MIN, fitX::BIN_MAX);
      mass->setVal(mm);
      dshh->add(*mass);
      hh->Fill(mm, 1);
    }
  fitX::printhist(hh);
  dshh->Print();
  wo->import(*dshh);
  wo->import(*dshmcp_a);
  wo->import(*dshmcp_b);
  std::map<std::string, fitX::fitXresult*> rt = fitX::fit(hh, 0, hmcp_a, hmcp_b,
                                                          dshh, dshmcp_a, dshmcp_b,
                                                          Form("plots/%s/idx", output.c_str()), 0, true, "projection", "Projection", "default", true,
                                                          "#scale[1.25]{#bf{CMS}}", "10 nb^{-1} (PbPb 5.02 TeV)"); // fix mean = false
  // float ysig_a = rt["unbinned"]->ysig_a();
  // float ysigerr_a = rt["unbinned"]->ysigerr_a();
  // float ysig_b = rt["unbinned"]->ysig_b();
  // float ysigerr_b = rt["unbinned"]->ysigerr_b();
  RooWorkspace* wt = rt["unbinned"]->ww();
  wt->Print();
  RooRealVar* par10 = wt->var("par10");
  RooAddPdf* pdf = (RooAddPdf*)wt->obj("pdf");
  RooAbsData* data = wt->data("dshh_ws");

  std::string outputname = Form("rootfiles/%s/projection.root", output.c_str());
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();

  xjjroot::setgstyle();
  TCanvas* cnll = new TCanvas("cnll", "", 600, 600);
  RooAbsReal* nll = pdf->createNLL(*data);
  RooMinuit(*nll).migrad();
  RooPlot* frame1 = par10->frame(RooFit::Bins(20), RooFit::Range(0., 1200.));
  xjjroot::sethempty(frame1);
  nll->plotOn(frame1, RooFit::ShiftToZero(), RooFit::LineColor(kBlack), RooFit::LineWidth(4), RooFit::Normalization(2., RooAbsReal::Raw));
  cnll->cd();
  frame1->SetMinimum(0);
  frame1->SetMaximum(140);
  frame1->GetXaxis()->SetTitle("N_{sig}^{(X3872)}");
  frame1->GetYaxis()->SetTitle("-2log#lambda");
  frame1->Draw();
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Projection}", 0.05, -0.08);
  xjjroot::drawCMSright("10 nb^{-1} (PbPb 5.02 TeV)");
  std::string outputnamenll = "plots/" + output + "/cnllscan.pdf";
  xjjroot::mkdir(outputnamenll);
  cnll->SaveAs(outputnamenll.c_str());

  gDirectory->cd("root:/");
  RooWorkspace* wwnll = new RooWorkspace("wwnll");
  wwnll->import(*pdf);
  wwnll->import(*data);
  wwnll->import(*par10);
  outf->cd();
  gDirectory->Add(wwnll);
  wwnll->Write();
  wwnll->Print();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==4) { projection(argv[1], argv[2], atof(argv[3])); return 0; }
  return 1;
}

