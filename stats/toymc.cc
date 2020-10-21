#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "xjjrootuti.h"
#include "fit.h"

const int nn = 1000; // 600

void toymc(std::string input, std::string output)
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

  std::string outputname = Form("rootfiles/%s/toymc.root", output.c_str());
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  TTree* nt = new TTree("nt", "");
  float ysig_a; nt->Branch("ysig_a", &ysig_a, "ysig_a/F");
  float ysigerr_a; nt->Branch("ysigerr_a", &ysigerr_a, "ysigerr_a/F");
  float ysig_b; nt->Branch("ysig_b", &ysig_b, "ysig_b/F");
  float ysigerr_b; nt->Branch("ysigerr_b", &ysigerr_b, "ysigerr_b/F");
  float msig_a; nt->Branch("msig_a", &msig_a, "msig_a/F");
  float msigerr_a; nt->Branch("msigerr_a", &msigerr_a, "msigerr_a/F");
  float msig_b; nt->Branch("msig_b", &msig_b, "msig_b/F");
  float msigerr_b; nt->Branch("msigerr_b", &msigerr_b, "msigerr_b/F");
  float minNll; nt->Branch("minNll", &minNll, "minNll/F");

  std::map<std::string, fitX::fitXresult*> result = fitX::fit(h, 0, hmcp_a, hmcp_b,
                                                              dsh, dshmcp_a, dshmcp_b,
                                                              Form("plots/%s/idx", output.c_str()), 0, true, "nominal", "default", "real-data", true); // fix mean = false
  ysig_a = result["unbinned"]->ysig_a();
  ysigerr_a = result["unbinned"]->ysigerr_a();
  ysig_b = result["unbinned"]->ysig_b();
  ysigerr_b = result["unbinned"]->ysigerr_b();
  msig_a = result["unbinned"]->msig_a();
  msigerr_a = result["unbinned"]->msigerr_a();
  msig_b = result["unbinned"]->msig_b();
  msigerr_b = result["unbinned"]->msigerr_b();
  minNll = result["unbinned"]->get("minNll");
  nt->Fill();

  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooWorkspace* wo = new RooWorkspace("wo");
  TF1* f = result["unbinned"]->f();
  int nentries = h->GetEntries();
  std::vector<TH1F*> hh(nn);
  std::vector<RooDataSet*> dshh(nn);
  for(int i=0; i<nn; i++)
    {
      dshh[i] = new RooDataSet(Form("dsh_%d", i), "", RooArgSet(*mass));
      // hh[i] = smear::smear(h, Form("_%d", i));
      hh[i] = new TH1F(Form("h_%d", i), ";;", fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX);
      for(int j=0; j<nentries; j++)
        {
          float mm = f->GetRandom(fitX::BIN_MIN, fitX::BIN_MAX);
          mass->setVal(mm);
          dshh[i]->add(*mass);
          hh[i]->Fill(mm, 1);
        }
      fitX::printhist(hh[i]);
      dshh[i]->Print();
      wo->import(*dshh[i]);
      std::map<std::string, fitX::fitXresult*> rt = fitX::fit(hh[i], 0, hmcp_a, hmcp_b,
                                                              dshh[i], dshmcp_a, dshmcp_b,
                                                              Form("plots/%s/idx", output.c_str()), 0, i%100==0, Form("-%d",i), "default", Form("pseudo-data (%d)", i), true); // fix mean = false
      ysig_a = rt["unbinned"]->ysig_a();
      ysigerr_a = rt["unbinned"]->ysigerr_a();
      ysig_b = rt["unbinned"]->ysig_b();
      ysigerr_b = rt["unbinned"]->ysigerr_b();
      msig_a = rt["unbinned"]->msig_a();
      msigerr_a = rt["unbinned"]->msigerr_a();
      msig_b = rt["unbinned"]->msig_b();
      msigerr_b = rt["unbinned"]->msigerr_b();
      minNll = rt["unbinned"]->get("minNll");
      nt->Fill();
    }
  outf->cd();
  nt->Write();

  // for(auto& ih : hh) { ih->Write(); }
  // outf->cd();
  // gDirectory->Add(wo);
  // wo->Write();
  // wo->Print();
  // outf->cd();
  // hysig_a->Write();
  // hysig_b->Write();
  // hmsig_a->Write();
  // hmsig_b->Write();
  // hminNll->Write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==3) { toymc(argv[1], argv[2]); return 0; }
  return 1;
}

