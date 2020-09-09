#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooMinuit.h>
#include <RooAddPdf.h>
#include <TH1F.h>
#include <string>

namespace fitX
{
  class xnll
  {
  public:
    xnll(RooWorkspace* ww, std::string name);
    ~xnll() {;}
    int nbin() { return fnbin; }
    float minbin() { return fminbin; }
    float maxbin() { return fmaxbin; }
    RooAbsReal* nll() { return fnll; }
    RooAbsReal* pll() { return fpll; }
    RooAddPdf* pdf() { return fpdf; }
    RooPolynomial* bkg() { return fbkg; }
    RooAbsReal* ndof() { return fndof; }
    TH1F* hnll() { return fhnll; }
    TH1F* hnllcorr1() { return fhnllcorr1; }
    TH1F* hnllcorr2() { return fhnllcorr2; }
    TH1F* hpll() { return fhpll; }
    double min_nll() { return fmin_nll; }
    double max_nll() { return fmax_nll; }
    double min_nllcorr1() { return fmin_nllcorr1; }
    double max_nllcorr1() { return fmax_nllcorr1; }
    double min_nllcorr2() { return fmin_nllcorr2; }
    double max_nllcorr2() { return fmax_nllcorr2; }
    double zero_nll() { return fzero_nll; } 
    double xmean_nll() { return fxmean_nll; }
    double min_pll() { return fmin_pll; }
    double max_pll() { return fmax_pll; }
    double zero_pll() { return fzero_pll; } 
    double xmean_pll() { return fxmean_pll; }
    double x68_d() { return fx68_d; }
    double x68_u() { return fx68_u; }
    double x95_d() { return fx95_d; }
    double x95_u() { return fx95_u; }
    double sig() { return fsig; }
    double nullpvalue() { return fnullpvalue->getVal(); }
    double nullsig() { return fnullsig->getVal(); }
    void createHist();

  private:
    std::string fname;
    RooAbsReal* fnll;
    RooAbsReal* fpll;
    RooAddPdf* fpdf;
    RooPolynomial* fbkg;
    RooAbsReal* fndof;
    TH1F* fhnll;
    TH1F* fhnllcorr1;
    TH1F* fhnllcorr2;
    TH1F* fhpll;
    RooAbsReal* fnullpvalue;
    RooAbsReal* fnullsig;

    double fmin_nll;
    double fmax_nll;
    double fmin_nllcorr1;
    double fmax_nllcorr1;
    double fmin_nllcorr2;
    double fmax_nllcorr2;
    double fzero_nll; 
    double fxmean_nll;
    double fmin_pll;
    double fmax_pll;
    double fzero_pll; 
    double fxmean_pll;
    double fx68_d;
    double fx68_u;
    double fx95_d;
    double fx95_u;
    double fsig;

    RooRealVar* fpar10;
    RooAbsData* fdata;

    const int fnbin = 400;
    const float fbinwid = 0.5;
    // const float fminbin = 0-fbinwid/2., fmaxbin = fminbin + fnbin*fbinwid;
    const float fminbin = 0-fbinwid/2.;
    const float fmaxbin = fminbin + fnbin*fbinwid;

  };
}

fitX::xnll::xnll(RooWorkspace* ww, std::string name) : fname(name)
{
  fpar10 = ww->var("par10");
  fndof = ww->var("ndof");
  fpdf = (RooAddPdf*)ww->obj("pdf");
  fbkg = (RooPolynomial*)ww->obj("bkg");
  fdata = ww->data("dsh_ws");
  fnll = fpdf->createNLL(*fdata);
  // RooMinuit(*fnll).migrad();
  fpll = fnll->createProfile(*fpar10);
  fnullpvalue = ww->var("nullpvalue");
  fnullsig = ww->var("nullsig");
}

void fitX::xnll::createHist()
{
  fhnll = (TH1F*)fnll->createHistogram(Form("hnll_%s", fname.c_str()), *fpar10, RooFit::Binning(fnbin, fminbin, fmaxbin), RooFit::Scaling(kFALSE));
  for(int j=0; j<fnbin; j++) { fhnll->SetBinContent(j+1, fhnll->GetBinContent(j+1)*2.); }
  fhnll->GetYaxis()->SetTitle("Projection of -2log(LL)");
  fhpll = (TH1F*)fpll->createHistogram(Form("hpll_%s", fname.c_str()), *fpar10, RooFit::Binning(fnbin, fminbin, fmaxbin), RooFit::Scaling(kFALSE));
  for(int j=0; j<fnbin; j++) { fhpll->SetBinContent(j+1, fhpll->GetBinContent(j+1)*2.); }
  fhpll->GetYaxis()->SetTitle("Projection of -2log(profileLL)");

  fhnllcorr1 = (TH1F*)fhnll->Clone(Form("hnllcorr1_%s", fname.c_str()));
  for(int j=0; j<fnbin; j++) { fhnllcorr1->SetBinContent(j+1, fhnll->GetBinContent(j+1)+fndof->getVal()); }
  fhnllcorr1->GetYaxis()->SetTitle(Form("%s + N_{par}", fhnll->GetYaxis()->GetTitle()));
  fhnllcorr2 = (TH1F*)fhnll->Clone(Form("hnllcorr2_%s", fname.c_str()));
  fhnllcorr2->GetYaxis()->SetTitle(Form("%s + 2N_{par}", fhnll->GetYaxis()->GetTitle()));
  for(int j=0; j<fnbin; j++) { fhnllcorr2->SetBinContent(j+1, fhnll->GetBinContent(j+1)+fndof->getVal()*2.); }

  fmin_nll = fhnll->GetBinContent(fhnll->GetMinimumBin());
  fmax_nll = fhnll->GetBinContent(fhnll->GetMaximumBin());
  fzero_nll = fhnll->GetBinContent(1);
  fxmean_nll = fnll->findRoot(*fpar10, fminbin, fmaxbin, fmin_nll/2.);

  fmin_pll = fhpll->GetBinContent(fhpll->GetMinimumBin());
  fmax_pll = fhpll->GetBinContent(fhpll->GetMaximumBin());
  fzero_pll = fhpll->GetBinContent(1);
  fxmean_pll = fhpll->GetBinCenter(fhpll->GetMinimumBin());

  fx68_d = fpll->findRoot(*fpar10, fminbin, fmin_pll, fmin_pll+0.50);
  fx68_u = fpll->findRoot(*fpar10, fmin_pll, fmaxbin, fmin_pll+0.50);
  fx95_d = fpll->findRoot(*fpar10, fminbin, fmin_pll, fmin_pll+1.96);
  fx95_u = fpll->findRoot(*fpar10, fmin_pll, fmaxbin, fmin_pll+1.96);

  fsig = sqrt(fzero_pll-fmin_pll);

  fmin_nllcorr1 = fhnllcorr1->GetBinContent(fhnllcorr1->GetMinimumBin());
  fmax_nllcorr1 = fhnllcorr1->GetBinContent(fhnllcorr1->GetMaximumBin());
  fmin_nllcorr2 = fhnllcorr2->GetBinContent(fhnllcorr2->GetMinimumBin());
  fmax_nllcorr2 = fhnllcorr2->GetBinContent(fhnllcorr2->GetMaximumBin());
}
