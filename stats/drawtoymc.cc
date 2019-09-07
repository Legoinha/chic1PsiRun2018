#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>

#include <string>

#include "xjjrootuti.h"
#include "fit.h"

int nn = 600;
void labelsmc(std::string label, double mean, double sigma, std::string var, int precision, std::string unit);
void drawtoymc(std::string input, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());

  // --> ysig

  TH1F* hysig_a = (TH1F*)inf->Get("hysig_a");
  hysig_a->GetXaxis()->SetNdivisions(505);
  hysig_a->GetXaxis()->SetTitle("N_{Signal}");
  hysig_a->SetMaximum(hysig_a->GetMaximum()*1.2);
  xjjroot::setthgrstyle(hysig_a, fitX::color_a, 20, 1., fitX::color_a, 1, 1);
  xjjroot::sethempty(hysig_a);
  TH1F* hysig_b = (TH1F*)inf->Get("hysig_b");
  hysig_b->GetXaxis()->SetNdivisions(505);
  hysig_b->GetXaxis()->SetTitle("N_{Signal}");
  hysig_b->SetMaximum(hysig_b->GetMaximum()*1.2);
  xjjroot::setthgrstyle(hysig_b, fitX::color_b, 20, 1., fitX::color_b, 1, 1);
  xjjroot::sethempty(hysig_b);

  TF1* fysig_a = new TF1("fysig_a", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fysig_a, fitX::color_a, 7, 4);
  fysig_a->SetParLimits(0, 0, 1.e+4);
  fysig_a->SetParameter(0, nn);
  fysig_a->SetParLimits(1, 0, 200);
  fysig_a->SetParameter(1, 100);
  fysig_a->SetParameter(2, 30);
  fysig_a->SetParLimits(2, 0, 60);
  TF1* fysig_b = new TF1("fysig_b", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fysig_b, fitX::color_b, 7, 4);
  fysig_b->SetParLimits(0, 0, 1.e+4);
  fysig_b->SetParameter(0, nn);
  fysig_b->SetParLimits(1, 0, 200);
  fysig_b->SetParameter(1, 100);
  fysig_b->SetParameter(2, 30);
  fysig_b->SetParLimits(2, 0, 60);

  xjjroot::setgstyle();
  TCanvas* cysig = new TCanvas("cysig", "", 1200, 600);
  cysig->Divide(2, 1);
  cysig->cd(1);
  hysig_a->Draw("pe");
  hysig_a->Fit("fysig_a", "Nq");
  hysig_a->Fit("fysig_a", "Nq");
  hysig_a->Fit("fysig_a", "LL");
  fitX::drawpull(hysig_a, fysig_a, fitX::color_a);
  labelsmc("#psi(2S)", fysig_a->GetParameter(1), fysig_a->GetParameter(2), "N", 2, "");
  xjjroot::drawcomment(output.c_str());
  cysig->cd(2);
  hysig_b->Draw("pe");
  hysig_b->Fit("fysig_b", "Nq");
  hysig_b->Fit("fysig_b", "Nq");
  hysig_b->Fit("fysig_b", "LL");
  fitX::drawpull(hysig_b, fysig_b, fitX::color_b);
  labelsmc("X(3872)", fysig_b->GetParameter(1), fysig_b->GetParameter(2), "N", 2, "");
  std::string outputname_ysig = Form("plots/%s/ctoymc_ysig.pdf", output.c_str());
  xjjroot::mkdir(outputname_ysig);
  cysig->SaveAs(outputname_ysig.c_str());
  std::cout<<fysig_a->GetParameter(2)<<" "<<fysig_b->GetParameter(2)<<std::endl;;

  // --> msig

  TH1F* hmsig_a = (TH1F*)inf->Get("hmsig_a");
  hmsig_a->GetXaxis()->SetNdivisions(505);
  hmsig_a->GetXaxis()->SetTitle("#bar{m} (GeV/c^{2})");
  hmsig_a->SetMaximum(hmsig_a->GetMaximum()*1.2);
  xjjroot::setthgrstyle(hmsig_a, fitX::color_a, 20, 1., fitX::color_a, 1, 1);
  xjjroot::sethempty(hmsig_a);
  TH1F* hmsig_b = (TH1F*)inf->Get("hmsig_b");
  hmsig_b->GetXaxis()->SetNdivisions(505);
  hmsig_b->GetXaxis()->SetTitle("#bar{m} (GeV/c^{2})");
  hmsig_b->SetMaximum(hmsig_b->GetMaximum()*1.2);
  xjjroot::setthgrstyle(hmsig_b, fitX::color_b, 20, 1., fitX::color_b, 1, 1);
  xjjroot::sethempty(hmsig_b);

  TF1* fmsig_a = new TF1("fmsig_a", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fmsig_a, fitX::color_a, 7, 4);
  fmsig_a->SetParLimits(0, 0, 1.e+4);
  fmsig_a->SetParameter(0, nn);
  fmsig_a->SetParLimits(1, 3.686, 3.689);
  fmsig_a->SetParameter(1, 3.686097);
  fmsig_a->SetParameter(2, 0.001);
  fmsig_a->SetParLimits(2, 0, 0.01);
  TF1* fmsig_b = new TF1("fmsig_b", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fmsig_b, fitX::color_b, 7, 4);
  fmsig_b->SetParLimits(0, 0, 1.e+4);
  fmsig_b->SetParameter(0, nn);
  fmsig_b->SetParLimits(1, 3.865, 3.870);
  fmsig_b->SetParameter(1, 3.867);
  fmsig_b->SetParameter(2, 0.002);
  fmsig_b->SetParLimits(2, 0, 0.01);

  xjjroot::setgstyle();
  TCanvas* cmsig = new TCanvas("cmsig", "", 1200, 600);
  cmsig->Divide(2, 1);
  cmsig->cd(1);
  hmsig_a->Draw("pe");
  hmsig_a->Fit("fmsig_a", "Nq");
  hmsig_a->Fit("fmsig_a", "Nq");
  hmsig_a->Fit("fmsig_a", "LL");
  fitX::drawpull(hmsig_a, fmsig_a, fitX::color_a);
  labelsmc("#psi(2S)", fmsig_a->GetParameter(1), fmsig_a->GetParameter(2), "m", 4, " GeV");
  xjjroot::drawcomment(output.c_str());
  cmsig->cd(2);
  hmsig_b->Draw("pe");
  hmsig_b->Fit("fmsig_b", "Nq");
  hmsig_b->Fit("fmsig_b", "Nq");
  hmsig_b->Fit("fmsig_b", "LL");
  fitX::drawpull(hmsig_b, fmsig_b, fitX::color_b);
  labelsmc("X(3872)", fmsig_b->GetParameter(1), fmsig_b->GetParameter(2), "m", 4, " GeV");
  std::string outputname_msig = Form("plots/%s/ctoymc_msig.pdf", output.c_str());
  xjjroot::mkdir(outputname_msig);
  cmsig->SaveAs(outputname_msig.c_str());
  std::cout<<fmsig_a->GetParameter(2)<<" "<<fmsig_b->GetParameter(2)<<std::endl;;

  // --> minNLL

  TH1F* hminNll = (TH1F*)inf->Get("hminNll");
  hminNll->GetXaxis()->SetNdivisions(505);
  hminNll->GetXaxis()->SetTitle("-log(L)");
  hminNll->SetMaximum(hminNll->GetMaximum()*1.2);
  xjjroot::setthgrstyle(hminNll, kGray+3, 20, 1., kGray+3, 1, 1);
  xjjroot::sethempty(hminNll);

  TF1* fminNll = new TF1("fminNll", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  float xmin = hminNll->GetXaxis()->GetXmin(), xmax = hminNll->GetXaxis()->GetXmax();
  xjjroot::settfstyle(fminNll, kGray+3, 7, 5);
  fminNll->SetParLimits(0, 0, 1.e+4);
  fminNll->SetParameter(0, nn);
  fminNll->SetParLimits(1, xmin, xmax);
  fminNll->SetParameter(1, (xmin+xmax)/2.);
  fminNll->SetParameter(2, (xmax-xmin)/10.);
  fminNll->SetParLimits(2, 0, (xmax-xmin)/2.);

  xjjroot::setgstyle(2);
  TCanvas* cminNll = new TCanvas("cminNll", "", 600, 600);
  hminNll->Draw("pe");
  hminNll->Fit("fminNll", "Nq");
  hminNll->Fit("fminNll", "Nq");
  hminNll->Fit("fminNll", "LL");
  fitX::drawpull(hminNll, fminNll, kGray+2);
  labelsmc("Minimum likelihood", fminNll->GetParameter(1), fminNll->GetParameter(2), "l", 0, "");
  xjjroot::drawcomment(output.c_str());
  std::string outputname_minNll = Form("plots/%s/ctoymc_minNll.pdf", output.c_str());
  xjjroot::mkdir(outputname_minNll);
  cminNll->SaveAs(outputname_minNll.c_str());

}

void labelsmc(std::string label, double mean, double sigma, std::string var, int precision, std::string unit)
{
  xjjroot::drawtex(0.22, 0.85, label.c_str(), 0.04, 12, 62);
  xjjroot::drawtex(0.22, 0.85-0.05, fitX::pttag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.22, 0.85-0.05*2, fitX::ytag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.22, 0.85-0.05*3, fitX::centtag().c_str(), 0.04, 12, 42);
  std::string mtitle(Form("#bar{%s} = %.*f%s", var.c_str(), precision, mean, unit.c_str()));
  xjjroot::drawtex(0.65, 0.86, mtitle.c_str(), 0.038);
  std::string stitle(Form("#sigma_{%s} = %.*f%s", var.c_str(), precision, sigma, unit.c_str()));
  xjjroot::drawtex(0.65, 0.86-0.05, stitle.c_str(), 0.038);
  xjjroot::drawCMS();
}

int main(int argc, char* argv[])
{
  if(argc==3) { drawtoymc(argv[1], argv[2]); return 0; }
  return 1;
}
