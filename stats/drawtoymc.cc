#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>

#include <string>

#include "xjjrootuti.h"
#include "fit.h"

void labelsmc(std::string label, double mean, double sigma);
void drawtoymc(std::string input, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());

  TH1F* hysig_a = (TH1F*)inf->Get("hysig_a");
  xjjroot::setthgrstyle(hysig_a, fitX::color_a, 20, 1., fitX::color_a, 1, 1);
  xjjroot::sethempty(hysig_a);
  TH1F* hysig_b = (TH1F*)inf->Get("hysig_b");
  xjjroot::setthgrstyle(hysig_b, fitX::color_b, 20, 1., fitX::color_b, 1, 1);
  xjjroot::sethempty(hysig_b);
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

  TF1* fmsig_a = new TF1("fmsig_a", "[0]*TMath::Gaus(x, [1], [2])/TMath::Sqrt(2*3.1415926*[2])");
  xjjroot::settfstyle(fmsig_a, fitX::color_a, 7, 2);
  fmsig_a->SetParLimits(0, 0, 1.e+4);
  fmsig_a->SetParameter(0, 400);
  fmsig_a->SetParLimits(1, 3.686, 3.689);
  fmsig_a->SetParameter(1, 3.686097);
  fmsig_a->SetParameter(2, 0.001);
  TF1* fmsig_b = new TF1("fmsig_b", "[0]*TMath::Gaus(x, [1], [2])/TMath::Sqrt(2*3.1415926*[2])");
  xjjroot::settfstyle(fmsig_b, fitX::color_b, 7, 2);
  fmsig_b->SetParLimits(0, 0, 1.e+4);
  fmsig_b->SetParameter(0, 400);
  fmsig_b->SetParLimits(1, 3.865, 3.870);
  fmsig_b->SetParameter(1, 3.867);
  fmsig_b->SetParameter(2, 0.002);
  fmsig_b->SetParLimits(2, 0, 0.01);

  xjjroot::setgstyle();
  TCanvas* cm = new TCanvas("cm", "", 1200, 600);
  cm->Divide(2, 1);
  cm->cd(1);
  hmsig_a->Draw("pe");
  hmsig_a->Fit("fmsig_a", "Nq");
  hmsig_a->Fit("fmsig_a", "Nq");
  hmsig_a->Fit("fmsig_a", "LL");
  fitX::drawpull(hmsig_a, fmsig_a, fitX::color_a);
  labelsmc("#psi(2S)", fmsig_a->GetParameter(1), fmsig_a->GetParameter(2));
  fitX::drawcomment(output.c_str());
  cm->cd(2);
  hmsig_b->Draw("pe");
  hmsig_b->Fit("fmsig_b", "Nq");
  hmsig_b->Fit("fmsig_b", "Nq");
  hmsig_b->Fit("fmsig_b", "LL");
  fitX::drawpull(hmsig_b, fmsig_b, fitX::color_b);
  labelsmc("X(3872)", fmsig_b->GetParameter(1), fmsig_b->GetParameter(2));
  std::string outputname = Form("plots/%s/ctoymc_msig.pdf", output.c_str());
  xjjroot::mkdir(outputname);
  cm->SaveAs(outputname.c_str());
  std::cout<<fmsig_a->GetParameter(2)<<" "<<fmsig_b->GetParameter(2)<<std::endl;;

}

void labelsmc(std::string label, double mean, double sigma)
{
  xjjroot::drawtex(0.22, 0.85, label.c_str(), 0.04, 12, 62);
  xjjroot::drawtex(0.22, 0.85-0.05, fitX::pttag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.22, 0.85-0.05*2, fitX::ytag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.22, 0.85-0.05*3, fitX::centtag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.65, 0.86, Form("#bar{m} = %.4f GeV", mean), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05, Form("#sigma = %.4f GeV", sigma), 0.038);
  xjjroot::drawCMS();
}

int main(int argc, char* argv[])
{
  if(argc==3) drawtoymc(argv[1], argv[2]);
  
}
