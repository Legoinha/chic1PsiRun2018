#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>

#include <string>

#include "xjjrootuti.h"
#include "fit.h"

void labelsmc(std::string label, double mean, double sigma, std::string var, int precision, std::string unit);
void drawtoymc(std::string input, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = new TFile(input.c_str());
  TTree* nt = (TTree*)inf->Get("nt");
  float ysig_a; nt->SetBranchAddress("ysig_a", &ysig_a);
  float ysigerr_a; nt->SetBranchAddress("ysigerr_a", &ysigerr_a);
  float ysig_b; nt->SetBranchAddress("ysig_b", &ysig_b);
  float ysigerr_b; nt->SetBranchAddress("ysigerr_b", &ysigerr_b);
  float msig_a; nt->SetBranchAddress("msig_a", &msig_a);
  float msigerr_a; nt->SetBranchAddress("msigerr_a", &msigerr_a);
  float msig_b; nt->SetBranchAddress("msig_b", &msig_b);
  float msigerr_b; nt->SetBranchAddress("msigerr_b", &msigerr_b);
  float minNll; nt->SetBranchAddress("minNll", &minNll);
  
  // nominal result
  nt->GetEntry(0);
  float no_ysig_a = ysig_a;
  float no_ysigerr_a = ysigerr_a;
  float no_ysig_b = ysig_b;
  float no_ysigerr_b = ysigerr_b;
  float no_msig_a = msig_a;
  float no_msigerr_a = msigerr_a;
  float no_msig_b = msig_b;
  float no_msigerr_b = msigerr_b;
  float no_minNll = minNll;

  // toy MC
  std::map<std::string, TH1F*> hh;
  hh["ysig_a"] = new TH1F("hysig_a", ";(N - #bar{N})/#sigma_{N};", 40, -5, 5);
  hh["ysig_b"] = new TH1F("hysig_b", ";(N - #bar{N})/#sigma_{N};", 40, -5, 5);
  hh["msig_a"] = new TH1F("hmsig_a", ";(m - #bar{m})/#sigma_{m};", 40, -5, 5);
  hh["msig_b"] = new TH1F("hmsig_b", ";(m - #bar{m})/#sigma_{m};", 40, -5, 5);
  hh["minNll"] = new TH1F("hminNll", ";-log(L);", 40, no_minNll - 100 , no_minNll + 100); //

  int nn = nt->GetEntries()-1;
  for(int i=1; i<=nn; i++)
    {
      nt->GetEntry(i);
      hh["ysig_a"]->Fill((ysig_a-no_ysig_a)/no_ysigerr_a);
      hh["ysig_b"]->Fill((ysig_b-no_ysig_b)/no_ysigerr_b);
      hh["msig_a"]->Fill((msig_a-no_msig_a)/no_msigerr_a);
      hh["msig_b"]->Fill((msig_b-no_msig_b)/no_msigerr_b);
      hh["minNll"]->Fill(minNll);
    }

  for(auto& h : hh)
    {
      h.second->GetXaxis()->SetNdivisions(505);
      h.second->SetMaximum(h.second->GetMaximum()*1.2);
      xjjroot::sethempty(h.second);
    }

  // --> ysig

  xjjroot::setthgrstyle(hh["ysig_a"], fitX::color_a, 20, 1., fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(hh["ysig_b"], fitX::color_b, 20, 1., fitX::color_b, 1, 1);

  TF1* fysig_a = new TF1("fysig_a", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fysig_a, fitX::color_a, 7, 4);
  fysig_a->SetParLimits(0, 0, 1.e+4);
  fysig_a->SetParameter(0, nn);
  fysig_a->SetParLimits(1, -5, 5); // 0-200
  fysig_a->SetParameter(1, 0); // 100
  fysig_a->SetParameter(2, 1); // 30
  fysig_a->SetParLimits(2, 0, 3); // 0-60
  TF1* fysig_b = new TF1("fysig_b", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fysig_b, fitX::color_b, 7, 4);
  fysig_b->SetParLimits(0, 0, 1.e+4);
  fysig_b->SetParameter(0, nn);
  fysig_b->SetParLimits(1, -5, 5); // 0-200
  fysig_b->SetParameter(1, 0); // 1
  fysig_b->SetParameter(2, 1); // 30
  fysig_b->SetParLimits(2, 0, 3); //0-60

  xjjroot::setgstyle();
  TCanvas* cysig = new TCanvas("cysig", "", 1200, 600);
  cysig->Divide(2, 1);
  cysig->cd(1);
  hh["ysig_a"]->Draw("pe");
  hh["ysig_a"]->Fit("fysig_a", "Nq");
  hh["ysig_a"]->Fit("fysig_a", "Nq");
  hh["ysig_a"]->Fit("fysig_a", "LL");
  fitX::drawpull(hh["ysig_a"], fysig_a, fitX::color_a);
  labelsmc("#psi(2S)", fysig_a->GetParameter(1), fysig_a->GetParameter(2), "N", 3, "");
  xjjroot::drawcomment(output.c_str());
  cysig->cd(2);
  hh["ysig_b"]->Draw("pe");
  hh["ysig_b"]->Fit("fysig_b", "Nq");
  hh["ysig_b"]->Fit("fysig_b", "Nq");
  hh["ysig_b"]->Fit("fysig_b", "LL");
  fitX::drawpull(hh["ysig_b"], fysig_b, fitX::color_b);
  labelsmc("X(3872)", fysig_b->GetParameter(1), fysig_b->GetParameter(2), "N", 3, "");
  std::string outputname_ysig = Form("plots/%s/ctoymc_ysig.pdf", output.c_str());
  xjjroot::mkdir(outputname_ysig);
  cysig->SaveAs(outputname_ysig.c_str());
  std::cout<<fysig_a->GetParameter(2)<<" "<<fysig_b->GetParameter(2)<<std::endl;;

  // --> msig

  xjjroot::setthgrstyle(hh["msig_a"], fitX::color_a, 20, 1., fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(hh["msig_b"], fitX::color_b, 20, 1., fitX::color_b, 1, 1);

  TF1* fmsig_a = new TF1("fmsig_a", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fmsig_a, fitX::color_a, 7, 4);
  fmsig_a->SetParLimits(0, 0, 1.e+4);
  fmsig_a->SetParameter(0, nn);
  fmsig_a->SetParLimits(1, -5, 5);
  fmsig_a->SetParameter(1, 0);
  fmsig_a->SetParameter(2, 1);
  fmsig_a->SetParLimits(2, 0, 3);
  TF1* fmsig_b = new TF1("fmsig_b", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  xjjroot::settfstyle(fmsig_b, fitX::color_b, 7, 4);
  fmsig_b->SetParLimits(0, 0, 1.e+4);
  fmsig_b->SetParameter(0, nn);
  fmsig_b->SetParLimits(1, -5, 5);
  fmsig_b->SetParameter(1, 0);
  fmsig_b->SetParameter(2, 1);
  fmsig_b->SetParLimits(2, 0, 3);

  xjjroot::setgstyle();
  TCanvas* cmsig = new TCanvas("cmsig", "", 1200, 600);
  cmsig->Divide(2, 1);
  cmsig->cd(1);
  hh["msig_a"]->Draw("pe");
  hh["msig_a"]->Fit("fmsig_a", "Nq");
  hh["msig_a"]->Fit("fmsig_a", "Nq");
  hh["msig_a"]->Fit("fmsig_a", "LL");
  fitX::drawpull(hh["msig_a"], fmsig_a, fitX::color_a);
  labelsmc("#psi(2S)", fmsig_a->GetParameter(1), fmsig_a->GetParameter(2), "m", 3, "");
  xjjroot::drawcomment(output.c_str());
  cmsig->cd(2);
  hh["msig_b"]->Draw("pe");
  hh["msig_b"]->Fit("fmsig_b", "Nq");
  hh["msig_b"]->Fit("fmsig_b", "Nq");
  hh["msig_b"]->Fit("fmsig_b", "LL");
  fitX::drawpull(hh["msig_b"], fmsig_b, fitX::color_b);
  labelsmc("X(3872)", fmsig_b->GetParameter(1), fmsig_b->GetParameter(2), "m", 3, "");
  std::string outputname_msig = Form("plots/%s/ctoymc_msig.pdf", output.c_str());
  xjjroot::mkdir(outputname_msig);
  cmsig->SaveAs(outputname_msig.c_str());
  std::cout<<fmsig_a->GetParameter(2)<<" "<<fmsig_b->GetParameter(2)<<std::endl;;

  // --> minNLL
  xjjroot::setthgrstyle(hh["minNll"], kGray+3, 20, 1., kGray+3, 1, 1);

  TF1* fminNll = new TF1("fminNll", "[0]*TMath::Gaus(x, [1], [2])/(TMath::Sqrt(2*3.1415926)*[2])");
  float xmin = hh["minNll"]->GetXaxis()->GetXmin(), xmax = hh["minNll"]->GetXaxis()->GetXmax();
  xjjroot::settfstyle(fminNll, kGray+3, 7, 5);
  fminNll->SetParLimits(0, 0, 1.e+4);
  fminNll->SetParameter(0, nn);
  fminNll->SetParLimits(1, xmin, xmax);
  fminNll->SetParameter(1, (xmin+xmax)/2.);
  fminNll->SetParameter(2, (xmax-xmin)/10.);
  fminNll->SetParLimits(2, 0, (xmax-xmin)/2.);

  xjjroot::setgstyle(2);
  TCanvas* cminNll = new TCanvas("cminNll", "", 600, 600);
  hh["minNll"]->Draw("pe");
  hh["minNll"]->Fit("fminNll", "Nq");
  hh["minNll"]->Fit("fminNll", "Nq");
  hh["minNll"]->Fit("fminNll", "LL");
  fitX::drawpull(hh["minNll"], fminNll, kGray+2);
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
  // std::string mtitle(Form("#bar{%s} = %.*f%s", var.c_str(), precision, mean, unit.c_str()));
  std::string mtitle(Form("mean = %.*f%s", precision, mean, unit.c_str()));
  xjjroot::drawtex(0.65, 0.86, mtitle.c_str(), 0.038);
  // std::string stitle(Form("#sigma_{%s} = %.*f%s", var.c_str(), precision, sigma, unit.c_str()));
  std::string stitle(Form("width = %.*f%s", precision, sigma, unit.c_str()));
  xjjroot::drawtex(0.65, 0.86-0.05, stitle.c_str(), 0.038);
  xjjroot::drawCMS();
}

int main(int argc, char* argv[])
{
  if(argc==3) { drawtoymc(argv[1], argv[2]); return 0; }
  return 1;
}
