#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include "ppref/CMS_2013_R.h"

#include "xjjrootuti.h"
#include "systematics.h"

void fitX_drawhist(std::string inputname, std::string output)
{
  TFile* inf = TFile::Open(Form("%s.root", inputname.c_str()));
  TH1F* hratio = (TH1F*)inf->Get("hratio");
  xjjroot::setthgrstyle(hratio, kBlack, 21, 1.2, kBlack, 1, 1);

  std::vector<float> xx, yy, xel, xeh, yel, yeh;
  for(int i=0; i<hratio->GetNbinsX(); i++)
    {
      std::cout<<hratio->GetBinCenter(i+1)<<" "<<hratio->GetBinContent(i+1)<<std::endl;
      xx.push_back(hratio->GetBinCenter(i+1));
      xel.push_back(1.);
      xeh.push_back(1.);
      yy.push_back(hratio->GetBinContent(i+1));
      yel.push_back(syst::getsyst(2, "d")*hratio->GetBinContent(i+1));
      yeh.push_back(syst::getsyst(2, "u")*hratio->GetBinContent(i+1));
    }
  TGraphAsymmErrors* gsyst = new TGraphAsymmErrors(hratio->GetNbinsX(), xx.data(), yy.data(), xel.data(), xeh.data(), yel.data(), yeh.data());
  gsyst->SetName("gr_ratio_syst");
  xjjroot::setthgrstyle(gsyst, kBlack, 21, 1.2, 0, 0, 0, kGray+1, 0.5, 1001);

  xjjroot::setgstyle(2);
  TCanvas* cratio = new TCanvas("cratio", "", 600, 600);
  cratio->SetLogy();
  TH2F* hempty = new TH2F("hempty", ";p_{T};R", 10, 10, 50, 10, 0.01, 100);
  xjjroot::sethempty(hempty, 0, 0);
  hempty->Draw();
  xjjroot::drawline(10, 1, 50, 1, kGray+2, 2, 2);
  ppref::CMS_2013_R();
  gsyst->Draw("same 5");
  hratio->Draw("same ple");
  xjjroot::drawCMSleft("", 0.05, -0.08);
  TLegend* leg = new TLegend(0.60, 0.87-4*0.05, 0.95, 0.87);
  xjjroot::setleg(leg, 0.042);
  leg->AddEntry((TObject*)0, "", NULL);
  leg->AddEntry(ppref::grae_syst, "|y| < 1.2", "pf");
  leg->AddEntry((TObject*)0, "", NULL);
  leg->AddEntry(gsyst, "|y| < 1.6", "pf");
  leg->Draw();
  xjjroot::drawtex(0.60+0.01, 0.87-0.01, "pp (7 TeV)", 0.042, 13);
  xjjroot::drawtex(0.60+0.01, 0.87-0.01-0.05*2, "PbPb (5.02 TeV)", 0.042, 13);
  cratio->SaveAs(Form("plots/%s/cratio.pdf", output.c_str()));

}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX_drawhist(argv[1], argv[2]); return 0; }
  return 1;
}
