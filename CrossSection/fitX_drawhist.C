#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include "ppref/HEPData-ins1219950-v1/CMS_2013_R.h"
#include "ppref/HEPData-ins1495026-v1/ppATLAS.h"

#include "xjjrootuti.h"
#include "systematics.h"
#include "fitX.h"

void fitX_drawhist(std::string inputname, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = TFile::Open(inputname.c_str());
  fitX::init(inf);
  TH1F* hratio = (TH1F*)inf->Get("hratio");
  xjjroot::setthgrstyle(hratio, kBlack, 21, 1.2, kBlack, 1, 1);
  ppref::ppATLAS pprefATLAS("ppref/HEPData-ins1495026-v1");

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
  TH2F* hempty = new TH2F("hempty", ";p_{T};R", 10, 10, 70, 10, 0.01, 100);
  xjjroot::sethempty(hempty, 0, 0);
  hempty->Draw();
  xjjroot::drawline(10, 1, 70, 1, kGray+2, 2, 2);
  ppref::CMS_2013_R();
  pprefATLAS.Draw();
  gsyst->Draw("same 5");
  hratio->Draw("same ple");
  TLegend* legpp = new TLegend(0.63, 0.87-5*0.045, 0.63+0.30, 0.87);
  xjjroot::setleg(legpp, 0.035);
  legpp->AddEntry((TObject*)0, "", NULL);
  legpp->AddEntry(ppref::grae_syst, "Inclusive", "pf");
  legpp->AddEntry((TObject*)0, "", NULL);
  legpp->AddEntry(pprefATLAS.gg["promptRatio"]["syst"], "Prompt", "pf");
  legpp->AddEntry(pprefATLAS.gg["nonpromptRatio"]["syst"], "Nonprompt", "pf");
  legpp->Draw();
  xjjroot::drawtex(0.63+0.01, 0.87-0.01, "pp (7 TeV) |y| < 1.2", 0.035, 13);
  xjjroot::drawtex(0.63+0.01, 0.87-0.01-0.045*2, "pp (8 TeV) |y| < 0.75", 0.035, 13);
  TLegend* legpbpb = new TLegend(0.22, 0.87-5*0.045, 0.22+0.30, 0.87-2*0.045);
  xjjroot::setleg(legpbpb, 0.035);
  legpbpb->AddEntry((TObject*)0, "", NULL);
  legpbpb->AddEntry((TObject*)0, "", NULL);
  legpbpb->AddEntry(gsyst, "Prompt", "pf");
  legpbpb->Draw();
  xjjroot::drawtex(0.22+0.01, 0.87-0.01-0.045*2, Form("PbPb (5.02 TeV) %s", fitX::ytag().c_str()), 0.035, 13);
  xjjroot::drawtex(0.22+0.01, 0.87-0.01-0.045*3, Form("Centrality %.0f-%.0f%s", fitX::centmincut, fitX::centmaxcut, "%"), 0.035, 13);
  xjjroot::drawCMSleft("", 0.05, -0.08);
  xjjroot::drawCMSright("1.5 nb^{-1} (2018 PbPb 5.02 TeV)");
  xjjroot::drawcomment(output.c_str(), "r");
  xjjroot::mkdir(Form("plots/%s/cratio.pdf", output.c_str()));
  cratio->SaveAs(Form("plots/%s/cratio.pdf", output.c_str()));
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  // fitX::init(TFile::Open(argv[1]));
  // std::string dirname = std::string(argv[2])+fitX::tagname();
  // if(argc==3) { fitX_drawhist(argv[1], dirname); return 0; }
  if(argc==3) { fitX_drawhist(argv[1], argv[2]); return 0; }
  return 1;
}
