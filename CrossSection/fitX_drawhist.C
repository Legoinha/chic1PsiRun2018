#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include "ppref/HEPData-ins1219950-v1/CMS_2013_Rprompt.h"
#include "ppref/HEPData-ins1219950-v1/CMS_2013_R.h"
#include "ppref/HEPData-ins1495026-v1/ppATLAS.h"

#include "xjjrootuti.h"
#include "systematics.h"
#include "fitX.h"
#include "results.h"

void print(TFile* inf);
void fitX_drawhist(std::string inputname, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  TFile* inf = TFile::Open(inputname.c_str());
  fitX::init(inf);
  TH1F* hratio = (TH1F*)inf->Get("hratio");
  xjjroot::setthgrstyle(hratio, xjjroot::mycolor_satmiddle["red"], 21, 1.2, xjjroot::mycolor_satmiddle["red"], 1, 1);
  ppref::ppATLAS pprefATLAS("ppref/HEPData-ins1495026-v1");
  ppref::CMS_2013_R pprefCMS;
  ppref::CMS_2013_Rprompt pprefCMSprompt;

  std::vector<float> xx, yy, xel, xeh, yel, yeh;
  for(int i=0; i<hratio->GetNbinsX(); i++)
    {
      // std::cout<<hratio->GetBinCenter(i+1)<<" "<<hratio->GetBinContent(i+1)<<std::endl;
      xx.push_back(hratio->GetBinCenter(i+1));
      xel.push_back(1.);
      xeh.push_back(1.);
      yy.push_back(hratio->GetBinContent(i+1));
      yel.push_back(syst::getsyst(2, "d")*hratio->GetBinContent(i+1));
      yeh.push_back(syst::getsyst(2, "u")*hratio->GetBinContent(i+1));
      // std::cout<<i<<" "<<xx[i]<<" "<<hratio->GetBinWidth(i+1)/2.<<" "<<yy[i]<<" "<<hratio->GetBinError(i+1)<<" "<<yel[i]<<" "<<yeh[i]<<std::endl;
    }
  TGraphAsymmErrors* gsyst = new TGraphAsymmErrors(hratio->GetNbinsX(), xx.data(), yy.data(), xel.data(), xeh.data(), yel.data(), yeh.data());
  gsyst->SetName("gr_ratio_syst");
  xjjroot::setthgrstyle(gsyst, xjjroot::mycolor_satmiddle["red"], 21, 1.2, 0, 0, 0, xjjroot::mycolor_middle["red"], 0.5, 1001);

  TH2F* hemptypaperlog = new TH2F("hemptypaperlog", ";p_{T};R = #frac{N^{X(3872)#rightarrowJ/#psi#pi#pi}}{N^{#psi(2S)#rightarrowJ/#psi#pi#pi}}", 10, 10, 70, 10, 0.04, 20);
  xjjroot::sethempty(hemptypaperlog, 0, 0);
  TH2F* hemptylog = new TH2F("hemptylog", ";p_{T};R = #frac{N^{X(3872)#rightarrowJ/#psi#pi#pi}}{N^{#psi(2S)#rightarrowJ/#psi#pi#pi}}", 10, 10, 70, 10, 0.02, 50);
  xjjroot::sethempty(hemptylog, 0, 0);
  TH2F* hemptylinear = new TH2F("hemptylinear", ";p_{T};R = #frac{N^{X(3872)#rightarrowJ/#psi#pi#pi}}{N^{#psi(2S)#rightarrowJ/#psi#pi#pi}}", 10, 10, 70, 10, 0, 2.5);
  xjjroot::sethempty(hemptylinear, 0, 0);

  // --> paper leg <--
  float linesp = 0.045, textsp = 0.035;
  TLegend* legpaper_pp = new TLegend(0.57, 0.82-4*linesp, 0.57+0.30, 0.82);
  xjjroot::setleg(legpaper_pp, textsp);
  legpaper_pp->AddEntry(pprefCMSprompt.grae_syst(), "#bf{pp} (7 TeV)", "pf");
  legpaper_pp->AddEntry((TObject*)0, "|y| < 1.2 (CMS)", NULL);
  legpaper_pp->AddEntry(pprefATLAS.gg["promptRatio"]["syst"], "#bf{pp} (8 TeV)", "pf");
  legpaper_pp->AddEntry((TObject*)0, "|y| < 0.75 (ATLAS)", NULL);
  TLegend* legpaper_pbpb = new TLegend(0.22, 0.82-4*linesp, 0.22+0.30, 0.82);
  xjjroot::setleg(legpaper_pbpb, textsp);
  legpaper_pbpb->AddEntry((TObject*)0, "", NULL);
  legpaper_pbpb->AddEntry((TObject*)0, "", NULL);
  legpaper_pbpb->AddEntry(gsyst, "#bf{PbPb} (5.02 TeV)", "pf");
  legpaper_pbpb->AddEntry((TObject*)0, "|y| < 1.6, 0-90\%", NULL);

  // --> paper canvas <--
  xjjroot::setgstyle(2);
  TCanvas* cratio = new TCanvas("cratiolog", "", 600, 600);
  cratio->SetLogy();
  hemptypaperlog->Draw();
  xjjroot::drawline(10, 1, 70, 1, kGray+2, 2, 2);
  pprefCMSprompt.Draw();
  pprefATLAS.Draw("prompt");
  gsyst->Draw("same 5");
  hratio->Draw("same ple");
  legpaper_pp->Draw();
  legpaper_pbpb->Draw();
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}}", 0.05, -0.08);
  xjjroot::drawCMSleft("#it{Prompt}", 0.05, -0.13);
  xjjroot::drawCMSright("1.7 nb^{-1} (PbPb 5.02 TeV)");
  xjjroot::mkdir(Form("plots/%s/cratiolog.pdf", output.c_str()));
  cratio->SaveAs(Form("plots/%s/cratiolog.pdf", output.c_str()));

  cratio = new TCanvas("cratiolinear", "", 600, 600);
  hemptylinear->Draw();
  xjjroot::drawline(10, 1, 70, 1, kGray+2, 2, 2);
  pprefCMSprompt.Draw();
  pprefATLAS.Draw("prompt");
  gsyst->Draw("same 5");
  hratio->Draw("same ple");
  legpaper_pp->Draw();
  legpaper_pbpb->Draw();
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}}", 0.05, -0.08);
  xjjroot::drawCMSleft("#it{Prompt}", 0.05, -0.13);
  xjjroot::drawCMSright("1.7 nb^{-1} (PbPb 5.02 TeV)");
  xjjroot::mkdir(Form("plots/%s/cratiolinear.pdf", output.c_str()));
  cratio->SaveAs(Form("plots/%s/cratiolinear.pdf", output.c_str()));

  // --> supp leg <--
  float linesuppsp = 0.042, textsuppsp = 0.033;
  TLegend* legsupp_pp = new TLegend(0.62, 0.88-8*linesuppsp, 0.62+0.30, 0.88);
  xjjroot::setleg(legsupp_pp, textsuppsp);
  legsupp_pp->AddEntry((TObject*)0, "", NULL);
  legsupp_pp->AddEntry((TObject*)0, "", NULL);
  // legsupp_pp->AddEntry(pprefCMSprompt.grae_syst(), "Prompt", "pf");
  legsupp_pp->AddEntry(pprefCMS.grae_syst(), "Inclusive", "pf");
  legsupp_pp->AddEntry((TObject*)0, "", NULL);
  legsupp_pp->AddEntry((TObject*)0, "", NULL);
  legsupp_pp->AddEntry((TObject*)0, "", NULL);
  legsupp_pp->AddEntry(pprefATLAS.gg["promptRatio"]["syst"], "Prompt", "pf");
  legsupp_pp->AddEntry(pprefATLAS.gg["nonpromptRatio"]["syst"], "Nonprompt", "pf");
  TLegend* legsupp_pbpb = new TLegend(0.22, 0.88-5*linesuppsp, 0.22+0.30, 0.88-2*linesuppsp);
  xjjroot::setleg(legsupp_pbpb, textsuppsp);
  legsupp_pbpb->AddEntry((TObject*)0, "", NULL);
  legsupp_pbpb->AddEntry((TObject*)0, "", NULL);
  legsupp_pbpb->AddEntry(gsyst, "Prompt", "pf");

  auto legsupp = [](float linesuppsp, float textsuppsp) // define lambda: like function inside function
  {
    xjjroot::drawtex(0.63+0.01, 0.88-0.01, "#bf{pp} (7 TeV, CMS)", textsuppsp, 13);
    xjjroot::drawtex(0.63+0.01, 0.88-0.01-linesuppsp, "|y| < 1.2", textsuppsp, 13);
    xjjroot::drawtex(0.63+0.01, 0.88-0.01-linesuppsp*4, "#bf{pp} (8 TeV, ATLAS)", textsuppsp, 13);
    xjjroot::drawtex(0.63+0.01, 0.88-0.01-linesuppsp*5, "|y| < 0.75", textsuppsp, 13);
    xjjroot::drawtex(0.22+0.01, 0.88-0.01-linesuppsp*2, Form("#bf{PbPb} (5.02 TeV, CMS)"), textsuppsp, 13);
    xjjroot::drawtex(0.22+0.01, 0.88-0.01-linesuppsp*3, Form("%s, %.0f-%.0f%s", fitX::ytag().c_str(), fitX::centmincut, fitX::centmaxcut, "%"), textsuppsp, 13);
  };

  // --> supp canvas <--
  cratio = new TCanvas("cratiosupplog", "", 600, 600);
  cratio->SetLogy();
  hemptylog->Draw();
  xjjroot::drawline(10, 1, 70, 1, kGray+2, 2, 2);
  // pprefCMSprompt.Draw();
  pprefCMS.Draw();
  pprefATLAS.Draw();
  gsyst->Draw("same 5");
  hratio->Draw("same ple");
  legsupp_pp->Draw();
  legsupp_pbpb->Draw();
  legsupp(linesuppsp, textsuppsp);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}}", 0.05, -0.08);
  xjjroot::drawCMSleft("#it{Supplementary}", 0.16, -0.08); // preliminary
  xjjroot::drawCMSright("1.7 nb^{-1} (PbPb 5.02 TeV)");
  xjjroot::mkdir(Form("plots/%s/cratiosupplog.pdf", output.c_str()));
  cratio->SaveAs(Form("plots/%s/cratiosupplog.pdf", output.c_str()));

  cratio = new TCanvas("cratiosupplinear", "", 600, 600);
  hemptylinear->Draw();
  xjjroot::drawline(10, 1, 70, 1, kGray+2, 2, 2);
  // pprefCMSprompt.Draw();
  pprefCMS.Draw();
  pprefATLAS.Draw();
  gsyst->Draw("same 5");
  hratio->Draw("same ple");
  legsupp_pp->Draw();
  legsupp_pbpb->Draw();
  legsupp(linesuppsp, textsuppsp);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}}", 0.05, -0.08);
  xjjroot::drawCMSleft("#it{Supplementary}", 0.16, -0.08); // preliminary
  xjjroot::drawCMSright("1.7 nb^{-1} (PbPb 5.02 TeV)");
  xjjroot::mkdir(Form("plots/%s/cratiosupplinear.pdf", output.c_str()));
  cratio->SaveAs(Form("plots/%s/cratiosupplinear.pdf", output.c_str()));

  std::cout<<std::endl;
  fitX::results rs(inf, output.c_str());
  rs.print(syst::getsyst(0, "u", "q"), syst::getsyst(1, "u", "q"), syst::getsyst(2, "u", "q"));
}

int main(int argc, char* argv[])
{
  // fitX::init(TFile::Open(argv[1]));
  // std::string dirname = std::string(argv[2])+fitX::tagname();
  // if(argc==3) { fitX_drawhist(argv[1], dirname); return 0; }
  if(argc==3) { fitX_drawhist(argv[1], argv[2]); return 0; }
  return 1;
}
