#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

#include "xjjrootuti.h"
#include "systematics.h"
#include "fitX.h"

// std::vector<std::string> fitopt = {"poly4", "poly3", "cheb4", "cheb3", "3gaus", "floatwidth"};
std::vector<std::string> fitopt = {"poly4", "poly3", "range", "3gaus", "floatwidth"};
void drawpercent(TH1F* hh, float y=0.6);

void pdf_drawhist(std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  fitX::init(TFile::Open(std::string("rootfiles/"+output+"/fitX_fithist.root").c_str()));
  std::map<std::string, TH1F*> hyieldpromptCorr, hratio;
  std::vector<int> cc;
  TH1F* hhyieldpromptCorr_a = new TH1F("hhyieldpromptCorr_a", ";;N_{signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", fitopt.size(), 0, fitopt.size());
  xjjroot::sethempty(hhyieldpromptCorr_a, 0, 0.4);
  hhyieldpromptCorr_a->GetXaxis()->SetLabelSize(hhyieldpromptCorr_a->GetXaxis()->GetLabelSize()*1.4);
  xjjroot::setthgrstyle(hhyieldpromptCorr_a, fitX::color_a, 20, 1.2, fitX::color_a, 2, 2);
  TH1F* hhyieldpromptCorr_b = new TH1F("hhyieldpromptCorr_b", ";;N_{signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", fitopt.size(), 0, fitopt.size());
  xjjroot::sethempty(hhyieldpromptCorr_b, 0, 0.4);
  hhyieldpromptCorr_b->GetXaxis()->SetLabelSize(hhyieldpromptCorr_b->GetXaxis()->GetLabelSize()*1.4);
  xjjroot::setthgrstyle(hhyieldpromptCorr_b, fitX::color_b, 20, 1.2, fitX::color_b, 2, 2);
  TH1F* hhratio = new TH1F("hhratio", ";;R", fitopt.size(), 0, fitopt.size());
  xjjroot::sethempty(hhratio, 0, 0.2);
  hhratio->GetXaxis()->SetLabelSize(hhratio->GetXaxis()->GetLabelSize()*1.4);
  xjjroot::setthgrstyle(hhratio, kBlack, 21, 1.2, kBlack, 1, 2);
  for(int i=0; i<fitopt.size(); i++)
    {
      auto tt = fitopt[i];
      TFile* inf = TFile::Open(std::string("rootfiles/"+output+"/"+tt+"/fitX_fithist.root").c_str());
      hyieldpromptCorr[tt] = (TH1F*)inf->Get("hyieldpromptCorr");
      hyieldpromptCorr[tt]->SetName(Form("hyieldpromptCorr_%s", tt.c_str()));
      hratio[tt] = (TH1F*)inf->Get("hratio");
      hratio[tt]->SetName(Form("hratio_%s", tt.c_str()));
      xjjroot::sethempty(hyieldpromptCorr[tt]);
      xjjroot::setthgrstyle(hyieldpromptCorr[tt], xjjroot::mycolor_middle[xjjroot::cc[i]], 20, 1.2, xjjroot::mycolor_middle[xjjroot::cc[i]], 2, 2);
      xjjroot::sethempty(hratio[tt]);
      xjjroot::setthgrstyle(hratio[tt], xjjroot::mycolor_middle[xjjroot::cc[i]], 21, 1.2, xjjroot::mycolor_middle[xjjroot::cc[i]], 1, 2);
      cc.push_back(xjjroot::mycolor_middle[xjjroot::cc[i]]);
      hhyieldpromptCorr_a->GetXaxis()->SetBinLabel(i+1, tt.c_str());
      hhyieldpromptCorr_a->SetBinContent(i+1, hyieldpromptCorr[tt]->GetBinContent(fitX::ibin_a));
      hhyieldpromptCorr_a->SetBinError(i+1, hyieldpromptCorr[tt]->GetBinError(fitX::ibin_a));
      hhyieldpromptCorr_b->GetXaxis()->SetBinLabel(i+1, tt.c_str());
      hhyieldpromptCorr_b->SetBinContent(i+1, hyieldpromptCorr[tt]->GetBinContent(fitX::ibin_b));
      hhyieldpromptCorr_b->SetBinError(i+1, hyieldpromptCorr[tt]->GetBinError(fitX::ibin_b));
      hhratio->GetXaxis()->SetBinLabel(i+1, tt.c_str());
      hhratio->SetBinContent(i+1, hratio[tt]->GetBinContent(1));
      hhratio->SetBinError(i+1, hratio[tt]->GetBinError(1));
    }
  hhyieldpromptCorr_a->SetMinimum(0);
  hhyieldpromptCorr_a->SetMaximum(hhyieldpromptCorr_a->GetMaximum()*2.5);
  hhyieldpromptCorr_b->SetMinimum(0);
  hhyieldpromptCorr_b->SetMaximum(hhyieldpromptCorr_b->GetMaximum()*2.5);
  hhratio->SetMinimum(0);
  hhratio->SetMaximum(hhratio->GetMaximum()*3);

  gROOT->cd();

  xjjroot::setgstyle(1);

  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hhyieldpromptCorr_a->Draw("pe");
  drawpercent(hhyieldpromptCorr_a, 0.60);
  fitX::drawkinematics();
  xjjroot::drawtex(0.24, 0.81, fitX::title_a.c_str(), 0.042, 0.12, 62);
  xjjroot::drawCMS();
  c->cd(2);
  hhyieldpromptCorr_b->Draw("pe");
  drawpercent(hhyieldpromptCorr_b, 0.60);
  fitX::drawkinematics();
  xjjroot::drawtex(0.24, 0.81, fitX::title_b.c_str(), 0.042, 0.12, 62);
  xjjroot::drawCMS();
  // c->cd(3);
  // hhratio->Draw("pe");
  // drawpercent(hhratio, 0.60);
  // fitX::drawkinematics();
  // xjjroot::drawCMS();
  std::string outputname("plots/"+output+"/pdfvar/cpdfvar.pdf");
  xjjroot::mkdir(outputname);
  c->SaveAs(outputname.c_str());

}

int main(int argc, char* argv[])
{
  if(argc==2) { pdf_drawhist(argv[1]); return 0; }
  return 1;
}

void drawpercent(TH1F* hh, float y)
{
  int nn = hh->GetXaxis()->GetNbins();
  float lm = gPad->GetLeftMargin(), rm = gPad->GetRightMargin(), dm = (1-rm-lm)/nn;
  for(int i=0; i<nn; i++) 
    {
      float per = 100*(hh->GetBinContent(i+1)/hh->GetBinContent(1)-1);
      std::string tper(i?Form("%s%.1f%s", (per>0?"+":"-"), fabs(per), "%"):"Nominal");
      xjjroot::drawtex(lm + dm/2. + dm*i, y, tper.c_str(), 0.038, 22, 62, kGray+1);
    }
}
