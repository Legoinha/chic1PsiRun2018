#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TEfficiency.h>

#include "xjjrootuti.h"

void calcfprompt(std::string input_data, std::string input_Benrich,
                 std::string input_temp)
{
  TFile* inf_data = TFile::Open(input_data.c_str());
  TH1F* hdata = (TH1F*)inf_data->Get("hyield"); hdata->SetName("hdata");
  TFile* inf_Benrich = TFile::Open(input_Benrich.c_str());
  TH1F* hBenrich = (TH1F*)inf_Benrich->Get("hyield"); hBenrich->SetName("hBenrich");
  TFile* inf_temp = TFile::Open(input_temp.c_str());
  TH1F* hfgt = (TH1F*)inf_temp->Get("hfgt");

  TH1F* hnonprompt = (TH1F*)hBenrich->Clone("hnonprompt");
  hnonprompt->Divide(hfgt);
  TH1F* hprompt = (TH1F*)xjjroot::histMinusCorr(hdata, hnonprompt, "hprompt");
  // TH1F* hfprompt = (TH1F*)hprompt->Clone("hfprompt");
  // hfprompt->Divide(hdata);
  TEfficiency* grfprompt = new TEfficiency(*hprompt, *hdata); grfprompt->SetName("grfprompt");
  TH1F* hdatalt = (TH1F*)xjjroot::histMinusCorr(hdata, hBenrich, "hdatalt");
  // TH1F* hfpromptlt = (TH1F*)hprompt->Clone("hfprompt");
  // hfpromptlt->Divide(hdatalt);
  TEfficiency* grfpromptlt = new TEfficiency(*hprompt, *hdatalt); grfpromptlt->SetName("grfpromptlt");

  xjjroot::setthgrstyle(grfprompt, kGray+3, 20, 1.2, kGray+3, 4, 3, kGray+3, 0.1, 1001);
  xjjroot::setthgrstyle(grfpromptlt, kGray, 21, 1.2, kGray, 2, 3, kGray, 0.1, 1001);

  TH2F* hempty = new TH2F("hempty", ";;f_{prompt} After Cuts", 5, 0, 5, 10, 0, 1);
  xjjroot::sethempty(hempty, 0, 0.3);
  hempty->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hempty->GetXaxis()->SetBinLabel(4, "X(3872)");
  hempty->GetXaxis()->SetLabelSize(hempty->GetXaxis()->GetLabelSize()*1.5);
  xjjroot::setgstyle();
  TCanvas* cfp = new TCanvas("cfp", "", 600, 600);
  hempty->Draw();
  grfprompt->Draw("pe same");
  grfpromptlt->Draw("pe same");
  // xjjroot::drawtex(0.42, ysig1/ymax + 0.2, Form("%.0f #pm %.0f", ysig1, ysig1err), 0.042, 22, 62, fitX::color_a);
  // xjjroot::drawtex(0.72, ysig2/ymax + 0.2, Form("%.0f #pm %.0f", ysig2, ysig2err), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawtex(0.92, 0.84-0.55, "|y| < 1.5", 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.77-0.55, "p_{T} > 15 GeV/c", 0.042, 32, 62);
  xjjroot::drawtex(0.22, 0.84, "No l_{xy} Cut", 0.035, 12, 62, kGray+3);
  xjjroot::drawtex(0.22, 0.78, "l_{xy} < 0.1 mm", 0.035, 12, 62, kGray);
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  cfp->SaveAs("plots/cfprompt.pdf");

}

int main(int argc, char* argv[])
{
  if(argc==4) { calcfprompt(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}
