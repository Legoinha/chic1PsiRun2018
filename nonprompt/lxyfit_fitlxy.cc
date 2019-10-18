#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>

#include <string>
#include <vector>

#include "var.h"
#include "fitX.h"

#include "xjjcuti.h"
#include "xjjrootuti.h"

float init_fprompt_a = 0.55;
float init_fprompt_b = 0.90;

TH1F *hhmcpdis_a, *hhmcnpdis_a;
TH1F *hhmcpdis_b, *hhmcnpdis_b;
Double_t funMix_a(Double_t* x_, Double_t* para)
{
  float x = x_[0];
  // float A = para[0];
  float R = para[0];
  float A = 1.;
  float ymcp = hhmcpdis_a->GetBinContent(hhmcpdis_a->GetXaxis()->FindBin(x));
  float ymcnp = hhmcnpdis_a->GetBinContent(hhmcnpdis_a->GetXaxis()->FindBin(x));
  return A*(R*ymcp+(1-R)*ymcnp);
}
Double_t funMix_b(Double_t* x_, Double_t* para)
{
  float x = x_[0];
  // float A = para[0];
  float R = para[0];
  float A = 1.;
  float ymcp = hhmcpdis_b->GetBinContent(hhmcpdis_b->GetXaxis()->FindBin(x));
  float ymcnp = hhmcnpdis_b->GetBinContent(hhmcnpdis_b->GetXaxis()->FindBin(x));
  return A*(R*ymcp+(1-R)*ymcnp);
}

void lxyfit_fitlxy(std::string input, std::string output, std::string type)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  output += (fitX::tagname()+"/"+type);
  datamc::var* vv = new datamc::var(type.c_str());
  if(!vv->valid()) return;

  hhmcpdis_a = (TH1F*)inf->Get("hmcpdis_a");
  hhmcnpdis_a = (TH1F*)inf->Get("hmcnpdis_a");
  TH1F* hdis_a = (TH1F*)inf->Get("hdis_a");
  hdis_a->SetMinimum(1.e-4);
  hdis_a->SetMaximum(hdis_a->GetMaximum()*1.e+2);
  hhmcpdis_b = (TH1F*)inf->Get("hmcpdis_b");
  hhmcnpdis_b = (TH1F*)inf->Get("hmcnpdis_b");
  TH1F* hdis_b = (TH1F*)inf->Get("hdis_b");
  hdis_b->SetMinimum(1.e-4);
  hdis_b->SetMaximum(hdis_b->GetMaximum()*1.e+2);

  TF1* fMix_a = new TF1("fMix_a", &funMix_a, vv->vars().front(), vv->vars().back(), 1);
  fMix_a->SetParameter(0, init_fprompt_a);
  fMix_a->SetParLimits(0, 0, 1);
  xjjroot::settfstyle(fMix_a, xjjroot::mycolor_middle["green"], 2, 3);
  TF1* fMix_b = new TF1("fMix_b", &funMix_b, vv->vars().front(), vv->vars().back(), 1);
  fMix_b->SetParameter(0, init_fprompt_b);
  fMix_b->SetParLimits(0, 0, 1);
  xjjroot::settfstyle(fMix_b, xjjroot::mycolor_middle["green"], 2, 3);

  xjjroot::setgstyle(1);
  TCanvas* c_a = new TCanvas("c_a", "", 600, 600);
  gPad->SetLogy();
  hdis_a->Draw("pe");
  TFitResultPtr fitResult_a = hdis_a->Fit("fMix_a","E S0Q", "", vv->vars().front(), vv->vars().back());
  fMix_a->Draw("same");
  float fprompt_a = fMix_a->GetParameter(0), errfprompt_a = fMix_a->GetParError(0);
  TH1F* hmcpdis_a = (TH1F*)hhmcpdis_a->Clone("hmcpdis_a");
  TH1F* hmcnpdis_a = (TH1F*)hhmcnpdis_a->Clone("hmcnpdis_a");
  hmcpdis_a->Scale(fprompt_a);
  hmcnpdis_a->Scale(1-fprompt_a);
  TH1F* hmcpnpdis_a = (TH1F*)hmcpdis_a->Clone("hmcpnpdis_a");
  hmcpnpdis_a->Add(hmcnpdis_a);
  xjjroot::setthgrstyle(hmcpnpdis_a, xjjroot::mycolor_middle["red"], 21, 0.5, xjjroot::mycolor_middle["red"], 1, 3, xjjroot::mycolor_light["red"], 1, 1001);
  xjjroot::setthgrstyle(hmcnpdis_a, xjjroot::mycolor_middle["azure"], 21, 0.5, xjjroot::mycolor_middle["azure"], 1, 3, xjjroot::mycolor_light["azure"], 1, 1001);
  TLegend* leg_a = new TLegend(0.60, 0.73-0.042*3, 0.89, 0.73);
  xjjroot::setleg(leg_a, 0.038);
  leg_a->AddEntry(hdis_a, "Signal Data", "p");
  leg_a->AddEntry(hmcpnpdis_a, "Prompt MC", "f");
  leg_a->AddEntry(hmcnpdis_a, "Nonprompt MC", "f");
  hmcpnpdis_a->Draw("hist same");
  hmcnpdis_a->Draw("hist same");
  hdis_a->Draw("pe same");  
  leg_a->Draw();
  xjjroot::drawCMS();
  fitX::drawkinematics();
  xjjroot::drawtex(0.23, 0.84, fitX::title_a.c_str(), 0.038, 12, 62);
  xjjroot::drawtex(0.23, 0.84-0.045, Form("f_{prompt} = %.3f #pm %.3f", fprompt_a, errfprompt_a), 0.038, 12, 62);
  gPad->RedrawAxis();
  std::string outputname_a(Form("plots/%s/cfitdis_a.pdf", output.c_str()));
  xjjroot::mkdir(outputname_a);
  c_a->SaveAs(outputname_a.c_str());

  TCanvas* c_b = new TCanvas("c_b", "", 600, 600);
  gPad->SetLogy();
  hdis_b->Draw("pe");
  TFitResultPtr fitResult_b = hdis_b->Fit("fMix_b","E S0Q", "", vv->vars().front(), vv->vars().back());
  fMix_b->Draw("same");
  float fprompt_b = fMix_b->GetParameter(0), errfprompt_b = fMix_b->GetParError(0);
  TH1F* hmcpdis_b = (TH1F*)hhmcpdis_b->Clone("hmcpdis_b");
  TH1F* hmcnpdis_b = (TH1F*)hhmcnpdis_b->Clone("hmcnpdis_b");
  hmcpdis_b->Scale(fprompt_b);
  hmcnpdis_b->Scale(1-fprompt_b);
  TH1F* hmcpnpdis_b = (TH1F*)hmcpdis_b->Clone("hmcpnpdis_b");
  hmcpnpdis_b->Add(hmcnpdis_b);
  xjjroot::setthgrstyle(hmcpnpdis_b, xjjroot::mycolor_middle["red"], 21, 0.5, xjjroot::mycolor_middle["red"], 1, 3, xjjroot::mycolor_light["red"], 1, 1001);
  xjjroot::setthgrstyle(hmcnpdis_b, xjjroot::mycolor_middle["azure"], 21, 0.5, xjjroot::mycolor_middle["azure"], 1, 3, xjjroot::mycolor_light["azure"], 1, 1001);
  TLegend* leg_b = new TLegend(0.60, 0.73-0.042*3, 0.89, 0.73);
  xjjroot::setleg(leg_b, 0.038);
  leg_b->AddEntry(hdis_b, "Signal Data", "p");
  leg_b->AddEntry(hmcpnpdis_b, "Prompt MC", "f");
  leg_b->AddEntry(hmcnpdis_b, "Nonprompt MC", "f");
  hmcpnpdis_b->Draw("hist same");
  hmcnpdis_b->Draw("hist same");
  hdis_b->Draw("pe same");  
  leg_b->Draw();
  xjjroot::drawCMS();
  fitX::drawkinematics();
  xjjroot::drawtex(0.23, 0.84, fitX::title_b.c_str(), 0.038, 12, 62);
  xjjroot::drawtex(0.23, 0.84-0.045, Form("f_{prompt} = %.3f #pm %.3f", fprompt_b, errfprompt_b), 0.038, 12, 62);
  gPad->RedrawAxis();
  std::string outputname_b(Form("plots/%s/cfitdis_b.pdf", output.c_str()));
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());

  std::cout<<"float fprompt_a = "<<fprompt_a<<", errfprompt_a = "<<fMix_a->GetParError(0)<<";"<<std::endl;
  std::cout<<"float fprompt_b = "<<fprompt_b<<", errfprompt_b = "<<fMix_b->GetParError(0)<<";"<<std::endl;

}

int main(int argc, char* argv[])
{
  if(argc==4) { lxyfit_fitlxy(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}
