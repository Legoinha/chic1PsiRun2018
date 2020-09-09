#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>

#include <string>
#include <vector>
#include <map>

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

std::map<std::string, float> fitlxy(std::string input, std::string output, std::string type, int ii, float dsmear)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  std::map<std::string, float> rs;

  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  output += (fitX::tagname()+"/"+type);
  datamc::var* vv = new datamc::var(type.c_str());
  if(!vv->valid()) return rs;

  hhmcpdis_a = (TH1F*)inf->Get(Form("hmcpdis_a_%d", ii));
  hhmcnpdis_a = (TH1F*)inf->Get(Form("hmcnpdis_a_%d", ii));
  TH1F* hdis_a = (TH1F*)inf->Get("hdis_a");
  hdis_a->SetMinimum(1.e-4);
  hdis_a->SetMaximum(hdis_a->GetMaximum()*1.e+2);
  hhmcpdis_b = (TH1F*)inf->Get(Form("hmcpdis_b_%d", ii));
  hhmcnpdis_b = (TH1F*)inf->Get(Form("hmcnpdis_b_%d", ii));
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
  xjjroot::drawtex(0.23, 0.84-0.045*2, Form("#chi^{2} Prob = %.1f%s", TMath::Prob(fMix_a->GetChisquare(), fMix_a->GetNDF())*100., "%"), 0.038, 12, 62);
  xjjroot::drawtex(0.23, 0.84-0.045*3, Form("#alpha = %.2f", 1.0+dsmear*ii), 0.038, 12, 62);
  gPad->RedrawAxis();
  std::string outputname_a(Form("plots/%s/smear/cfitdis_a_%d.pdf", output.c_str(), ii));
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
  xjjroot::drawtex(0.23, 0.84-0.045*2, Form("#chi^{2} Prob = %.1f%s", TMath::Prob(fMix_b->GetChisquare(), fMix_b->GetNDF())*100., "%"), 0.038, 12, 62);
  xjjroot::drawtex(0.23, 0.84-0.045*3, Form("#alpha = %.2f", 1.0+dsmear*ii), 0.038, 12, 62);
  gPad->RedrawAxis();
  std::string outputname_b(Form("plots/%s/smear/cfitdis_b_%d.pdf", output.c_str(), ii));
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());

  std::cout<<Form("#alpha = %.2f", 1.0+dsmear*ii)<<std::endl;
  std::cout<<Form("_a: #chi^{2} Prob = %.3f%s", TMath::Prob(fMix_a->GetChisquare(), fMix_a->GetNDF())*100., "%")<<std::endl;
  std::cout<<"float fprompt_a = "<<fprompt_a<<", errfprompt_a = "<<fMix_a->GetParError(0)<<";"<<std::endl;
  std::cout<<Form("_b: #chi^{2} Prob = %.3f%s", TMath::Prob(fMix_b->GetChisquare(), fMix_b->GetNDF())*100., "%")<<std::endl;
  std::cout<<"float fprompt_b = "<<fprompt_b<<", errfprompt_b = "<<fMix_b->GetParError(0)<<";"<<std::endl;

  rs["alpha"] = 1.0+dsmear*ii;
  rs["prob_a"] = TMath::Prob(fMix_a->GetChisquare(), fMix_a->GetNDF());
  rs["prob_b"] = TMath::Prob(fMix_b->GetChisquare(), fMix_b->GetNDF());
  rs["fprompt_a"] = fprompt_a;
  rs["fprompt_b"] = fprompt_b;

  return rs;
}

void lxyfit_fitlxy(std::string input, std::string output, std::string type, int nsmear, float dsmear)
{
  std::vector<float> valpha, vprob_a, vprob_b, vfprompt_a, vfprompt_b;
  for(int i=0; i<nsmear; i++)
    {
      std::map<std::string, float> rs = fitlxy(input, output, type, i, dsmear);
      valpha.push_back(rs["alpha"]);
      vprob_a.push_back(rs["prob_a"]);
      vprob_b.push_back(rs["prob_b"]);
      vfprompt_a.push_back(rs["fprompt_a"]);
      vfprompt_b.push_back(rs["fprompt_b"]);
    }
  TGraph* gr_prob_alpha_a = new TGraph(nsmear, valpha.data(), vprob_a.data());
  gr_prob_alpha_a->SetName("gr_prob_alpha_a");
  xjjroot::setthgrstyle(gr_prob_alpha_a, fitX::color_a, 21, 0.9, fitX::color_a, 1, 2);
  TGraph* gr_prob_alpha_b = new TGraph(nsmear, valpha.data(), vprob_b.data());
  gr_prob_alpha_b->SetName("gr_prob_alpha_b");
  xjjroot::setthgrstyle(gr_prob_alpha_b, fitX::color_b, 21, 0.9, fitX::color_b, 1, 2);

  TGraph* gr_prob_fprompt_a = new TGraph(nsmear, vfprompt_a.data(), vprob_a.data());
  gr_prob_fprompt_a->SetName("gr_prob_fprompt_a");
  xjjroot::setthgrstyle(gr_prob_fprompt_a, fitX::color_a, 21, 0.9, fitX::color_a, 1, 2);
  TGraph* gr_prob_fprompt_b = new TGraph(nsmear, vfprompt_b.data(), vprob_b.data());
  gr_prob_fprompt_b->SetName("gr_prob_fprompt_b");
  xjjroot::setthgrstyle(gr_prob_fprompt_b, fitX::color_b, 21, 0.9, fitX::color_b, 1, 2);

  TH2F* hempty_alpha = new TH2F("hempty_alpha", ";Smear parameter #alpha;#chi^{2} probability", 10, 0.9, 2.3, 10, 0.5, 1.0);
  xjjroot::sethempty(hempty_alpha, 0, 0.2);
  TH2F* hempty_fprompt = new TH2F("hempty_fprompt", ";f_{prompt};#chi^{2} probability", 10, 0.3, 1.0, 10, 0.5, 1.0);
  xjjroot::sethempty(hempty_fprompt, 0, 0.2);

  TLegend* leg = new TLegend(0.50, 0.22, 0.80, 0.22+0.035*2);
  xjjroot::setleg(leg, 0.03);
  leg->AddEntry(gr_prob_alpha_a, fitX::title_a.c_str(), "pe");
  leg->AddEntry(gr_prob_alpha_b, fitX::title_b.c_str(), "pe");

  xjjroot::setgstyle(1);
  TCanvas* c_alpha = new TCanvas("c_alpha", "", 600, 600);
  hempty_alpha->Draw();
  gr_prob_alpha_a->Draw("pl same");
  gr_prob_alpha_b->Draw("pl same");
  leg->Draw();
  xjjroot::drawCMS();
  std::string outputname_alpha(Form("plots/%s/cfitsmear_alpha.pdf", output.c_str()));
  xjjroot::mkdir(outputname_alpha);
  c_alpha->SaveAs(outputname_alpha.c_str());  

  TCanvas* c_fprompt = new TCanvas("c_fprompt", "", 600, 600);
  hempty_fprompt->Draw();
  gr_prob_fprompt_a->Draw("pl same");
  gr_prob_fprompt_b->Draw("pl same");
  leg->Draw();
  xjjroot::drawCMS();
  std::string outputname_fprompt(Form("plots/%s/cfitsmear_fprompt.pdf", output.c_str()));
  xjjroot::mkdir(outputname_fprompt);
  c_fprompt->SaveAs(outputname_fprompt.c_str());  
}

int main(int argc, char* argv[])
{
  if(argc==6) { lxyfit_fitlxy(argv[1], argv[2], argv[3], atoi(argv[4]), atof(argv[5])); return 0; }
  return 1;
}
