#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooBinning.h>

#include <string>
#include <vector>

#include "var.h"
#include "fit.h"
// #include "fit_a.h"

#include "xjjcuti.h"
#include "xjjrootuti.h"

void drawkinematics();
void drawdatamc(std::string input, std::string output, std::string type)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  output += (fitX::tagname()+"/"+type);

  datamc::var* vv = new datamc::var(type.c_str());
  if(!vv->valid()) return;
  int nn = vv->n()-1;
  std::vector<TH1F*> h(nn, 0);
  for(int i=0; i<nn; i++) { h[i] = (TH1F*)inf->Get(Form("h%d", i)); }
  std::vector<TF1*> ff(nn, 0);
  for(int i=0; i<nn; i++) { ff[i] = (TF1*)inf->Get(Form("f--%d-unbinned", i)); }
  std::vector<std::string> tt(nn);
  std::vector<Color_t> cc(nn);
  for(int i=0; i<nn; i++) 
    { int icut = vv->gt()?i:i+1;
      tt[i] = Form("%s < %s < %s %s",
                   xjjc::number_remove_zero(vv->gt()?vv->vars()[icut]:vv->vars().front()).c_str(),
                   vv->title().c_str(),
                   xjjc::number_remove_zero(vv->gt()?vv->vars().back():vv->vars()[icut]).c_str(),
                   vv->unit().c_str()); 
      cc[i] = h[i]->GetLineColor(); }

  // distribution
  TH1F* hmcdis_a = (TH1F*)inf->Get("hmcdis_a");
  TH1F* hmcdis_b = (TH1F*)inf->Get("hmcdis_b");
  TH1F* gmcdis_a = (TH1F*)inf->Get("gmcdis_a");
  TH1F* gmcdis_b = (TH1F*)inf->Get("gmcdis_b");
  TH1F* hbkgdis_a = (TH1F*)inf->Get("hbkgdis_a");
  TH1F* hbkgdis_b = (TH1F*)inf->Get("hbkgdis_b");
  TH1F* hdis_a = (TH1F*)inf->Get("hdis_a");
  TH1F* hdis_b = (TH1F*)inf->Get("hdis_b");
  TH1F* gdis_a = (TH1F*)inf->Get("gdis_a");
  TH1F* gdis_b = (TH1F*)inf->Get("gdis_b");
  TH1F* hratiodis_a = (TH1F*)inf->Get("hratiodis_a");
  TH1F* hratiodis_b = (TH1F*)inf->Get("hratiodis_b");
  TH1F* gratiodis_a = (TH1F*)inf->Get("gratiodis_a");
  TH1F* gratiodis_b = (TH1F*)inf->Get("gratiodis_b");

  // chi2 test
  std::map<std::string, double> chi2_a = xjjroot::chi2test(hdis_a, hmcdis_a, "UW");
  std::map<std::string, double> chi2_b = xjjroot::chi2test(hdis_b, hmcdis_b, "UW");

  // Draw
  TLegend* leg_a = new TLegend(0.22, 0.81-0.047*3, 0.64, 0.81);
  xjjroot::setleg(leg_a, 0.042);
  leg_a->AddEntry(gdis_a, "Data signal", "pl");
  leg_a->AddEntry(hmcdis_a, "MC", "pl");
  leg_a->AddEntry(hbkgdis_a, "Sideband", "pl");
  TLegend* leg_b = new TLegend(0.22, 0.81-0.047*3, 0.64, 0.81);
  xjjroot::setleg(leg_b, 0.042);
  leg_b->AddEntry(gdis_b, "Data signal", "pl");
  leg_b->AddEntry(hmcdis_b, "MC", "pl");
  leg_b->AddEntry(hbkgdis_b, "Sideband", "pl");

  float ymax = std::max(h.front()->GetMaximum(), h.back()->GetMaximum())*1.2;
  TH2F* hempty_a = new TH2F("hempty_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});%s", Form("Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3)), 
                            fitX::NBIN/2, fitX::BIN_MIN, fitX::BIN_MIN+fitX::BIN_WIDTH*fitX::NBIN/2, 10, 0, ymax);
  // fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX, 10, 0, ymax);
  hempty_a->GetXaxis()->SetNdivisions(505);
  xjjroot::sethempty(hempty_a, 0, 0);
  TH2F* hempty_b = new TH2F("hempty_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});%s", Form("Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3)), 
                            fitX::NBIN/2, fitX::BIN_MAX-fitX::BIN_WIDTH*fitX::NBIN/2, fitX::BIN_MAX, 10, 0, ymax);
  hempty_b->GetXaxis()->SetNdivisions(505);
  xjjroot::sethempty(hempty_b, 0, 0);
  
  xjjroot::setgstyle(1);

  TCanvas* c_a = new TCanvas("c_a", "", 1800, 600);
  c_a->Divide(3, 1);
  c_a->cd(1);
  hempty_a->Draw();
  for(int i=0;i<nn;i++) { h[i]->Draw("pe same"); ff[i]->Draw("same"); }
  xjjroot::drawtexgroup(0.89, 0.86, tt, 1, 0.5, 0.038, 33, 62, cc);
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.038, 12, 62, kBlack);
  xjjroot::drawCMS();
  xjjroot::drawcomment(output);
  c_a->cd(2);
  hmcdis_a->Draw("hist e");
  hbkgdis_a->Draw("hist e same");
  gdis_a->Draw("pe same");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.038, 12, 62, fitX::color_a);
  leg_a->Draw();
  xjjroot::drawCMS();
  c_a->cd(3);
  hratiodis_a->Draw("AXIS");
  xjjroot::drawline(vv->vars().front(), 1., vv->vars().back(), 1., fitX::color_a, 1, 2, 0.6);
  gratiodis_a->Draw("pe same");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.055, Form("#chi^{2} prob = %.1f%s", chi2_a["chi2prob"]*1.e+2, "%"), 0.04, 12, 42, kBlack);
  xjjroot::drawCMS();
  std::string outputname_a(Form("plots/%s/cdis_a.pdf", output.c_str()));
  xjjroot::mkdir(outputname_a);
  c_a->SaveAs(outputname_a.c_str());

  TCanvas* c_b = new TCanvas("c_b", "", 1800, 600);
  c_b->Divide(3, 1);
  c_b->cd(1);
  hempty_b->Draw();
  for(int i=0;i<nn;i++) { h[i]->Draw("pe same"); ff[i]->Draw("same"); }
  xjjroot::drawtexgroup(0.89, 0.86, tt, 1, 0.5, 0.038, 33, 62, cc);
  xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.038, 12, 62, kBlack);
  xjjroot::drawCMS();
  xjjroot::drawcomment(output);
  c_b->cd(2);
  hmcdis_b->Draw("hist e");
  hbkgdis_b->Draw("hist e same");
  gdis_b->Draw("pe same");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.038, 12, 62, fitX::color_b);
  leg_b->Draw();
  xjjroot::drawCMS();
  c_b->cd(3);
  hratiodis_b->Draw("AXIS");
  xjjroot::drawline(vv->vars().front(), 1., vv->vars().back(), 1., fitX::color_b, 1, 2, 0.6);
  gratiodis_b->Draw("pe same");
  drawkinematics();
  xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.038, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.24, 0.84-0.055, Form("#chi^{2} prob = %.1f%s", chi2_b["chi2prob"]*1.e+2, "%"), 0.04, 12, 42, kBlack);
  xjjroot::drawCMS();
  std::string outputname_b(Form("plots/%s/cdis_b.pdf", output.c_str()));
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());

}

int main(int argc, char* argv[])
{
  if(argc==4) { drawdatamc(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}

void drawkinematics()
{
  xjjroot::drawtex(0.90, 0.84, fitX::pttag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 62);
}

