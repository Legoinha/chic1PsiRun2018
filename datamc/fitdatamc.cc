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

#include "xjjcuti.h"
#include "xjjrootuti.h"

void drawkinematics();
void fitdatamc(std::string input, std::string output, std::string type)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  output += (fitX::tagname()+"/"+type);

  datamc::var* vv = new datamc::var(type.c_str());
  if(!vv->valid()) return;
  // mass
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  std::vector<RooDataSet*> dsh(vv->n()-1);
  for(int i=0; i<vv->n()-1; i++) { dsh[i] = (RooDataSet*)ww->data(Form("dsh%d", i)); }
  RooDataSet* dshmcp_a = (RooDataSet*)ww->data("dshmcp_a");
  RooDataSet* dshmcp_b = (RooDataSet*)ww->data("dshmcp_b");
  std::vector<TH1F*> h(vv->n()-1);
  for(int i=0; i<vv->n()-1; i++) { h[i] = (TH1F*)inf->Get(Form("h%d", i)); }
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  hmcp_a->Scale(hmcp_a->GetEntries()/hmcp_a->Integral());
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
  hmcp_b->Scale(hmcp_b->GetEntries()/hmcp_b->Integral());

  // distribution
  TH1F* hmcdis_a = (TH1F*)inf->Get("hmcdis_a");
  hmcdis_a->GetXaxis()->SetNdivisions(505);
  TH1F* hmcdis_b = (TH1F*)inf->Get("hmcdis_b");
  hmcdis_b->GetXaxis()->SetNdivisions(505);
  TH1F* hmcfinedis_a = (TH1F*)inf->Get("hmcfinedis_a");
  hmcfinedis_a->GetXaxis()->SetNdivisions(505);
  TH1F* hmcfinedis_b = (TH1F*)inf->Get("hmcfinedis_b");
  hmcfinedis_b->GetXaxis()->SetNdivisions(505);
  TH1F* hbkgdis_a = (TH1F*)inf->Get("hbkgdis_a");
  hbkgdis_a->GetXaxis()->SetNdivisions(505);
  TH1F* hbkgdis_b = (TH1F*)inf->Get("hbkgdis_b");
  hbkgdis_b->GetXaxis()->SetNdivisions(505);
  TH1F* hdis_a = new TH1F("hdis_a", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());
  hdis_a->GetXaxis()->SetNdivisions(505);
  TH1F* hdis_b = new TH1F("hdis_b", Form(";%s %s;Probability", vv->title().c_str(), vv->unit().c_str()), vv->n()-1, vv->vars().data());
  hdis_b->GetXaxis()->SetNdivisions(505);
  RooRealVar* varr = new RooRealVar(vv->formula().c_str(), "varr", vv->vars().front(), vv->vars().back());

  // fit
  std::vector<TF1*> ff(vv->n()-1, 0);
  std::vector<std::string> tt(vv->n()-1);
  std::vector<Color_t> cc(vv->n()-1);
  std::vector<float> ysig_a(vv->n()-1), ysigerr_a(vv->n()-1), ysig_b(vv->n()-1), ysigerr_b(vv->n()-1);
  RooDataSet* dsh_ws;
  RooWorkspace* wws = new RooWorkspace("wws");
  for(int i=0;i<vv->n()-1;i++)
    {
      int icut = vv->gt()?i:i+1;
      std::string label(Form("%s < %s < %s %s", 
                             xjjc::number_remove_zero(vv->gt()?vv->vars()[icut]:vv->vars().front()).c_str(), 
                             vv->title().c_str(),
                             xjjc::number_remove_zero(vv->gt()?vv->vars().back():vv->vars()[icut]).c_str(),
                             vv->unit().c_str()));
      std::map<std::string, fitX::fitXresult*> result = fitX::fit(h[i], 0, hmcp_a, hmcp_b,
                                                                  dsh[i], dshmcp_a, dshmcp_b,
                                                                  Form("plots/%s/idx", output.c_str()), 0, true, Form("-%d", i), label); // fix mean = false
      ysig_a[i] = result["unbinned"]->ysig_a();
      ysigerr_a[i] = result["unbinned"]->ysigerr_a();
      ysig_b[i] = result["unbinned"]->ysig_b();
      ysigerr_b[i] = result["unbinned"]->ysigerr_b();
      ff[i] = result["unbinned"]->f();
      xjjroot::setthgrstyle(h[i], xjjroot::mycolor_middle[xjjroot::cc[i]], 20, 0.9, xjjroot::mycolor_middle[xjjroot::cc[i]], 1, 1);
      xjjroot::settfstyle(ff[i], xjjroot::mycolor_middle[xjjroot::cc[i]], 7, 2);
      tt[i] = label;
      cc[i] = xjjroot::mycolor_middle[xjjroot::cc[i]];
      if((vv->gt() && i==0) || (!vv->gt() && i==(vv->n()-2)))
        {
          RooWorkspace* w = result["unbinned"]->ww();
          dsh_ws = (RooDataSet*)w->data(Form("%s_ws", dsh[i]->GetName()));
        }
    }
  // RooArgSet* variables = roopdf->getVariables();
  // RooArgSet* params = roopdf->getParameters((*variables)["Bmass"]);
  // std::cout<<((RooRealVar*)(params->find("par5")))->getVal()<<" "<<((RooRealVar*)(params->find("par10")))->getVal()<<" "<<((RooRealVar*)(params->find("nbkg")))->getVal()<<std::endl;
  // RooStats::SPlot* sData = new RooStats::SPlot("sData", "An SPlot", *roodsh, roopdf, RooArgList(*(params->find("par5")), *(params->find("par10")), *(params->find("nbkg"))));
  // RooStats::SPlot* sData = new RooStats::SPlot("sData", "An SPlot", *roodsh, roopdf, RooArgList(*(params->find("par5")), *(params->find("nbkg"))));
  wws->import(*dsh_ws);
  dsh_ws->Print("v");
  RooDataSet* dsh_ws_par5 = new RooDataSet(dsh_ws->GetName(), dsh_ws->GetTitle(), dsh_ws, *dsh_ws->get(), 0, "par5_sw");
  RooDataSet* dsh_ws_par10 = new RooDataSet(dsh_ws->GetName(), dsh_ws->GetTitle(), dsh_ws, *dsh_ws->get(), 0, "par10_sw");

  for(int i=0; i<vv->n()-1; i++)
    {
      float yysig_a, yysig_b, yysigerr_a, yysigerr_b;
      int thisi = vv->gt()?(vv->n()-2-i):i;
      if(i)
        {
          int lasti = vv->gt()?(thisi+1):(thisi-1);
          yysig_a = std::max(ysig_a[thisi] - ysig_a[lasti], (float)0.);
          yysig_b = std::max(ysig_b[thisi] - ysig_b[lasti], (float)0.);
          yysigerr_a = TMath::Sqrt(TMath::Abs(ysigerr_a[thisi]*ysigerr_a[thisi] - ysigerr_a[lasti]*ysigerr_a[lasti]));
          yysigerr_b = TMath::Sqrt(TMath::Abs(ysigerr_b[thisi]*ysigerr_b[thisi] - ysigerr_b[lasti]*ysigerr_b[lasti]));
        }
      else
        {
          yysig_a = ysig_a[thisi];
          yysig_b = ysig_b[thisi];
          yysigerr_a = ysigerr_a[thisi];
          yysigerr_b = ysigerr_b[thisi];
        }
      hdis_a->SetBinContent(thisi+1, yysig_a);
      hdis_a->SetBinError(thisi+1, yysigerr_a);
      hdis_b->SetBinContent(thisi+1, yysig_b);
      hdis_b->SetBinError(thisi+1, yysigerr_b);
    }
  hmcdis_a->Scale(1./hmcdis_a->Integral(), "width");
  hmcdis_b->Scale(1./hmcdis_b->Integral(), "width");
  hmcfinedis_a->Scale(1./hmcfinedis_a->Integral(), "width");
  hmcfinedis_b->Scale(1./hmcfinedis_b->Integral(), "width");
  hbkgdis_a->Scale(1./hbkgdis_a->Integral(), "width");
  hbkgdis_b->Scale(1./hbkgdis_b->Integral(), "width");
  hdis_a->Scale(1./hdis_a->Integral(), "width");
  hdis_b->Scale(1./hdis_b->Integral(), "width");

  // Get bin weighted center
  std::vector<double> xw_a(vv->n()-1), yw_a(vv->n()-1), xel_a(vv->n()-1), xeh_a(vv->n()-1), ye_a(vv->n()-1);
  std::vector<double> xw_b(vv->n()-1), yw_b(vv->n()-1), xel_b(vv->n()-1), xeh_b(vv->n()-1), ye_b(vv->n()-1);
  for(int i=0; i<vv->n()-1; i++)
    {
      int nbinl = hmcfinedis_a->FindBin(vv->vars()[i]);
      int nbinh = hmcfinedis_a->FindBin(vv->vars()[i+1])-1;
      yw_a[i] = hmcfinedis_a->IntegralAndError(nbinl, nbinh, ye_a[i], "width");
      yw_b[i] = hmcfinedis_b->IntegralAndError(nbinl, nbinh, ye_b[i], "width");
      xw_a[i] = 0;
      xw_b[i] = 0;
      for(int ibin=nbinl; ibin<=nbinh; ibin++)
        {
          xw_a[i] += (hmcfinedis_a->GetBinContent(ibin)*hmcfinedis_a->GetBinCenter(ibin)*hmcfinedis_a->GetBinWidth(ibin));
          xw_b[i] += (hmcfinedis_b->GetBinContent(ibin)*hmcfinedis_b->GetBinCenter(ibin)*hmcfinedis_b->GetBinWidth(ibin));
        }
      xw_a[i] /= yw_a[i];
      xw_b[i] /= yw_b[i];
      yw_a[i] /= hmcdis_a->GetBinWidth(i+1);
      yw_b[i] /= hmcdis_b->GetBinWidth(i+1);
      ye_a[i] /= hmcdis_a->GetBinWidth(i+1);
      ye_b[i] /= hmcdis_b->GetBinWidth(i+1);
      xel_a[i] = xw_a[i] - vv->vars()[i];
      xeh_a[i] = vv->vars()[i+1] - xw_a[i];
      xel_b[i] = xw_b[i] - vv->vars()[i];
      xeh_b[i] = vv->vars()[i+1] - xw_b[i];
    }
  TGraphAsymmErrors* gmcdis_a = new TGraphAsymmErrors(vv->n()-1, xw_a.data(), yw_a.data(), xel_a.data(), xeh_a.data(), ye_a.data(), ye_a.data());
  gmcdis_a->SetName("gmcdis_a");
  TGraphAsymmErrors* gmcdis_b = new TGraphAsymmErrors(vv->n()-1, xw_b.data(), yw_b.data(), xel_b.data(), xeh_b.data(), ye_b.data(), ye_b.data());
  gmcdis_b->SetName("gmcdis_b");

  TH1F* hratiodis_a = (TH1F*)hdis_a->Clone("hratiodis_a");
  hratiodis_a->GetYaxis()->SetTitle("Data / MC");
  hratiodis_a->Divide(hmcdis_a);
  hratiodis_a->SetMinimum(0);
  hratiodis_a->SetMaximum(3);
  TH1F* hratiodis_b = (TH1F*)hdis_b->Clone("hratiodis_b");
  hratiodis_b->GetYaxis()->SetTitle("Data / MC");
  hratiodis_b->Divide(hmcdis_b);
  hratiodis_b->SetMinimum(0);
  hratiodis_b->SetMaximum(3);

  TGraphAsymmErrors* gdis_a = xjjroot::setwcenter(hdis_a, xw_a, "gdis_a");
  TGraphAsymmErrors* gdis_b = xjjroot::setwcenter(hdis_b, xw_b, "gdis_b");
  TGraphAsymmErrors* gratiodis_a = xjjroot::setwcenter(hratiodis_a, xw_a, "gratiodis_a");
  TGraphAsymmErrors* gratiodis_b = xjjroot::setwcenter(hratiodis_b, xw_b, "gratiodis_b");

  // Draw
  hmcdis_a->SetMinimum(0);
  hmcdis_a->SetMaximum(std::max(hmcdis_a->GetMaximum(), hdis_a->GetMaximum())*2);
  hmcfinedis_a->SetMinimum(0);
  hmcfinedis_a->SetMaximum(hmcfinedis_a->GetMaximum()*1.2);
  xjjroot::sethempty(hmcdis_a, 0, 0);
  xjjroot::setthgrstyle(hmcdis_a, fitX::color_a, 21, 0.5, fitX::color_a, 1, 1, fitX::color_a, 0.6, 3005);
  xjjroot::setthgrstyle(gmcdis_a, fitX::color_a, 41, 1.9, fitX::color_a, 1, 1);
  xjjroot::sethempty(hmcfinedis_a, 0, 0);
  xjjroot::setthgrstyle(hmcfinedis_a, fitX::color_a, 0, 0., fitX::color_a, 3, 1, fitX::color_a, 0.3, 1001, 0.3);
  xjjroot::sethempty(hbkgdis_a, 0, 0);
  xjjroot::setthgrstyle(hbkgdis_a, kGray+1, 20, 1.1, kGray+1, 2, 1);
  xjjroot::sethempty(hdis_a, 0, 0);
  xjjroot::setthgrstyle(hdis_a, kBlack, 47, 1.9, kBlack, 1, 1);
  xjjroot::setthgrstyle(gdis_a, kBlack, 47, 1.9, kBlack, 1, 1);
  xjjroot::sethempty(hratiodis_a, 0, 0);
  xjjroot::setthgrstyle(hratiodis_a, kBlack, 47, 1.9, kBlack, 1, 1);
  xjjroot::setthgrstyle(gratiodis_a, kBlack, 47, 1.9, kBlack, 1, 1);
  hmcdis_b->SetMinimum(0);
  hmcdis_b->SetMaximum(std::max(hmcdis_b->GetMaximum(), hdis_b->GetMaximum())*2);
  hmcfinedis_b->SetMinimum(0);
  hmcfinedis_b->SetMaximum(hmcfinedis_b->GetMaximum()*1.2);
  xjjroot::sethempty(hmcdis_b, 0, 0);
  xjjroot::setthgrstyle(hmcdis_b, fitX::color_b, 21, 0.5, fitX::color_b, 1, 1, fitX::color_b, 0.6, 3005);
  xjjroot::setthgrstyle(gmcdis_b, fitX::color_b, 41, 1.9, fitX::color_b, 1, 1);
  xjjroot::sethempty(hmcfinedis_b, 0, 0);
  xjjroot::setthgrstyle(hmcfinedis_b, fitX::color_b, 0, 0., fitX::color_b, 3, 1, fitX::color_b, 0.3, 1001, 0.3);
  xjjroot::sethempty(hbkgdis_b, 0, 0);
  xjjroot::setthgrstyle(hbkgdis_b, kGray+1, 20, 1.1, kGray+1, 2, 1);
  xjjroot::sethempty(hdis_b, 0, 0);
  xjjroot::setthgrstyle(hdis_b, kBlack, 47, 1.9, kBlack, 1, 1);
  xjjroot::setthgrstyle(gdis_b, kBlack, 47, 1.9, kBlack, 1, 1);
  xjjroot::sethempty(hratiodis_b, 0, 0);
  xjjroot::setthgrstyle(hratiodis_b, kBlack, 47, 1.9, kBlack, 1, 1);
  xjjroot::setthgrstyle(gratiodis_b, kBlack, 47, 1.9, kBlack, 1, 1);

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
  hempty_a->GetXaxis()->SetNdivisions(505);
  xjjroot::sethempty(hempty_a, 0, 0);
  TH2F* hempty_b = new TH2F("hempty_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});%s", Form("Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3)), 
                            fitX::NBIN/2, fitX::BIN_MAX-fitX::BIN_WIDTH*fitX::NBIN/2, fitX::BIN_MAX, 10, 0, ymax);
  hempty_b->GetXaxis()->SetNdivisions(505);
  xjjroot::sethempty(hempty_b, 0, 0);
  
  xjjroot::setgstyle(1);

  TCanvas* c_wc = new TCanvas("c_wc", "", 1200, 600);
  c_wc->Divide(2, 1);
  c_wc->cd(1);
  hmcfinedis_a->Draw("hist");
  hmcdis_a->Draw("hist same");
  gmcdis_a->Draw("pe same");
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.038, 12, 62, kBlack);
  drawkinematics();
  xjjroot::drawCMS("Simulation");
  c_wc->cd(2);
  hmcfinedis_b->Draw("hist");
  hmcdis_b->Draw("hist same");
  gmcdis_b->Draw("pe same");
  xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.038, 12, 62, kBlack);
  drawkinematics();
  xjjroot::drawCMS("Simulation");
  c_wc->SaveAs(Form("plots/%s/cdis_wc.pdf", output.c_str()));

  TCanvas* c_ws = new TCanvas("c_ws", "", 1200, 600);
  c_ws->Divide(2, 1);
  c_ws->cd(1);
  RooPlot* frame_a = varr->frame();
  dsh_ws_par5->plotOn(frame_a, RooFit::Binning(RooBinning(vv->n()-1, vv->vars().data())), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(21), RooFit::LineColor(1), RooFit::LineWidth(1));
  frame_a->Draw();
  c_ws->cd(2);
  RooPlot* frame_b = varr->frame();
  dsh_ws_par10->plotOn(frame_b, RooFit::Binning(RooBinning(vv->n()-1, vv->vars().data())), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(21), RooFit::LineColor(1), RooFit::LineWidth(1));
  frame_b->Draw();
  c_ws->SaveAs(Form("plots/%s/cdis_ws.pdf", output.c_str()));

  TCanvas* c_a = new TCanvas("c_a", "", 1800, 600);
  c_a->Divide(3, 1);
  c_a->cd(1);
  hempty_a->Draw();
  for(int i=0;i<vv->n()-1;i++) { h[i]->Draw("pe same"); ff[i]->Draw("same"); }
  xjjroot::drawtexgroup(0.89, 0.86, tt, 1, 0.5, 0.038, 33, 62, cc);
  xjjroot::drawtex(0.24, 0.84, fitX::title_a.c_str(), 0.038, 12, 62, kBlack);
  xjjroot::drawCMS();
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
  xjjroot::drawCMS();
  std::string outputname_a(Form("plots/%s/cdis_a.pdf", output.c_str()));
  xjjroot::mkdir(outputname_a);
  c_a->SaveAs(outputname_a.c_str());

  TCanvas* c_b = new TCanvas("c_b", "", 1800, 600);
  c_b->Divide(3, 1);
  c_b->cd(1);
  hempty_b->Draw();
  for(int i=0;i<vv->n()-1;i++) { h[i]->Draw("pe same"); ff[i]->Draw("same"); }
  xjjroot::drawtexgroup(0.89, 0.86, tt, 1, 0.5, 0.038, 33, 62, cc);
  xjjroot::drawtex(0.24, 0.84, fitX::title_b.c_str(), 0.038, 12, 62, kBlack);
  xjjroot::drawCMS();
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
  xjjroot::drawCMS();
  std::string outputname_b(Form("plots/%s/cdis_b.pdf", output.c_str()));
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());

  std::string outputname = std::string("rootfiles/"+output+"/datamc_fithist.root");
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hmcdis_a->Write();
  hmcdis_b->Write();
  gmcdis_a->Write();
  gmcdis_b->Write();
  hmcfinedis_a->Write();
  hmcfinedis_b->Write();
  hdis_a->Write();
  hdis_b->Write();
  gdis_a->Write();
  gdis_b->Write();
  hratiodis_a->Write();
  hratiodis_b->Write();
  gratiodis_a->Write();
  gratiodis_b->Write();
  outf->cd();
  gDirectory->Add(wws);
  wws->Write();
  wws->Print();
  outf->cd();
  fitX::write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==4) { fitdatamc(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}

void drawkinematics()
{
  xjjroot::drawtex(0.90, 0.84, fitX::pttag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04, fitX::ytag().c_str(), 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.04*2, fitX::centtag().c_str(), 0.038, 32, 62);
}

