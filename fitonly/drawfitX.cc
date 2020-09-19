#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
// https://root.cern.ch/doc/master/rs102__hypotestwithshapes_8C.html
// https://root.cern.ch/doc/master/StandardProfileLikelihoodDemo_8C.html
#include <RooMinuit.h>
#include <RooStats/ProfileLikelihoodCalculator.h> 
#include <RooStats/HypoTestResult.h>
#include <RooStats/HypoTestPlot.h>
#include <RooStats/ConfInterval.h>
#include <RooStats/LikelihoodInterval.h>
#include <RooStats/LikelihoodIntervalPlot.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "fit.h"

float runfit(TH1F* h, TH1F* hmcp_a, TH1F* hmcp_b,                         // binned fits
             RooDataSet* dsh, RooDataSet* dshmcp_a, RooDataSet* dshmcp_b,  // unbinned fits
             std::string fitopt, float mm, TH1F* hyield,                   // fitting details
             std::string vname, std::string vtitle, std::string outputdir); // labels

void drawfitX(std::string input, std::string output, std::string fitopt="")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if(fitopt!="") { output += ("/"+fitopt); }
  if(fitopt=="range") { fitX::NBIN = 30; fitX::BIN_MAX = 3.92; }
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  RooDataSet* dsh = (RooDataSet*)ww->data("dsh");
  RooDataSet* dshBenr = (RooDataSet*)ww->data("dshBenr");
  RooDataSet* dshmcp_a = (RooDataSet*)ww->data("dshmcp_a");
  RooDataSet* dshmcp_b = (RooDataSet*)ww->data("dshmcp_b");
  TH1F* h = (TH1F*)inf->Get("h");
  TH1F* hBenr = (TH1F*)inf->Get("hBenr");
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");

  // fit + yield
  TH1F* hyield = new TH1F("hyield", ";;Raw Yield", 5, 0, 5);
  xjjroot::setthgrstyle(hyield, kBlack, 21, 1.2, kBlack, 2, 3, kBlack, 0.1, 1001);
  fitX::setaxis(hyield);
  TH1F* hBenryield = new TH1F("hBenryield", "Raw Yield", 5, 0, 5);
  xjjroot::setthgrstyle(hBenryield, kBlack, 25, 1.2, kBlack, 2, 3, kBlack, 0.1, 1001);
  fitX::setaxis(hBenryield);

  float mm = 0;
  xjjroot::setgstyle();
  // inclusive --> 
  mm = runfit(h, hmcp_a, hmcp_b,
              dsh, dshmcp_a, dshmcp_b,
              fitopt, mm, hyield,
              "th", Form("#splitline{Inclusive}{%s}", fitopt.c_str()),  Form("%s", output.c_str()));

  // // B-enr --> 
  // mm = runfit(hBenr, hmcp_a, hmcp_b,
  //             dshBenr, dshmcp_a, dshmcp_b,
  //             fitopt, mm, hBenryield,
  //             "thBenr", Form("#splitline{B-enriched (l_{xy} > 0.1 mm)}{%s}", fitopt.c_str()),  Form("plots/%s", output.c_str()));

}

int main(int argc, char* argv[])
{
  if(argc==4) { drawfitX(argv[1], argv[2], argv[3]); return 0; }
  if(argc==3) { drawfitX(argv[1], argv[2]); return 0; }
  return 1;
}

float runfit(TH1F* h, TH1F* hmcp_a, TH1F* hmcp_b,                         // binned fits
             RooDataSet* dsh, RooDataSet* dshmcp_a, RooDataSet* dshmcp_b,  // unbinned fits
             std::string fitopt, float mm, TH1F* hyield,                   // fitting details
             std::string vname, std::string vtitle, std::string outputdir) // labels
{
  // ====>
  std::map<std::string, fitX::fitXresult*> result = fitX::fit(h, 0, hmcp_a, hmcp_b, 
                                                              dsh, dshmcp_a, dshmcp_b,
                                                              std::string("plots/"+outputdir).c_str(), mm, true, "_"+vname, vtitle, fitopt);
  // <====

  // xjjroot::setgstyle();
  float ysig_a = result["unbinned"]->ysig_a();
  float ysigerr_a = result["unbinned"]->ysigerr_a();
  float ysig_b = result["unbinned"]->ysig_b();
  float ysigerr_b = result["unbinned"]->ysigerr_b();

  // yield
  hyield->SetBinContent(fitX::ibin_a, ysig_a);
  hyield->SetBinError(fitX::ibin_a, ysigerr_a);
  hyield->SetBinContent(fitX::ibin_b, ysig_b);
  hyield->SetBinError(fitX::ibin_b, ysigerr_b);

  std::cout<<"\e[36;1m";

  RooWorkspace* ww = result["unbinned"]->ww();
  RooRealVar* par10 = ww->var("par10");
  RooAddPdf* pdf = (RooAddPdf*)ww->obj("pdf");
  RooAbsData* data = ww->data("dsh_ws");
  RooRealVar* ndof = ww->var("ndof");

  std::string rootoutputname = "rootfiles/" + outputdir + "/rootnll.root";
  xjjroot::mkdir(rootoutputname);
  TFile* outf = new TFile(rootoutputname.c_str(), "recreate");
  outf->cd();

  // xjjroot::setgstyle();
  // TCanvas* cnll = new TCanvas("cnll", "", 600, 600);
  // RooAbsReal* nll = pdf->createNLL(*data); // <---------- switch
  // RooAbsReal* pll = nll->createProfile(*par10);

  // RooPlot* frame1 = par10->frame(RooFit::Bins(20), RooFit::Range(0., 200.));
  // xjjroot::sethempty(frame1);
  // nll->plotOn(frame1, RooFit::ShiftToZero(), RooFit::LineColor(kBlack), RooFit::LineWidth(4), RooFit::Normalization(2., RooAbsReal::Raw));
  // cnll->cd();
  // frame1->SetMinimum(0);
  // frame1->SetMaximum(20);
  // frame1->GetXaxis()->SetTitle("N_{sig}^{(X3872)}");
  // frame1->GetYaxis()->SetTitle("-2log#lambda");
  // frame1->Draw();
  // xjjroot::drawCMSleft("", 0.05, -0.08);
  // xjjroot::drawCMSright("1.7 nb^{-1} (2018 PbPb 5.02 TeV)");
  // std::string outputnamenll = "plots/" + outputdir + "/cnllscan_" + vname + ".pdf";
  // xjjroot::mkdir(outputnamenll);
  // cnll->SaveAs(outputnamenll.c_str());

  gDirectory->cd("root:/");
  RooWorkspace* wwnll = new RooWorkspace("wwnll");
  wwnll->import(*pdf);
  wwnll->import(*data);
  wwnll->import(*par10);
  wwnll->import(*ndof);

  // ----------> Likelihood scan using createNLL <----------
  // 
  std::cout<<"\e[36;2m";

  RooArgSet poi(*par10); // <---------- switch
  RooStats::ModelConfig* mcpdf = new RooStats::ModelConfig();
  mcpdf->SetWorkspace(*ww);
  mcpdf->SetPdf("pdf");
  mcpdf->SetParametersOfInterest(poi);
  mcpdf->SetSnapshot(poi);

  RooStats::ProfileLikelihoodCalculator* plc = new RooStats::ProfileLikelihoodCalculator();
  plc->SetData(*(ww->data("dsh_ws")));
  plc->SetModel(*mcpdf);

  // ----------> Null hypo test <----------

  RooArgSet* nullParams = (RooArgSet*)poi.snapshot();
  nullParams->setRealValue("par10", 0);
  plc->SetNullParameters(*nullParams);
  RooStats::HypoTestResult* htr = plc->GetHypoTest();
  // float nullpvalue = htr->NullPValue();
  // float significance = htr->Significance();
  std::cout<<"\e[0m"<<"\e[36;1m";
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "The p-value for the null is " << htr->NullPValue() << std::endl;
  std::cout << "Corresponding to a significance of " << htr->Significance() << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout<<"\e[0m";

  RooRealVar nullpvalue("nullpvalue", "", htr->NullPValue());
  RooRealVar nullsig("nullsig", "", htr->Significance());
  wwnll->import(nullpvalue);
  wwnll->import(nullsig);
  gDirectory->Add(wwnll);
  wwnll->Write();
  wwnll->Print();
  outf->Close();

  return 0;

  // ----------> Likelihood scan <----------

  // plc->SetConfidenceLevel(0.683);
  // RooStats::ConfInterval* interval1 = plc->GetInterval();
  // double lowerLimit1 = ((RooStats::LikelihoodInterval*)interval1)->LowerLimit(*par10);
  // double upperLimit1 = ((RooStats::LikelihoodInterval*)interval1)->UpperLimit(*par10);

  // plc->SetConfidenceLevel(0.95);
  // RooStats::ConfInterval* interval2 = plc->GetInterval();
  // double lowerLimit2 = ((RooStats::LikelihoodInterval*)interval2)->LowerLimit(*par10);
  // double upperLimit2 = ((RooStats::LikelihoodInterval*)interval2)->UpperLimit(*par10);
  // std::cout<<"\e[0m"<<"\e[36;1m";
  // std::cout << "-------------------------------------------------" << std::endl;
  // std::cout << "Profile lower limit (68%) on nsig = " << lowerLimit1 << std::endl;
  // std::cout << "Profile upper limit (68%) on nsig = " << upperLimit1 << std::endl;
  // std::cout << "Profile lower limit (95%) on nsig = " << lowerLimit2 << std::endl;
  // std::cout << "Profile upper limit (95%) on nsig = " << upperLimit2 << std::endl;
  // std::cout << "-------------------------------------------------\n" << std::endl;
  // std::cout<<"\e[0m"<<"\e[36;2m";

  // RooStats::LikelihoodIntervalPlot* plotInt = new RooStats::LikelihoodIntervalPlot((RooStats::LikelihoodInterval*)interval1);
  // plotInt->SetRange(0, 210);
  // plotInt->Draw("tf1");
  // TH1D* hplotInt = (TH1D*)((TH1D*)plotInt->GetPlottedObject())->Clone("hplotInt");
  // for(int i=0;i<hplotInt->GetXaxis()->GetNbins();i++) { hplotInt->SetBinContent(i+1, hplotInt->GetBinContent(i+1)*2); }
  // hplotInt->SetTitle(";N_{sig}^{(X3872)};-2log#lambda(N_{sig}^{(X3872)})");
  // xjjroot::sethempty(hplotInt, 0, 0);
  // hplotInt->SetMinimum(0.);
  // hplotInt->SetMaximum(20.);
  // float hxmin = hplotInt->GetXaxis()->GetXmin();
  // float hxmax = hplotInt->GetXaxis()->GetXmax();
  // hplotInt->SetLineColor(kBlack);
  // hplotInt->SetLineWidth(3);
  // xjjroot::setgstyle();
  // TCanvas* cll = new TCanvas("cll", "", 600, 600);
  // hplotInt->Draw("L");
  // xjjroot::drawline(hxmin, 0.50*2, hxmax, 0.50*2, xjjroot::mycolor_middle["red"], 2, 3);
  // xjjroot::drawline(hxmin, 1.96*2, hxmax, 1.96*2, xjjroot::mycolor_middle["red"], 2, 3);
  // xjjroot::drawline(lowerLimit1, 0, lowerLimit1, 0.50*2, xjjroot::mycolor_middle["red"], 2, 3);
  // xjjroot::drawline(upperLimit1, 0, upperLimit1, 0.50*2, xjjroot::mycolor_middle["red"], 2, 3);
  // xjjroot::drawline(lowerLimit2, 0, lowerLimit2, 1.96*2, xjjroot::mycolor_middle["red"], 2, 3);
  // xjjroot::drawline(upperLimit2, 0, upperLimit2, 1.96*2, xjjroot::mycolor_middle["red"], 2, 3);
  // xjjroot::drawCMSleft("", 0.05, -0.08);
  // xjjroot::drawCMSright("1.7 nb^{-1} (2018 PbPb 5.02 TeV)");

  // std::string outputname = "plots/" + outputdir + "/cllscan_" + vname + ".pdf";
  // xjjroot::mkdir(outputname);

  // // ww->Print();

  // std::cout<<"\e[0m"<<std::endl;
  // cll->SaveAs(outputname.c_str());
  // std::cout<<std::endl;
  // return result["unbinned"]->msig_b();
}


// hyield->SetMaximum(hyield->GetMaximum()*2.5);
// hyield->Draw("ple");
// fitX::drawkinematics();
// xjjroot::drawtex(0.22, 0.84, Form("#splitline{Inclusive}{%s}", fitopt.c_str()), 0.042, 12, 62);
// xjjroot::drawCMS();

// auto data = ww->data("dsh_ws");
// RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*)ww->obj("ModelConfig");
// sbModel->SetName("S+B Model");
// RooRealVar* poi = (RooRealVar*)sbModel->GetParametersOfInterest()->first();
// poi->setVal(50);
// sbModel->SetSnapshot(RooArgSet(*poi));
// auto bModel = sbModel->Clone();
// bModel->SetName("B Model");
// poi->setVal(0);
// bModel->SetSnapshot(RooArgSet(*poi));

// std::cout<<"\e[0m"<<"\e[36;2m";

// RooStats::ProfileLikelihoodCalculator plc;
// plc.SetData(*(ww->data("dsh_ws")));
// RooArgSet* nullParams = (RooArgSet*)poi.snapshot();
// nullParams->setRealValue("par10", 0);
// plc.SetModel(*mcpdf);
// plc.SetNullParameters(*nullParams);

// RooStats::HypoTestResult* htr = plc.GetHypoTest();

