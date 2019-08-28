#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TCut.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

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
  TH1F* hdis_a = new TH1F("hdis_a", Form(";%s;Probability", vv->title().c_str()), vv->n()-1, vv->vars().data());
  hdis_a->GetXaxis()->SetNdivisions(505);
  TH1F* hdis_b = new TH1F("hdis_b", Form(";%s;", vv->title().c_str()), vv->n()-1, vv->vars().data());
  hdis_b->GetXaxis()->SetNdivisions(505);
  TH1F* hdisratio_a = (TH1F*)hdis_a->Clone("hdisratio_a");
  hdisratio_a->Divide(hmcdis_a);

  // fit
  std::vector<float> ysig_a(vv->n()-1), ysigerr_a(vv->n()-1), ysig_b(vv->n()-1), ysigerr_b(vv->n()-1);
  for(int i=0;i<vv->n()-1;i++)
    {
      std::map<std::string, fitX::fitXresult*> result = fitX::fit(h[i], 0, hmcp_a, hmcp_b,
                                                                  dsh[i], dshmcp_a, dshmcp_b,
                                                                  Form("plots/%s/idx", output.c_str()), false, true, Form("-%d", i), "default"); // fix mean = false
      ysig_a[i] = result["unbinned"]->ysig_a();
      ysigerr_a[i] = result["unbinned"]->ysigerr_a();
      ysig_b[i] = result["unbinned"]->ysig_b();
      ysigerr_b[i] = result["unbinned"]->ysigerr_b();
    }
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
          // std::cout<<thisi<<" "<<ysig_a[thisi]<<" "<<ysigerr_a[thisi]<<" "<<lasti<<" "<<ysig_a[lasti]<<" "<<ysigerr_a[lasti]<<std::endl;
          // std::cout<<thisi<<" "<<ysig_b[thisi]<<" "<<ysigerr_b[thisi]<<" "<<lasti<<" "<<ysig_b[lasti]<<" "<<ysigerr_b[lasti]<<std::endl;
          // std::cout<<yysig_a<<" "<<yysigerr_a<<std::endl;
          // std::cout<<yysig_b<<" "<<yysigerr_b<<std::endl;
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
  hdis_a->Scale(1./hdis_a->Integral(), "width");
  hdis_b->Scale(1./hdis_b->Integral(), "width");

  hmcdis_a->SetMinimum(0);
  hmcdis_a->SetMaximum(std::max(hmcdis_a->GetMaximum(), hdis_a->GetMaximum())*1.5);
  xjjroot::sethempty(hmcdis_a, 0, 0);
  xjjroot::setthgrstyle(hmcdis_a, fitX::color_a, 21, 1., fitX::color_a, 1, 2);
  xjjroot::setthgrstyle(hdis_a, kBlack, 20, 1., kBlack, 1, 2);

  hmcdis_b->SetMinimum(0);
  hmcdis_b->SetMaximum(std::max(hmcdis_b->GetMaximum(), hdis_b->GetMaximum())*1.5);
  xjjroot::sethempty(hmcdis_b, 0, 0);
  xjjroot::setthgrstyle(hmcdis_b, fitX::color_b, 21, 1., fitX::color_b, 1, 2);
  xjjroot::setthgrstyle(hdis_b, kBlack, 20, 1., kBlack, 1, 2);

  xjjroot::setgstyle(1);
  TCanvas* c_a = new TCanvas("c_a", "", 1200, 600);
  c_a->Divide(2, 1);
  c_a->cd(1);
  hmcdis_a->Draw("hist e");
  hdis_a->Draw("pe same");
  drawkinematics();
  xjjroot::drawCMS();
  std::string outputname_a(Form("plots/%s/cdis_a.pdf", output.c_str(), vv->type().c_str()));
  xjjroot::mkdir(outputname_a);
  c_a->SaveAs(outputname_a.c_str());

  TCanvas* c_b = new TCanvas("c_b", "", 1200, 600);
  c_b->Divide(2, 1);
  c_b->cd(1);
  hmcdis_b->Draw("hist e");
  hdis_b->Draw("pe same");
  drawkinematics();
  xjjroot::drawCMS();
  std::string outputname_b(Form("plots/%s/cdis_b.pdf", output.c_str(), vv->type().c_str()));
  xjjroot::mkdir(outputname_b);
  c_b->SaveAs(outputname_b.c_str());
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
