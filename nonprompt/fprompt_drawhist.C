#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "lxydis.h"
#include "fitX.h"

std::vector<std::string> lxyvars = {"lxy", "lxyz"};
float fprompt_a = 0.388365, errfprompt_a = 0.0778002;
void fprompt_drawhist(std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::map<std::string, std::vector<TEfficiency*>> grfprompt_a, grfprompt_b;
  int nn = 0; for(auto& lxyvar : lxyvars) { nn += lxydis::lxycut[lxyvar].size(); }
  nn++;
  
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";;f_{prompt}", nn, 0, nn, 10, 0.1, 1.3);
  xjjroot::sethempty(hemptyfprompt, 0, 0.2);
  hemptyfprompt->GetXaxis()->SetLabelSize(hemptyfprompt->GetXaxis()->GetLabelSize()*1.2);
  float xmin = hemptyfprompt->GetXaxis()->GetXmin(), xmax = hemptyfprompt->GetXaxis()->GetXmax(),
    ymin = hemptyfprompt->GetYaxis()->GetXmin(), ymax = hemptyfprompt->GetYaxis()->GetXmax();

  std::vector<float> x(nn), ex(nn), y_a(nn), eyl_a(nn), eyh_a(nn), y_b(nn), eyl_b(nn), eyh_b(nn);
  for(int j=0, n=0; j<lxyvars.size(); j++)
    {
      std::string lxyvar = lxyvars[j];
      // -->
      std::cout<<("rootfiles/"+output+"/"+lxyvar+"/fprompt_fithist.root")<<std::endl;
      TFile* inf = TFile::Open(std::string("rootfiles/"+output+"/"+lxyvar+"/fprompt_fithist.root").c_str());
      // <--
      fitX::init(inf);
      grfprompt_a[lxyvar].resize(lxydis::lxycut[lxyvar].size());
      grfprompt_b[lxyvar].resize(lxydis::lxycut[lxyvar].size());
      for(int i=0; i<lxydis::lxycut[lxyvar].size(); i++)
        {
          grfprompt_a[lxyvar][i] = (TEfficiency*)inf->Get(Form("grfprompt_a-%d", i));
          grfprompt_a[lxyvar][i]->SetName(Form("grfprompt_a-%s-%d", lxyvar.c_str(), i));
          grfprompt_b[lxyvar][i] = (TEfficiency*)inf->Get(Form("grfprompt_b-%d", i));
          grfprompt_b[lxyvar][i]->SetName(Form("grfprompt_b-%s-%d", lxyvar.c_str(), i));
          x[n]     = n+0.5;
          ex[n]    = 0.5;
          y_a[n]   = grfprompt_a[lxyvar][i]->GetEfficiency(fitX::ibin_a);
          eyl_a[n] = grfprompt_a[lxyvar][i]->GetEfficiencyErrorLow(fitX::ibin_a);
          eyh_a[n] = grfprompt_a[lxyvar][i]->GetEfficiencyErrorUp(fitX::ibin_a);
          y_b[n]   = grfprompt_b[lxyvar][i]->GetEfficiency(fitX::ibin_b);
          eyl_b[n] = grfprompt_b[lxyvar][i]->GetEfficiencyErrorLow(fitX::ibin_b);
          eyh_b[n] = grfprompt_b[lxyvar][i]->GetEfficiencyErrorUp(fitX::ibin_b);
          
          hemptyfprompt->GetXaxis()->SetBinLabel(n+1, xjjc::str_replaceall(lxydis::vars[lxyvar], "mm", xjjc::number_remove_zero(lxydis::lxycut[lxyvar][i])).c_str());
          n++;
        }
    }
  x[nn-1]     = nn-1+0.5;
  ex[nn-1]    = 0.5;
  y_a[nn-1]   = fprompt_a;
  eyl_a[nn-1] = errfprompt_a;
  eyh_a[nn-1] = errfprompt_a;
  hemptyfprompt->GetXaxis()->SetBinLabel(nn, "l_{xy} Fit");
      
  TGraphAsymmErrors* gfprompt_a = new TGraphAsymmErrors(nn, x.data(), y_a.data(), ex.data(), ex.data(), eyl_a.data(), eyh_a.data());
  gfprompt_a->SetName("gfprompt_a");
  xjjroot::setthgrstyle(gfprompt_a, fitX::color_a, 25, 1.4, fitX::color_a, 1, 3);
  TGraphAsymmErrors* gfprompt_b = new TGraphAsymmErrors(nn, x.data(), y_b.data(), ex.data(), ex.data(), eyl_b.data(), eyh_b.data());
  gfprompt_b->SetName("gfprompt_b");
  xjjroot::setthgrstyle(gfprompt_b, fitX::color_b, 25, 1.4, fitX::color_b, 1, 3);

  TLegend* leg = new TLegend(0.22, 0.86-0.049*2, 0.58, 0.86);
  xjjroot::setleg(leg, 0.042);
  leg->AddEntry(gfprompt_a, fitX::title_a.c_str(), "pl");
  leg->AddEntry(gfprompt_b, fitX::title_b.c_str(), "pl");

  xjjroot::setgstyle(2);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  hemptyfprompt->Draw();
  xjjroot::drawline(xmin, 1, xmax, 1, kGray+1, 9, 2);
  gfprompt_a->Draw("pe same");
  gfprompt_b->Draw("pe same");
  xjjroot::drawline(xmin, y_a[0], xmax, y_a[0], fitX::color_a, 2, 5, 0.6);
  xjjroot::drawline(xmin, y_b[0], xmax, y_b[0], fitX::color_b, 2, 5, 0.6);
  leg->Draw();
  fitX::drawkinematics();
  xjjroot::drawcomment(output.c_str(), "r");
  xjjroot::drawCMS();
  xjjroot::mkdir(Form("plots/%s/cfpromptsyst.pdf", output.c_str()));
  c->SaveAs(Form("plots/%s/cfpromptsyst.pdf", output.c_str()));

  float syst_a = 0, syst_b = 0, syst_r = 0;
  for(int i=0; i<nn; i++)
    {
      if(fabs(y_a[i]-y_a[0]) > syst_a) { syst_a = fabs(y_a[i]-y_a[0]); }
      if(i==nn-1) continue;
      if(fabs(y_b[i]-y_b[0]) > syst_b) { syst_b = fabs(y_b[i]-y_b[0]); }
    }
  syst_a /= y_a[0];
  syst_b /= y_b[0];
  syst_r = TMath::Sqrt(syst_a*syst_a + syst_b*syst_b);
  
  std::cout<<"Prompt Fraction & "<<Form("%.1f", syst_a*1.e+2)<<"\\% & "<<Form("%.1f", syst_b*1.e+2)<<"\\% & "<<Form("%.1f", syst_r*1.e+2)<<"\\% \\\\"<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==2) { fprompt_drawhist(argv[1]); return 0; }
  return 1;
}

