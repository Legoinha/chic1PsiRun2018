#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <string>

#include <iostream>

#include "xjjrootuti.h"
#include "lxydis.h"

void drawlxydis(std::string inputdata, std::string inputsamesign,
                std::string inputpromptpsi, std::string inputnonpromptpsi,
                std::string inputpromptx, std::string inputnonpromptx)
{
  xjjroot::setgstyle(1);
  lxydis::setupdraw();
  std::map<std::string, std::vector<float>> xbins = lxydis::setupbins();

  // >>>
  TFile* infyield = TFile::Open("rootfilesCrossSection/root_fitXmc_nodls.root");
  TH1F* hyield = (TH1F*)infyield->Get("hyield");
  TH1F* hbkg = (TH1F*)infyield->Get("hbkg");
  TH1F* htotal = (TH1F*)hyield->Clone("htotal");
  htotal->Add(hbkg);
  TH1F* hbkgfrac = (TH1F*)hbkg->Clone("hbkgfrac");
  hbkgfrac->Divide(htotal);
  float ybkgfrac_a = hbkgfrac->GetBinContent(2);
  float ybkgfrac_b = hbkgfrac->GetBinContent(4);
  std::cout<<ybkgfrac_a<<" "<<ybkgfrac_b<<std::endl;
  // <<<

  TFile* infdata = TFile::Open(inputdata.c_str());
  TFile* infsamesign = TFile::Open(inputsamesign.c_str());
  TFile* infpromptpsi = TFile::Open(inputpromptpsi.c_str());
  TFile* infnonpromptpsi = TFile::Open(inputnonpromptpsi.c_str());
  TFile* infpromptx = TFile::Open(inputpromptx.c_str());
  TFile* infnonpromptx = TFile::Open(inputnonpromptx.c_str());

  std::map<std::string, TFile*> infs = 
    {
      std::pair<std::string, TFile*>("data", infdata),
      std::pair<std::string, TFile*>("samesign", infsamesign),
      std::pair<std::string, TFile*>("promptpsi", infpromptpsi),
      std::pair<std::string, TFile*>("nonpromptpsi", infnonpromptpsi),
      std::pair<std::string, TFile*>("promptx", infpromptx),
      std::pair<std::string, TFile*>("nonpromptx", infnonpromptx),
    };

  std::vector<TH1F*> hBlxyBdtgPromptL(lxydis::mvalist.size()), hBlxyBdtgPromptH(lxydis::mvalist.size()),
    hBlxyBdtgNonpromptL(lxydis::mvalist.size()), hBlxyBdtgNonpromptH(lxydis::mvalist.size());
  for(int k=0; k<lxydis::mvalist.size(); k++)
    {
      hBlxyBdtgPromptL[k] = (TH1F*)infs["promptpsi"]->Get(Form("hBlxyBdtgPromptL_%d", k));
      hBlxyBdtgPromptL[k]->Scale(1./hBlxyBdtgPromptL[k]->Integral(), "width");
      hBlxyBdtgPromptH[k] = (TH1F*)infs["promptx"]->Get(Form("hBlxyBdtgPromptH_%d", k));
      hBlxyBdtgPromptH[k]->Scale(1./hBlxyBdtgPromptH[k]->Integral(), "width");
      hBlxyBdtgNonpromptL[k] = (TH1F*)infs["nonpromptpsi"]->Get(Form("hBlxyBdtgNonpromptL_%d", k));
      hBlxyBdtgNonpromptL[k]->Scale(1./hBlxyBdtgNonpromptL[k]->Integral(), "width");
      hBlxyBdtgNonpromptH[k] = (TH1F*)infs["nonpromptx"]->Get(Form("hBlxyBdtgNonpromptH_%d", k));
      hBlxyBdtgNonpromptH[k]->Scale(1./hBlxyBdtgNonpromptH[k]->Integral(), "width");

      xjjroot::setthgrstyle(hBlxyBdtgPromptL[k],
                            lxydis::drawset["SignalRegion_prompt"]->mcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_prompt"]->mstyle, lxydis::drawset["SignalRegion_prompt"]->msize,
                            lxydis::drawset["SignalRegion_prompt"]->lcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_prompt"]->lstyle, lxydis::drawset["SignalRegion_prompt"]->lwidth);
      xjjroot::setthgrstyle(hBlxyBdtgPromptH[k],
                            lxydis::drawset["SignalRegion_prompt"]->mcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_prompt"]->mstyle, lxydis::drawset["SignalRegion_prompt"]->msize,
                            lxydis::drawset["SignalRegion_prompt"]->lcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_prompt"]->lstyle, lxydis::drawset["SignalRegion_prompt"]->lwidth);
      xjjroot::setthgrstyle(hBlxyBdtgNonpromptL[k],
                            lxydis::drawset["SignalRegion_nonprompt"]->mcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_nonprompt"]->mstyle, lxydis::drawset["SignalRegion_nonprompt"]->msize,
                            lxydis::drawset["SignalRegion_nonprompt"]->lcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_nonprompt"]->lstyle, lxydis::drawset["SignalRegion_nonprompt"]->lwidth);
      xjjroot::setthgrstyle(hBlxyBdtgNonpromptH[k],
                            lxydis::drawset["SignalRegion_nonprompt"]->mcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_nonprompt"]->mstyle, lxydis::drawset["SignalRegion_nonprompt"]->msize,
                            lxydis::drawset["SignalRegion_nonprompt"]->lcolor+k-lxydis::thismva, lxydis::drawset["SignalRegion_nonprompt"]->lstyle, lxydis::drawset["SignalRegion_nonprompt"]->lwidth);
    }

  std::map<std::string, std::map<std::string,TH1F*>> hBlxy, hBlxySignalRegionL, hBlxySignalRegionH, hBlxySidebandL, hBlxySidebandH;
  for(auto& vv : lxydis::vars)
    {
      std::string var = vv.first, vartitle = vv.second;
      for(auto& input : infs)
        {
          std::string type = input.first;
          hBlxy[type][var] = (TH1F*)infs[type]->Get(Form("hB%s", var.c_str())); 
          hBlxy[type][var]->SetName(Form("%s_%s", hBlxy[type][var]->GetName(), type.c_str()));
          hBlxy[type][var]->Scale(1./hBlxy[type][var]->Integral(), "width");
          hBlxySignalRegionL[type][var] = (TH1F*)infs[type]->Get(Form("hB%sSignalRegionL", var.c_str())); 
          hBlxySignalRegionL[type][var]->SetName(Form("%s_%s", hBlxySignalRegionL[type][var]->GetName(), type.c_str()));
          hBlxySignalRegionL[type][var]->Scale(1./hBlxySignalRegionL[type][var]->Integral(), "width");
          hBlxySignalRegionH[type][var] = (TH1F*)infs[type]->Get(Form("hB%sSignalRegionH", var.c_str())); 
          hBlxySignalRegionH[type][var]->SetName(Form("%s_%s", hBlxySignalRegionH[type][var]->GetName(), type.c_str()));
          hBlxySignalRegionH[type][var]->Scale(1./hBlxySignalRegionH[type][var]->Integral(), "width");
          hBlxySidebandL[type][var] = (TH1F*)infs[type]->Get(Form("hB%sSidebandL", var.c_str())); 
          hBlxySidebandL[type][var]->SetName(Form("%s_%s", hBlxySidebandL[type][var]->GetName(), type.c_str()));
          hBlxySidebandL[type][var]->Scale(1./hBlxySidebandL[type][var]->Integral(), "width");
          hBlxySidebandH[type][var] = (TH1F*)infs[type]->Get(Form("hB%sSidebandH", var.c_str())); 
          hBlxySidebandH[type][var]->SetName(Form("%s_%s", hBlxySidebandH[type][var]->GetName(), type.c_str()));
          hBlxySidebandH[type][var]->Scale(1./hBlxySidebandH[type][var]->Integral(), "width");
        }
      hBlxySignalRegionL["samesign"][var]->Scale(ybkgfrac_a);
      hBlxySidebandL["data"][var]->Scale(ybkgfrac_a);
      hBlxySignalRegionH["samesign"][var]->Scale(ybkgfrac_b);
      hBlxySidebandH["data"][var]->Scale(ybkgfrac_b);
      xjjroot::setthgrstyle(hBlxySignalRegionL["data"][var], 
                            lxydis::drawset["SignalRegion_data"]->mcolor, lxydis::drawset["SignalRegion_data"]->mstyle, lxydis::drawset["SignalRegion_data"]->msize, 
                            lxydis::drawset["SignalRegion_data"]->lcolor, lxydis::drawset["SignalRegion_data"]->lstyle, lxydis::drawset["SignalRegion_data"]->lwidth);
      xjjroot::setthgrstyle(hBlxySidebandL["data"][var], 
                            lxydis::drawset["Sideband_data"]->mcolor, lxydis::drawset["Sideband_data"]->mstyle, lxydis::drawset["Sideband_data"]->msize, 
                            lxydis::drawset["Sideband_data"]->lcolor, lxydis::drawset["Sideband_data"]->lstyle, lxydis::drawset["Sideband_data"]->lwidth);
      xjjroot::setthgrstyle(hBlxySignalRegionL["samesign"][var], 
                            lxydis::drawset["SignalRegion_samesign"]->mcolor, lxydis::drawset["SignalRegion_samesign"]->mstyle, lxydis::drawset["SignalRegion_samesign"]->msize, 
                            lxydis::drawset["SignalRegion_samesign"]->lcolor, lxydis::drawset["SignalRegion_samesign"]->lstyle, lxydis::drawset["SignalRegion_samesign"]->lwidth);
      xjjroot::setthgrstyle(hBlxySignalRegionL["promptpsi"][var], 
                            lxydis::drawset["SignalRegion_prompt"]->mcolor, lxydis::drawset["SignalRegion_prompt"]->mstyle, lxydis::drawset["SignalRegion_prompt"]->msize, 
                            lxydis::drawset["SignalRegion_prompt"]->lcolor, lxydis::drawset["SignalRegion_prompt"]->lstyle, lxydis::drawset["SignalRegion_prompt"]->lwidth);
      xjjroot::setthgrstyle(hBlxySignalRegionL["nonpromptpsi"][var], 
                            lxydis::drawset["SignalRegion_nonprompt"]->mcolor, lxydis::drawset["SignalRegion_nonprompt"]->mstyle, lxydis::drawset["SignalRegion_nonprompt"]->msize, 
                            lxydis::drawset["SignalRegion_nonprompt"]->lcolor, lxydis::drawset["SignalRegion_nonprompt"]->lstyle, lxydis::drawset["SignalRegion_nonprompt"]->lwidth);

      xjjroot::setthgrstyle(hBlxySignalRegionH["data"][var], 
                            lxydis::drawset["SignalRegion_data"]->mcolor, lxydis::drawset["SignalRegion_data"]->mstyle, lxydis::drawset["SignalRegion_data"]->msize, 
                            lxydis::drawset["SignalRegion_data"]->lcolor, lxydis::drawset["SignalRegion_data"]->lstyle, lxydis::drawset["SignalRegion_data"]->lwidth);
      xjjroot::setthgrstyle(hBlxySidebandH["data"][var], 
                            lxydis::drawset["Sideband_data"]->mcolor, lxydis::drawset["Sideband_data"]->mstyle, lxydis::drawset["Sideband_data"]->msize, 
                            lxydis::drawset["Sideband_data"]->lcolor, lxydis::drawset["Sideband_data"]->lstyle, lxydis::drawset["Sideband_data"]->lwidth);
      xjjroot::setthgrstyle(hBlxySignalRegionH["samesign"][var], 
                            lxydis::drawset["SignalRegion_samesign"]->mcolor, lxydis::drawset["SignalRegion_samesign"]->mstyle, lxydis::drawset["SignalRegion_samesign"]->msize, 
                            lxydis::drawset["SignalRegion_samesign"]->lcolor, lxydis::drawset["SignalRegion_samesign"]->lstyle, lxydis::drawset["SignalRegion_samesign"]->lwidth);
      xjjroot::setthgrstyle(hBlxySignalRegionH["promptx"][var], 
                            lxydis::drawset["SignalRegion_prompt"]->mcolor, lxydis::drawset["SignalRegion_prompt"]->mstyle, lxydis::drawset["SignalRegion_prompt"]->msize, 
                            lxydis::drawset["SignalRegion_prompt"]->lcolor, lxydis::drawset["SignalRegion_prompt"]->lstyle, lxydis::drawset["SignalRegion_prompt"]->lwidth);
      xjjroot::setthgrstyle(hBlxySignalRegionH["nonpromptx"][var], 
                            lxydis::drawset["SignalRegion_nonprompt"]->mcolor, lxydis::drawset["SignalRegion_nonprompt"]->mstyle, lxydis::drawset["SignalRegion_nonprompt"]->msize, 
                            lxydis::drawset["SignalRegion_nonprompt"]->lcolor, lxydis::drawset["SignalRegion_nonprompt"]->lstyle, lxydis::drawset["SignalRegion_nonprompt"]->lwidth);

    }

  TCanvas* ccm = new TCanvas("ccm", "", 1200, 1200);
  ccm->Divide(2, 2);
  for(int pp=1; pp<=4; pp++)
    {
      TPad* p = (TPad*)(ccm->cd(pp));
      p->SetLogy();      
    }
  int nbinprompt = hBlxyBdtgPromptL[lxydis::thismva]->GetXaxis()->FindBin(0.1)-1;
  int nbinnonprompt = hBlxyBdtgNonpromptL[lxydis::thismva]->GetXaxis()->FindBin(0.1)-1;
  // std::cout<<"hBlxyBdtgPromptL "<<hBlxyBdtgPromptL[lxydis::thismva]->GetXaxis()->FindBin(0.1)<<std::endl; 
  // std::cout<<"hBlxyBdtgNonpromptL "<<hBlxyBdtgNonpromptL[lxydis::thismva]->GetXaxis()->FindBin(0.1)<<std::endl;
  double fltLerr, fgtLerr, fltHerr, fgtHerr;
  double fltL = hBlxyBdtgNonpromptL[lxydis::thismva]->IntegralAndError(1,               nbinnonprompt,       fltLerr, "width");
  double fgtL = hBlxyBdtgNonpromptL[lxydis::thismva]->IntegralAndError(nbinnonprompt+1, xbins["lxynonprompt"].size(), fgtLerr, "width");
  double fltH = hBlxyBdtgNonpromptH[lxydis::thismva]->IntegralAndError(1,               nbinnonprompt,       fltHerr, "width");
  double fgtH = hBlxyBdtgNonpromptH[lxydis::thismva]->IntegralAndError(nbinnonprompt+1, xbins["lxynonprompt"].size(), fgtHerr, "width");
  TH1F* hfgt = new TH1F("hfgt", "", 5, 0, 5);
  hfgt->SetBinContent(2, fgtL);
  hfgt->SetBinError(2, fgtLerr);
  hfgt->SetBinContent(4, fgtH);
  hfgt->SetBinError(4, fgtHerr);

  for(int k=0; k<lxydis::mvalist.size(); k++)
    {
      if(!k) 
        {
          xjjroot::sethempty(hBlxyBdtgPromptL[k], 0, 0);
          xjjroot::sethempty(hBlxyBdtgPromptH[k], 0, 0);
          xjjroot::sethempty(hBlxyBdtgNonpromptL[k], 0, 0);
          xjjroot::sethempty(hBlxyBdtgNonpromptH[k], 0, 0);
        }
      std::string opt = k?"histe same":"histe";
      ccm->cd(1);
      hBlxyBdtgPromptL[k]->Draw(opt.c_str());
      xjjroot::drawline(0.1, hBlxyBdtgPromptL[k]->GetMinimum(), 0.1, hBlxyBdtgPromptL[k]->GetMaximum(), kGray, 2, 2);
      ccm->cd(2);
      hBlxyBdtgPromptH[k]->Draw(opt.c_str());
      xjjroot::drawline(0.1, hBlxyBdtgPromptH[k]->GetMinimum(), 0.1, hBlxyBdtgPromptH[k]->GetMaximum(), kGray, 2, 2);
      ccm->cd(3);
      hBlxyBdtgNonpromptL[k]->Draw(opt.c_str());
      xjjroot::drawline(0.1, hBlxyBdtgNonpromptL[k]->GetMinimum(), 0.1, hBlxyBdtgNonpromptL[k]->GetMaximum(), kGray, 2, 2);
      ccm->cd(4);
      hBlxyBdtgNonpromptH[k]->Draw(opt.c_str());
      xjjroot::drawline(0.1, hBlxyBdtgNonpromptH[k]->GetMinimum(), 0.1, hBlxyBdtgNonpromptH[k]->GetMaximum(), kGray, 2, 2);
    }
  for(int pp=1; pp<=4; pp++)
    {
      ccm->cd(pp);
      xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
      xjjroot::drawCMSright();
      std::string channel = pp%2==1?"#psi(2S)":"X(3872)";
      std::string pnp = (pp-1)/2==0?"Prompt":"Nonprompt";
      xjjroot::drawtex(0.91, 0.85, channel.c_str(), 0.04, 32, 62);
      xjjroot::drawtex(0.91, 0.80, pnp.c_str(), 0.04, 32, 62);
      if(pp == 3)
        {
          xjjroot::drawtex(0.20, 0.30, Form("%.2f%s", (1-fgtL)*100., "%"), 0.04, 12, 62);
          xjjroot::drawtex(0.50, 0.30, Form("%.2f%s", fgtL*100., "%"), 0.04, 12, 62);
        }
      if(pp == 4)
        {
          xjjroot::drawtex(0.20, 0.30, Form("%.2f%s", (1-fgtH)*100., "%"), 0.04, 12, 62);
          xjjroot::drawtex(0.50, 0.30, Form("%.2f%s", fgtH*100, "%"), 0.04, 12, 62);
        }
    }
  ccm->SaveAs("plots/clxydis_bdtg.pdf");

  TLegend* leg = new TLegend(0.54, 0.88-5*0.04, 0.94, 0.88);
  xjjroot::setleg(leg, 0.035);
  leg->AddEntry(hBlxySignalRegionL["data"]["lxy"], "data signal-region", "pl");
  leg->AddEntry(hBlxySidebandL["data"]["lxy"], "data side-band", "l");
  leg->AddEntry(hBlxySignalRegionL["samesign"]["lxy"], "data same-sign", "l");
  leg->AddEntry(hBlxySignalRegionL["promptpsi"]["lxy"], "MC prompt", "l");
  leg->AddEntry(hBlxySignalRegionL["nonpromptpsi"]["lxy"], "MC nonprompt", "l");

  for(auto& vv : lxydis::vars)
    {
      std::string var = vv.first;
      
      TCanvas* c = new TCanvas("c", "", 1200, 600);
      c->Divide(2, 1);

      TPad* p1 = (TPad*)(c->cd(1));
      p1->SetLogy();
      xjjroot::sethempty(hBlxySignalRegionL["data"][var], 0, 0);
      hBlxySignalRegionL["data"][var]->Draw("pe"); 
      hBlxySidebandL["data"][var]->Draw("histe same");  
      hBlxySignalRegionL["samesign"][var]->Draw("histe same");
      hBlxySignalRegionL["promptpsi"][var]->Draw("histe same");
      hBlxySignalRegionL["nonpromptpsi"][var]->Draw("histe same");
      hBlxySignalRegionL["data"][var]->Draw("pe same");  
      leg->Draw();
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();

      TPad* p2 = (TPad*)(c->cd(2));
      p2->SetLogy();
      xjjroot::sethempty(hBlxySignalRegionH["data"][var], 0, 0);
      hBlxySignalRegionH["data"][var]->Draw("pe");  
      hBlxySidebandH["data"][var]->Draw("histe same");  
      hBlxySignalRegionH["samesign"][var]->Draw("histe same");
      hBlxySignalRegionH["promptx"][var]->Draw("histe same");
      hBlxySignalRegionH["nonpromptx"][var]->Draw("histe same");
      hBlxySignalRegionH["data"][var]->Draw("pe same");  
      leg->Draw();
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();

      c->SaveAs(Form("plots/clxydis_%s.pdf", var.c_str()));

      delete p2;
      delete p1;
      delete c;
    }

  TFile* outf = new TFile("rootfiles/root_lxydis.root", "recreate");
  outf->cd();
  hfgt->Write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc == 7)
    {
      drawlxydis(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);      
    }
}
