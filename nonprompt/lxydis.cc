#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>

#include <string>
#include <iostream>
#include <map>

#include "ntuple.h"
#include "packtree.h"
#include "xjjcuti.h"

#include "lxydis.h"

void lxydisfun(std::string inputname, std::string type, float mvaval, std::string outputname)
{
  std::map<std::string, std::vector<float>> xbins = lxydis::setupbins();
  
  //
  std::map<std::string, TH1F*> hBlxy, hBlxySignalRegionL, hBlxySignalRegionH, hBlxySidebandL, hBlxySidebandH;
  for(auto& vv : lxydis::vars) 
    {
      std::string var = vv.first, vartitle = vv.second;
      hBlxy[var]              = new TH1F(Form("hB%s",              var.c_str()), Form(";%s;Probability", vartitle.c_str()), xbins[var].size()-1, xbins[var].data()); hBlxy[var]->Sumw2();
      hBlxySignalRegionL[var] = new TH1F(Form("hB%sSignalRegionL", var.c_str()), Form(";%s;Probability", vartitle.c_str()), xbins[var].size()-1, xbins[var].data()); hBlxySignalRegionL[var]->Sumw2();
      hBlxySignalRegionH[var] = new TH1F(Form("hB%sSignalRegionH", var.c_str()), Form(";%s;Probability", vartitle.c_str()), xbins[var].size()-1, xbins[var].data()); hBlxySignalRegionH[var]->Sumw2();
      hBlxySidebandL[var]     = new TH1F(Form("hB%sSidebandL",     var.c_str()), Form(";%s;Probability", vartitle.c_str()), xbins[var].size()-1, xbins[var].data()); hBlxySidebandL[var]->Sumw2();
      hBlxySidebandH[var]     = new TH1F(Form("hB%sSidebandH",     var.c_str()), Form(";%s;Probability", vartitle.c_str()), xbins[var].size()-1, xbins[var].data()); hBlxySidebandH[var]->Sumw2();
    }
  bool ismc = !(xjjc::str_contains(type, "data") || xjjc::str_contains(type, "samesign"));

  std::vector<TH1F*> hBlxyBdtgPromptL(lxydis::mvalist.size()), hBlxyBdtgPromptH(lxydis::mvalist.size()), 
    hBlxyBdtgNonpromptL(lxydis::mvalist.size()), hBlxyBdtgNonpromptH(lxydis::mvalist.size());
  for(int k=0; k<lxydis::mvalist.size(); k++)
    {
      hBlxyBdtgPromptL[k] = new TH1F(Form("hBlxyBdtgPromptL_%d", k), ";l_{xy} (mm);Probability", xbins["lxyprompt"].size()-1, xbins["lxyprompt"].data()); 
      hBlxyBdtgPromptL[k]->Sumw2();
      hBlxyBdtgPromptH[k] = new TH1F(Form("hBlxyBdtgPromptH_%d", k), ";l_{xy} (mm);Probability", xbins["lxyprompt"].size()-1, xbins["lxyprompt"].data()); 
      hBlxyBdtgPromptH[k]->Sumw2();
      hBlxyBdtgNonpromptL[k] = new TH1F(Form("hBlxyBdtgNonpromptL_%d", k), ";l_{xy} (mm);Probability", xbins["lxynonprompt"].size()-1, xbins["lxynonprompt"].data()); 
      hBlxyBdtgNonpromptL[k]->Sumw2();
      hBlxyBdtgNonpromptH[k] = new TH1F(Form("hBlxyBdtgNonpromptH_%d", k), ";l_{xy} (mm);Probability", xbins["lxynonprompt"].size()-1, xbins["lxynonprompt"].data()); 
      hBlxyBdtgNonpromptH[k]->Sumw2();
    }

  //
  TFile* inf = TFile::Open(inputname.c_str()); 
  xjjroot::packtree* ptr = new xjjroot::packtree(inf, "Bfinder/ntmix", type);
  mytmva::ntuple* ntp = ptr->ntp;
  // float mva; ntp->getnt()->SetBranchAddress("", &mva); // !!
  if(ismc && !ntp->isweight())
    { std::cout<<__FUNCTION__<<": error: weight is not correctly placed in MC sample."<<std::endl; return; }

  int nentries = ptr->getentries();
  for(int i=0; i<nentries; i++)
    {
      ptr->getentry(i);
      if(i%1000 == 0) xjjc::progressbar(i, nentries);

      if(!ntp->passedevtfil()) continue;

      float weight = ntp->isweight()?(ntp->pthatweight*ntp->Ncoll):1;

      for(int j=0; j<ntp->Bsize; j++)
        {
          if(ismc && !(ntp->Bgen[j]==23333 && ntp->BgencollisionId[j]==0)) continue;
          if(!ntp->passedpre(j)) continue;
          
          if(ntp->Bpt[j] < 15 || TMath::Abs(ntp->By[j]) > 1.5) continue;

          //
          std::map<std::string, float> fillval;
          fillval[ "lxy" ] = ntp->Blxy[j];
          fillval[ "dls" ] = ntp->BsvpvDistance[j] / ntp->BsvpvDisErr[j];
          fillval[ "dca" ] = ntp->BsvpvDistance[j] * TMath::Sin(ntp->Balpha[j]);
          
          for(int k=0; k<lxydis::mvalist.size(); k++)
            {
              if(ntp->BDTG[j] > lxydis::mvalist[k])
                {
                  hBlxyBdtgPromptL[k]->Fill(fillval["lxy"], weight);
                  hBlxyBdtgPromptH[k]->Fill(fillval["lxy"], weight);
                  hBlxyBdtgNonpromptL[k]->Fill(fillval["lxy"], weight);
                  hBlxyBdtgNonpromptH[k]->Fill(fillval["lxy"], weight);
                }
            }

          if(ntp->BDTG[j] < mvaval) continue; // !!

          // 
          for(auto& vv : lxydis::vars) 
            { 
              hBlxy[vv.first]->Fill(fillval[vv.first], weight);
              // SignalRegionL
              if(ntp->signalregionl(j)) // 
                { hBlxySignalRegionL[vv.first]->Fill(fillval[vv.first], weight); }
              // SignalRegionH
              if(ntp->signalregionh(j)) // 
                { hBlxySignalRegionH[vv.first]->Fill(fillval[vv.first], weight); }
              // SidebandL
              if(ntp->sidebandl(j)) // 
                { hBlxySidebandL[vv.first]->Fill(fillval[vv.first], weight); }
              // SidebandH
              if(ntp->sidebandh(j)) // 
                { hBlxySidebandH[vv.first]->Fill(fillval[vv.first], weight); }
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  inf->Close();

  TFile* outf = new TFile(Form("%s.root", outputname.c_str()), "recreate");
  outf->cd();
  for(auto& vv : lxydis::vars)
    {
      hBlxy[vv.first]->Write();
      hBlxySignalRegionL[vv.first]->Write();
      hBlxySignalRegionH[vv.first]->Write();
      hBlxySidebandL[vv.first]->Write();
      hBlxySidebandH[vv.first]->Write();
    }
  for(int k=0; k<lxydis::mvalist.size(); k++)
    {
      hBlxyBdtgPromptL[k]->Write();
      hBlxyBdtgPromptH[k]->Write();
      hBlxyBdtgNonpromptL[k]->Write();
      hBlxyBdtgNonpromptH[k]->Write();
    }
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==5) { lxydisfun(argv[1], argv[2], atof(argv[3]), argv[4]); return 0; }
  return 1;
}
