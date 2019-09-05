#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TH1F.h>
#include "packtree.h"
#include "ntuple.h"
#include "xjjcuti.h"

#include "HEPData-ins1495026-v1/ppATLAS.h"

void fillptshape(std::string inputname, std::string type, std::string outputname, std::string inputdir="HEPData-ins1495026-v1")
{
  std::string outdir = xjjc::str_replaceall(outputname, xjjc::str_divide(outputname, "/").back(), "");
  gSystem->Exec(Form("mkdir -p %s", outdir.c_str()));

  ppref::ppATLAS getpp(inputdir);

  std::cout<<"==> Opening files"<<std::endl;
  std::cout<<"input: "<<inputname<<std::endl;
  std::cout<<"output: "<<outputname<<std::endl;

  std::cout<<"==> Building histograms"<<std::endl;
  TH1F* hmc = new TH1F("hmc", ";Gen p_{T} (GeV/c);", 24, 10, 70); //
  hmc->Sumw2();
  TH1F* hmcweight = new TH1F("hmcweight", ";Gen p_{T} (GeV/c);", 24, 10, 70); //
  hmcweight->Sumw2();

  TFile* inf = TFile::Open(inputname.c_str());
  xjjroot::packtree* pt = new xjjroot::packtree(inf, "Bfinder/ntmix", "mcp_ptshape", "Bfinder/ntGen");
  mytmva::ntuple* ntp = pt->ntp();

  std::cout<<"==> Scaning file"<<std::endl;
  int nentries = ntp->gnt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      pt->getentry(i);
      if(i%1000 == 0) { xjjc::progressbar(i, nentries); }

      float weight = ntp->pthatweight * ntp->Ncoll;
      for(int j=0; j<ntp->Gsize; j++)
        {
          // if(xjjc::str_contains(inputname, "PromptXRho") && ntp->Gpt[j] < 15.) continue; // !!
          if(ntp->Gpt[j] < 15.) continue; // !!
          if(!(ntp->GisSignal[j]==7 && ntp->GcollisionId[j]==0)) continue;
          if(!(TMath::Abs(ntp->Gy[j]) < 0.75)) continue; //
          hmc->Fill(ntp->Gpt[j], weight);
          hmcweight->Fill(ntp->Gpt[j], weight*getpp.getweight(ntp->Gpt[j], type));
        }
    }
  xjjc::progressbar_summary(nentries);

  TFile* outf = new TFile(outputname.c_str(), "recreate");
  hmc->Write();
  hmcweight->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==5) { fillptshape(argv[1], argv[2], argv[3], argv[4]); return 0; }
  if(argc==4) { fillptshape(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}
