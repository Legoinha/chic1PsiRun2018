#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

int evtmatching(std::string inputname, std::string inputnamemuon, std::string outfile)
{
  TFile* inf = TFile::Open(inputname.c_str());
  TFile* infmuon = TFile::Open(inputnamemuon.c_str());

  TTree* hi = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  TTree* muoninfo = (TTree*)infmuon->Get("Bfinder/root");
  TTree* hi_mu = (TTree*)infmuon->Get("hiEvtAnalyzer/HiTree");

  UInt_t run; hi->SetBranchAddress("run", &run);
  ULong64_t evt; hi->SetBranchAddress("evt", &evt);
  UInt_t lumi; hi->SetBranchAddress("lumi", &lumi);
  float sample; hi->SetBranchAddress("sample", &sample);

  UInt_t run_mu; hi_mu->SetBranchAddress("run", &run_mu);
  ULong64_t evt_mu; hi_mu->SetBranchAddress("evt", &evt_mu);
  UInt_t lumi_mu;hi_mu->SetBranchAddress("lumi", &lumi_mu);
  float sample_mu; hi_mu->SetBranchAddress("sample", &sample_mu);

  TFile* outf = new TFile(outfile.c_str(), "recreate");
  TDirectory* dmuoninfo = outf->mkdir("muoninfo", "");
  dmuoninfo->cd();
  TTree* new_muoninfo = muoninfo->CloneTree(0);

  std::string outname_unmatch = xjjc::str_replaceall(xjjc::str_divide_lastel(outfile, "/"), ".root", "");
  outname_unmatch = "unmatch/un_" + xjjc::str_replaceallspecial(outname_unmatch) + ".txt";
  xjjroot::mkdir(outname_unmatch);
  std::ofstream fout_unmatch(outname_unmatch.c_str());

  std::cout<<"---- Event reading"<<std::endl;
  Long64_t nentries = hi->GetEntries();
  Long64_t nentries_mu = hi_mu->GetEntries();
  std::vector<UInt_t> hi_mu_run;
  std::vector<ULong64_t> hi_mu_evt;
  std::vector<UInt_t> hi_mu_lumi;
  std::vector<float> hi_mu_sample;
  float current_sample = -1;
  std::vector<Long64_t> ksample;
  for(Long64_t k = 0; k<nentries_mu; k++)
    {
      if(k % 10000 == 0 || k == nentries_mu-1) xjjc::progressbar(k, nentries_mu);
      hi_mu->GetEntry(k); //
      hi_mu_run.push_back(run_mu);
      hi_mu_evt.push_back(evt_mu);
      hi_mu_lumi.push_back(lumi_mu);      
      hi_mu_sample.push_back(sample_mu);
      
      if(sample_mu != current_sample) { ksample.push_back(k); current_sample = sample_mu; }
    }
  xjjc::progressbar_summary(nentries_mu);
  ksample.push_back(nentries_mu);

  for(auto& ks : ksample) std::cout<<ks<<" ";
  std::cout<<std::endl;

  std::cout<<std::endl<<"---- Extracted event info"<<std::endl;
  int countmatch = 0;
  std::vector<Long64_t> matchingtable;
  current_sample = -1;
  int i_ksample = -1;
  for(Long64_t j = 0; j<nentries; j++)
    {
      if(j % 1000 == 0 || j == nentries-1) xjjc::progressbar(j, nentries);

      hi->GetEntry(j); //

      if(sample != current_sample) { i_ksample++; current_sample = sample; }
      Long64_t k = ksample[i_ksample];
      for(; k < ksample[i_ksample+1]; k++)
        {
          if(hi_mu_evt[k]==evt && hi_mu_lumi[k]==lumi && hi_mu_run[k]==run && hi_mu_sample[k]==sample)
            break;
        }
      if(k>=ksample[i_ksample+1]) k = -1;
      matchingtable.push_back(k);
      if(k>=0) { countmatch++; }
      else { fout_unmatch << j << std::endl; }
    }
  xjjc::progressbar_summary(nentries);
  std::cout<<"countmatch = "<<countmatch<<" / "<<nentries<<std::endl;
  if(countmatch < nentries) { std::cout<<"error: Not all event matched."<<std::endl; /* return 1; */ }

  std::cout<<"---- Writing to tree"<<std::endl;
  for(Long64_t j = 0; j<nentries; j++)
    {
      if(j % 10000 == 0 || j == nentries-1) xjjc::progressbar(j, nentries);
      Long64_t k = matchingtable[j];
      if(k<0) continue;

      muoninfo->GetEntry(k); //

      dmuoninfo->cd();
      new_muoninfo->Fill();
    }
  xjjc::progressbar_summary(nentries);

  outf->Write();
  std::cout<<"---- Writing tree done"<<std::endl;
  outf->Close();
  std::cout<<outfile<<std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==4) { return evtmatching(argv[1], argv[2], argv[3]); }
  return 1;
}
