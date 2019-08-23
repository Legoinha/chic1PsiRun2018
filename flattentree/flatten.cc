#include "tree_convert.h"
#include <iostream>

int main()
{
  tree_flatten ntf("/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190806_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_skimBpt10_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root");
  ntf.flatten();
  TFile* outf = new TFile("outtest.root", "recreate");
  outf->cd();
  ntf.outnt->Write();
  std::cout<<ntf.outnt->GetEntries()<<std::endl;;
  outf->Close();
  return 0;
}
