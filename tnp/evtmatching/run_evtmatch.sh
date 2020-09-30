#!/bin/bash

# inputname=/export/d00/scratch/jwang/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
# inputname_mu=/export/d00/scratch/jwang/BntupleRun2018/crab_MuonOnly_20200825_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_addSamplePthat_noweight.root
# outputname=/export/d00/scratch/jwang/BntupleRun2018/mu_evtmatch_20200825_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0.root

inputname=/export/d00/scratch/jwang/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputname_mu=/export/d00/scratch/jwang/BntupleRun2018/crab_MuonOnly_20200825_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_addSamplePthat_noweight.root
outputname=/export/d00/scratch/jwang/BntupleRun2018/mu_evtmatch_20200825_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0.root

g++ evtmatching.C $(root-config --libs --cflags) -I$HOME -g -o evtmatching.exe || exit 1

[[ $1 -eq 1 ]] && { ./evtmatching.exe $inputname $inputname_mu $outputname ; }

rm evtmatching.exe