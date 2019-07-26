#!/bin/bash

BDTG=0.70
input=/export/d00/scratch/jwang/BntupleRun2018/mva_output/ntmix_20190711_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimBpt10_skimhltBsize_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
input_ss=/export/d00/scratch/jwang/BntupleRun2018/mva_output_old/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
inputmc_a_prompt=/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
inputmc_b_prompt=/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
inputmc_a_nonprompt=/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
inputmc_b_nonprompt=/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root

cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && hiBin >=0 && hiBin < 180"
cut="$cut && mvapref && BDTG>${BDTG}"
cutgen="hiBin >=0 && hiBin < 180"

name=trainX_sideband_tktk0p2_10p0_inf_0-10-1-2-9_officialMC
outputdir=rootfiles/$name
mkdir -p $outputdir

g++ lxydis.cc $(root-config --libs --cflags) -g -o lxydis.exe || exit 1
g++ drawlxydis.cc $(root-config --libs --cflags) -g -o drawlxydis.exe || { rm *.exe ; exit 1 ; }

[[ ${1:-0} -eq 1 ]] && ./lxydis.exe $input $inputmc_a_prompt $inputmc_b_prompt $inputmc_a_nonprompt $inputmc_b_nonprompt "$cut" "$outputdir/lxydis"
[[ ${2:-0} -eq 1 ]] && ./drawlxydis.exe "$outputdir/lxydis" "$outputdir/drawlxydis"

rm drawlxydis.exe
rm lxydis.exe
