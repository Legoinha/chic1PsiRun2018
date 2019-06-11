#!/bin/bash

arguments=(0 1 2 3 4 5)
mvaval=0.76

input=(
 "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root                 data          $mvaval rootfiles/root_data"          # 0
 "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root samesign      $mvaval rootfiles/root_samesign"      # 1
 "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root                 psi_prompt    $mvaval rootfiles/root_psi_prompt"    # 2
 "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root              psi_nonprompt $mvaval rootfiles/root_psi_nonprompt" # 3
 "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root                  x_prompt      $mvaval rootfiles/root_x_prompt"      # 4
 "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root               x_nonprompt   $mvaval rootfiles/root_x_nonprompt"   # 5
)
g++ lxydis.cc $(root-config --libs --cflags) -g -o lxydis.exe || exit 1
g++ drawlxydis.cc $(root-config --libs --cflags) -g -o drawlxydis.exe || exit 1
g++ calcfprompt.cc $(root-config --libs --cflags) -g -o calcfprompt.exe || exit 1

mkdir -p rootfiles plots

for i in ${arguments[@]}
do
    [[ ${1:-0} -eq 1 ]] && ./lxydis.exe ${input[i]}
done

[[ ${2:-0} -eq 1 ]] && ./drawlxydis.exe rootfiles/root_data.root rootfiles/root_samesign.root rootfiles/root_psi_prompt.root rootfiles/root_psi_nonprompt.root rootfiles/root_x_prompt.root rootfiles/root_x_nonprompt.root

[[ ${3:-0} -eq 1 ]] && ./calcfprompt.exe rootfilesCrossSection/root_fitXmc_nodls.root rootfilesCrossSection/root_fitXmc_lxygt0p1.root rootfiles/root_lxydis.root

rm calcfprompt.exe
rm drawlxydis.exe
rm lxydis.exe

