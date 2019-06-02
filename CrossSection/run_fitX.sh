#!/bin/bash

input=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
input_ss=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
inputmc_a=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
inputmc_b=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root

# cut="1" ; name="test" ;
cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && hiBin >=0 && hiBin < 180"
# cut=$cut" && Bpt>15 && mvapref && BDTG>0.76 && BsvpvDistance/BsvpvDisErr>0.8 && TMath::Abs(By)<1.5" ; name="dls0p8" ;
cut=$cut" && Bpt>15 && mvapref && BDTG>0.76 && TMath::Abs(By)<1.5" ; name="nodls" ;
cutgen="Gpt>15&&TMath::Abs(Gy)<1.5&&GisSignal==7&&GcollisionId==0"

mkdir -p plots

# g++ fitX.C $(root-config --libs --cflags) -g -o fitX.exe || exit 1
# ./fitX.exe $input "$cut"
# rm fitX.exe

g++ fitXmc.C $(root-config --libs --cflags) -g -o fitXmc.exe || exit 1
g++ MCefficiency.C $(root-config --libs --cflags) -g -o MCefficiency.exe || exit 1


[[ ${1:-0} -eq 1 ]] && ./fitXmc.exe $input $input_ss $inputmc_a $inputmc_b "${cut}" $name
[[ ${2:-0} -eq 1 ]] && ./MCefficiency.exe $inputmc_a $inputmc_b "${cut}" "$cutgen" $name


rm MCefficiency.exe
rm fitXmc.exe

