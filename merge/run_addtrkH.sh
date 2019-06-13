#!/bin/bash

inputs=(
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntmix_Xpt10.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root
)

g++ addbranch.C $(root-config --libs --cflags) -g -o addbranch.exe || exit 1

for ii in ${inputs[@]}
do
    yes y | ./addbranch.exe $ii Bfinder/ntmix
done

rm addbranch.exe