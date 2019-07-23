#!/bin/bash

inputpars=(
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-10_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 10."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-15_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 15."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-30_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 30."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-50_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 50."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-5_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 5."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-10_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 10."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-15_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 15."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-30_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 30."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-50_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 50."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-5_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 5."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-10_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 10."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-15_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 15."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-30_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 30."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-50_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 50."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-5_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 5."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-10_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 10."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-15_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 15."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-30_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 30."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-50_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 50."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-5_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 5."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-10_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 10."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-15_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 15."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-30_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 30."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-50_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 50."
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-5_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix_addSamplePthat.root hiEvtAnalyzer/HiTree sample 5."
)

g++ addbranch.C $(root-config --libs --cflags) -g -o addbranch.exe || exit 1

for ii in "${inputpars[@]}"
do
    echo $ii
    pars=($ii)
    echo $pars
    cpfile=${pars[0]%%_addSamplePthat.root}.root
    set -x
    cp $cpfile ${pars[0]}

    ./addbranch.exe $ii
done

rm addbranch.exe
set +x
