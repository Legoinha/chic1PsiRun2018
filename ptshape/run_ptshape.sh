#!/bin/bash

BDTG=0.70
inputdirname=HEPData-ins1495026-v1-csv

idxi=(0 1)
input=(
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root promptPsi"
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root promptX"
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root nonpromptPsi"
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root nonpromptX"
)

cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && hiBin >=0 && hiBin < 180"
cut="$cut && mvapref && BDTG>${BDTG}"
cutgen="hiBin >=0 && hiBin < 180"

g++ fillptshape.cc $(root-config --libs --cflags) -g -o fillptshape.exe || exit 1
g++ fitptshape.cc $(root-config --libs --cflags) -g -o fitptshape.exe || exit 1
g++ compeff.cc $(root-config --libs --cflags) -g -o compeff.exe || exit 1
g++ drawcompeff.cc $(root-config --libs --cflags) -g -o drawcompeff.exe || exit 1

for ii in "${idxi[@]}"
do
    iis=(${input[ii]})
    outputname="rootfiles/$inputdirname/${iis[1]}.root"
    [[ ${1:-0} -eq 1 ]] && {
        set -x
        ./fillptshape.exe ${input[ii]} $outputname $inputdirname
        set +x
    }
    [[ ${2:-0} -eq 1 ]] && {
        set -x
        ./fitptshape.exe $outputname ${iis[1]} $inputdirname
        set +x
    }
    [[ ${3:-0} -eq 1 ]] && {
        set -x
        ./compeff.exe ${input[ii]} "$cut" "$cutgen" "rootfiles/$inputdirname/compeff_${iis[1]}.root" $inputdirname
        set +x
    }
done

[[ ${4:-0} -eq 1 ]] && {
    set -x
    ./drawcompeff.exe "rootfiles/$inputdirname/compeff_promptPsi.root" "rootfiles/$inputdirname/compeff_promptX.root" "plots/$inputdirname/ccompeff.pdf"
    set +x
}

rm drawcompeff.exe
rm compeff.exe
rm fitptshape.exe
rm fillptshape.exe

