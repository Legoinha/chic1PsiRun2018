#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

inputdirname=HEPData-ins1495026-v1

idxi=(0 1 2 3)
input=(
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root promptPsi"
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root promptX"
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root nonpromptPsi"
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root nonpromptX"
)

cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
cut="$cut && mvapref && BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13"
cutgen="1"

##
RUN_FILLPT=${1:-0}
RUN_FITPT=${2:-0}
RUN_EFF=${3:-0}
RUN_DRAW=${4:-0}

echo -- fillptshape.cc
[[ $RUN_FILLPT -eq 1 || $# == 0 ]] && { g++ fillptshape.cc $(root-config --libs --cflags) -g -o fillptshape.exe || { rm *.exe 2>/dev/null ; exit 1 ; } ; }
echo -- fitptshape.cc
[[ $RUN_FITPT -eq 1 || $# == 0 ]] && { g++ fitptshape.cc $(root-config --libs --cflags) -g -o fitptshape.exe || { rm *.exe 2>/dev/null ; exit 1 ; } ; }
echo -- compeff.cc
[[ $RUN_EFF -eq 1 || $# == 0 ]] && { g++ compeff.cc $(root-config --libs --cflags) -lRooFitCore -lRooFit -lRooStats -g -o compeff.exe || { rm *.exe 2>/dev/null ; exit 1 ; } ; }
echo -- drawcompeff.cc
[[ $RUN_DRAW -eq 1 || $# == 0 ]] && { g++ drawcompeff.cc $(root-config --libs --cflags) -g -o drawcompeff.exe || { rm *.exe 2>/dev/null ; exit 1 ; } ; }

for ii in "${idxi[@]}"
do
    iis=(${input[ii]})
    outputname="rootfiles/$inputdirname/${iis[1]}.root"
    [[ $RUN_FILLPT -eq 1 ]] && {
        set -x
        ./fillptshape.exe ${input[ii]} $outputname $inputdirname
        set +x
    }
    [[ $RUN_FITPT -eq 1 ]] && {
        set -x
        ./fitptshape.exe $outputname ${iis[1]} $inputdirname
        set +x
    }
    [[ $RUN_EFF -eq 1 ]] && {
        set -x
        ./compeff.exe ${input[ii]} "$cut" "$cutgen" "rootfiles/$inputdirname/compeff_${iis[1]}.root" $inputdirname $ptmin $ptmax $centmin $centmax $ymin $ymax
        set +x
    }
done

[[ $RUN_DRAW -eq 1 ]] && {
    set -x
    ./drawcompeff.exe "rootfiles/$inputdirname/compeff_promptPsi.root" "rootfiles/$inputdirname/compeff_promptX.root" "plots/$inputdirname/ccompeff.pdf"
    set +x
}

rm drawcompeff.exe 2>/dev/null
rm compeff.exe 2>/dev/null
rm fitptshape.exe 2>/dev/null
rm fillptshape.exe 2>/dev/null
