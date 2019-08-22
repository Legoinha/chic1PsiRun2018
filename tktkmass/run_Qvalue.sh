#!/bin/bash

inputmc=(
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root 0"
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root 1"
    "/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXPiPi_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root 2"
)
name=nocut

RUN_FILL=${1:-0}
RUN_DRAW=${2:-0}

g++ checkQvalue.cc $(root-config --libs --cflags) -g -o checkQvalue.exe || exit 1
g++ drawQvalue.cc $(root-config --libs --cflags) -g -o drawQvalue.exe || exit 1

[[ $RUN_FILL -eq 1 ]] && {

    for i in "${inputmc[@]}"
    do
        ./checkQvalue.exe $i $name
    done
}

[[ $RUN_DRAW -eq 1 ]] && ./drawQvalue.exe $name

rm drawQvalue.exe
rm checkQvalue.exe

