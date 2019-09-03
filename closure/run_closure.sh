#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

counts=(1)
optcuts=(
    "&& BDT > 0.07 && (Bmass-3.096916-Btktkmass) < 0.13"
)
tags=(
    "BDTQvalue"
)

inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root

RUN_SAVEHIST=${1:-0}
RUN_FITHIST=${2:-0}
RUN_DRAWHIST=${3:-0}

tmp=$(date +%y%m%d%H%M%S)

g++ getfname.cc $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || exit 1
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe

echo -e "\e[32;1mcompiling...\e[0m"

[[ $RUN_SAVEHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitX_flatten.C
    g++ fitX_flatten.C $(root-config --libs --cflags) -g -o fitX_flatten_${tmp}.exe || { rm *_${tmp}.exe 2>/dev/null ; exit 1 ; }
    echo " -- "closure_savehist.cc
    g++ closure_savehist.cc $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o closure_savehist_${tmp}.exe || { rm *_${tmp}.exe 2>/dev/null ; exit 1 ; }
}

[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "closure_fithist.cc
    g++ closure_fithist.cc $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o closure_fithist_${tmp}.exe || { rm *_${tmp}.exe 2>/dev/null ; exit 1 ; }
}

[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "closure_drawhist.cc
    g++ closure_drawhist.cc $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o closure_drawhist_${tmp}.exe || { rm *_${tmp}.exe 2>/dev/null ; exit 1 ; }
}

[[ $RUN_SAVEHIST -eq 1 ]] && {
    flatinputmc_a_prompt=${inputmc_a_prompt%%.root}_flatten.root
    [[ -f $flatinputmc_a_prompt ]] || ./fitX_flatten_${tmp}.exe $inputmc_a_prompt $flatinputmc_a_prompt
    flatinputmc_b_prompt=${inputmc_b_prompt%%.root}_flatten.root
    [[ -f $flatinputmc_b_prompt ]] || ./fitX_flatten_${tmp}.exe $inputmc_b_prompt $flatinputmc_b_prompt
}

echo
for count in ${counts[@]}
do
    cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
    cut="$cut && mvapref ${optcuts[count]}"
    cutgen="1"

    ##
    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${tags[count]}
    rootdir=rootfiles/$name$kinematic/
    echo -e "----------------------------------------"
    echo -e "==> File directory: \e[4m$rootdir\e[0m"
    echo -e "----------------------------------------"

    [[ $RUN_SAVEHIST -eq 1 ]] && {
        ./closure_savehist_${tmp}.exe $inputmc_a_prompt $inputmc_b_prompt "$cut" "$cutgen" "$name" $ptmin $ptmax $centmin $centmax $ymin $ymax
    }
    [[ $RUN_FITHIST -eq 1 ]] && {
        ./closure_fithist_${tmp}.exe $rootdir/closure_savehist.root $name$kinematic
    }
    [[ $RUN_DRAWHIST -eq 1 ]] && {
        ./closure_drawhist_${tmp}.exe $rootdir/closure_fithist.root $name$kinematic
    }

done

rm *_${tmp}.exe 2>/dev/null
