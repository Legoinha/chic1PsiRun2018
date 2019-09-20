#!/bin/bash

##
tcounts=(0)
types=(
    "pt"     # 0
)
#
counts=(0 1)
inputmc=(
    /raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
    /raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
)
inputtags=("a" "b")

##
RUN_FIT=${1:-0}
RUN_EFF=${2:-0}
RUN_DRAW=${3:-0}

outputdir=rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_pt15-50_cent090_y0p0-1p6
output=${outputdir##*/}

set -x
[[ $RUN_FIT -eq 1 || $# == 0 ]] && { g++ eff_fit.cc -I"../includes/" $(root-config --libs --cflags) -g -o eff_fit.exe || { rm *.exe 2>/dev/null ; exit 1; } } 
[[ $RUN_EFF -eq 1 || $# == 0 ]] && { g++ eff_eff.cc -I"../includes/" $(root-config --libs --cflags) -g -o eff_eff.exe || { rm *.exe 2>/dev/null ; exit 1; } } 
[[ $RUN_DRAW -eq 1 || $# == 0 ]] && { g++ eff_draw.cc -I"../includes/" $(root-config --libs --cflags) -g -o eff_draw.exe || { rm *.exe 2>/dev/null ; exit 1; } } 
set +x

for i in ${tcounts[@]}
do
    thisoutput=$output/${types[i]}
    [[ $RUN_FIT -eq 1 ]] && {
        ./eff_fit.exe $outputdir/${types[i]}/datamc_fithist.root $thisoutput ${types[i]}
    }
done

[[ $RUN_EFF -eq 1 ]] && {
    cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && mvapref"
    cut=${cut}" && BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13"
    cutgen="1"
    ./eff_eff.exe $outputdir/${types[0]}/eff_fit.root ${inputmc[0]} ${inputmc[1]} "$cut" "$cutgen" $thisoutput 
}

[[ $RUN_DRAW -eq 1 ]] && {
    ./eff_draw.exe $outputdir/${types[0]}/eff_eff.root  $output
}

rm eff_draw.exe 2>/dev/null
rm eff_eff.exe 2>/dev/null
rm eff_fit.exe 2>/dev/null

