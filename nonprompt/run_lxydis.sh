#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

optcuts="&& BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13"
tags="BDTQvalue"
lxyvar="lxy" ## lxy, lxyz

input=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190806_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_skim_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root

##

RUN_SAVEHIST=${1:-0}
RUN_FITHIST=${2:-0}
RUN_DRAWHIST=${3:-0}

tmp=$(date +%y%m%d%H%M%S)

##
g++ getfname.cc $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe
echo -e "\e[32;1mcompiling...\e[0m"
[[ $RUN_SAVEHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fprompt_savehist.C
    g++ fprompt_savehist.C $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o fprompt_savehist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; } 
}
[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fprompt_fithist.C
    g++ fprompt_fithist.C $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o fprompt_fithist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}
[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fprompt_drawhist.C
    g++ fprompt_drawhist.C $(root-config --libs --cflags) -g -o fprompt_drawhist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

##
echo
cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
cut="$cut && mvapref ${optcuts}"

name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${tags}
echo -e "----------------------------------------"
echo -e "==> File directory: \e[4m$name/$lxyvar\e[0m"
echo -e "----------------------------------------"

[[ $RUN_SAVEHIST -eq 1 ]] && {
    ./fprompt_savehist_${tmp}.exe $input $inputmc_a_prompt $inputmc_b_prompt $inputmc_a_nonprompt $inputmc_b_nonprompt "$cut" "$name" "$lxyvar" $ptmin $ptmax $centmin $centmax $ymin $ymax
}

rootdir=rootfiles/$name$kinematic/

[[ $RUN_FITHIST -eq 1 ]] && {
    ./fprompt_fithist_${tmp}.exe "$rootdir/$lxyvar/fprompt_savehist.root" $name "$lxyvar"
}
[[ $RUN_DRAWHIST -eq 1 ]] && ./fprompt_drawhist_${tmp}.exe "$name$kinematic"

rm *_${tmp}.exe 2> /dev/null

