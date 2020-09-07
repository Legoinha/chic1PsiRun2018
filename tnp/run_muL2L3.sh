#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

name=trainX_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3_mbin ; ptbins='float ptbins[] = {15., 17., 19., 21., 23., 25., 27., 30., 33., 36., 40., 45., 50.};' ;
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_mutrg_20190808_20200830rmevt_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_pthatweight_trainX_20190808ptdep.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_mutrg_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_pthatweight_trainX_20190808ptdep.root
optcutntuples="ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13"

tmp=$(date +%y%m%d%H%M%S)

RUN_SAVEHIST=${1:-0}
RUN_DRAWHIST=${2:-0}

set -x
g++ getfname.cc -I"../includes/" $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe
echo -e "\e[32;1mcompiling...\e[0m"

cp tnpcc.h tnpcc_tmp.h
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc_tmp.h
cp muL2L3_savehist.cc muL2L3_savehist_tmp_${tmp}.cc
sed -i "s/__CUTINPUT__/${optcutntuples}/g" muL2L3_savehist_tmp_${tmp}.cc

[[ $RUN_SAVEHIST -eq 1 || $# == 0 ]] && { g++ muL2L3_savehist_tmp_${tmp}.cc -I"../includes/" $(root-config --libs --cflags) -g -o muL2L3_savehist_${tmp}.exe || { rm *_${tmp}.* ; exit 1 ; } }
[[ $RUN_DRAWHIST -eq 1 || $# == 0 ]] && { g++ muL2L3_drawhist.cc -I"../includes/" $(root-config --libs --cflags) -g -o muL2L3_drawhist_${tmp}.exe || { rm *_${tmp}.* ; exit 1 ; } }

rm tnpcc_tmp.h
rm muL2L3_savehist_tmp_${tmp}.cc
set +x

[[ $RUN_SAVEHIST -eq 1 ]] && {
    ./muL2L3_savehist_${tmp}.exe $inputmc_a_prompt $name "_a" "noweight" $ptmin $ptmax $centmin $centmax $ymin $ymax
    ./muL2L3_savehist_${tmp}.exe $inputmc_b_prompt $name "_b" "noweight" $ptmin $ptmax $centmin $centmax $ymin $ymax
}
rootdir=rootfiles/$name$kinematic/

[[ $RUN_DRAWHIST -eq 1 ]] && {
    ./muL2L3_drawhist_${tmp}.exe "$rootdir/muL2L3_a.root" $name "_a"
    ./muL2L3_drawhist_${tmp}.exe "$rootdir/muL2L3_b.root" $name "_b"
}

rm muL2L3_savehist_${tmp}.exe 2> /dev/null
rm muL2L3_drawhist_${tmp}.exe 2> /dev/null

