#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

name=trainX_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_BDTQvalue_1bin ; ptbins='float ptbins[] = {15., 50.};' ;
# name=trainX_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_BDTQvalue_mbin ; ptbins='float ptbins[] = {15., 20., 25., 30., 35., 40., 45., 50.};' ;
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
optcutntuples="ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13"

tmp=$(date +%y%m%d%H%M%S)

RUN_CONVER=${1:-0}
RUN_DRAW=${2:-0}
RUN_DRAWRATIO=${3:-0}

g++ getfname.cc $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe
echo -e "\e[32;1mcompiling...\e[0m"

cp tnpcc.h tnpcc_tmp.h
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc_tmp.h
cp tnp_converter.cc tnp_converter_tmp_${tmp}.cc
sed -i "s/__CUTINPUT__/${optcutntuples}/g" tnp_converter_tmp_${tmp}.cc

[[ $RUN_CONVER -eq 1 || $# == 0 ]] && { g++ tnp_converter_tmp_${tmp}.cc $(root-config --libs --cflags) -g -o tnp_converter_${tmp}.exe || exit 1 ; }
[[ $RUN_DRAW -eq 1 || $# == 0 ]] && { g++ draw_tnp.cc $(root-config --libs --cflags) -g -o draw_tnp_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ;  } ; }
[[ $RUN_DRAWRATIO -eq 1 || $# == 0 ]] && { g++ draw_tnp_ratio.cc $(root-config --libs --cflags) -g -o draw_tnp_ratio_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ;  } ; }

rm tnpcc_tmp.h
rm tnp_converter_tmp_${tmp}.cc

[[ $RUN_CONVER -eq 1 ]] && {
    ./tnp_converter_${tmp}.exe $inputmc_a_prompt $name "_a" $ptmin $ptmax $centmin $centmax $ymin $ymax
    ./tnp_converter_${tmp}.exe $inputmc_b_prompt $name "_b" $ptmin $ptmax $centmin $centmax $ymin $ymax
}
rootdir=rootfiles/$name$kinematic/

[[ $RUN_DRAW -eq 1 ]] && {
    ./draw_tnp_${tmp}.exe "$rootdir/tnp_a.root" $name "_a"
    ./draw_tnp_${tmp}.exe "$rootdir/tnp_b.root" $name "_b"
}

[[ $RUN_DRAWRATIO -eq 1 ]] && {
    ./draw_tnp_ratio_${tmp}.exe "$rootdir/drawtnp_a.root" "$rootdir/drawtnp_b.root" $name
}

rm draw_tnp_ratio_${tmp}.exe 2> /dev/null
rm draw_tnp_${tmp}.exe 2> /dev/null
rm tnp_converter_${tmp}.exe 2> /dev/null

