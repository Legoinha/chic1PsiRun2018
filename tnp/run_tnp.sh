#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

# name=trainX_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_BDTQvalue_1bin ; ptbins='float ptbins[] = {15., 50.};' ;
name=trainX_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_BDTQvalue_mbin ; ptbins='float ptbins[] = {15., 17., 19., 21., 23., 25., 27., 30., 33., 36., 40., 45., 50.};' ;
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
optcutntuples="ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13"

tmp=$(date +%y%m%d%H%M%S)

RUN_CONVER=${1:-0}
RUN_DRAW=${2:-0}
RUN_DRAWRATIO=${3:-0}
RUN_PTWEIGHT=${4:-0}

set -x
g++ getfname.cc -I"../includes/" $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe
echo -e "\e[32;1mcompiling...\e[0m"

cp tnpcc.h tnpcc_tmp.h
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc_tmp.h
cp tnp_converter.cc tnp_converter_tmp_${tmp}.cc
sed -i "s/__CUTINPUT__/${optcutntuples}/g" tnp_converter_tmp_${tmp}.cc

[[ $RUN_CONVER -eq 1 || $# == 0 ]] && { g++ tnp_converter_tmp_${tmp}.cc -I"../includes/" $(root-config --libs --cflags) -g -o tnp_converter_${tmp}.exe || { rm *_${tmp}.* ; exit 1 ; } }
[[ $RUN_DRAW -eq 1 || $# == 0 ]] && { g++ draw_tnp.cc -I"../includes/" $(root-config --libs --cflags) -g -o draw_tnp_${tmp}.exe || { rm *_${tmp}.* ; exit 1 ;  } ; 
    g++ draw_mu.cc -I"../includes/" $(root-config --libs --cflags) -g -o draw_mu_${tmp}.exe || { rm *_${tmp}.* ; exit 1 ;  } ; }
[[ $RUN_DRAWRATIO -eq 1 || $# == 0 ]] && { g++ draw_tnp_ratio.cc -I"../includes/" $(root-config --libs --cflags) -g -o draw_tnp_ratio_${tmp}.exe || { rm *_${tmp}.* ; exit 1 ;  } ; }

rm tnpcc_tmp.h
rm tnp_converter_tmp_${tmp}.cc
set +x

[[ $RUN_CONVER -eq 1 ]] && {
    ./tnp_converter_${tmp}.exe $inputmc_a_prompt $name "_a" "noweight" $ptmin $ptmax $centmin $centmax $ymin $ymax
    ./tnp_converter_${tmp}.exe $inputmc_b_prompt $name "_b" "noweight" $ptmin $ptmax $centmin $centmax $ymin $ymax
    [[ $RUN_PTWEIGHT -eq 1 ]] && {
        ./tnp_converter_${tmp}.exe $inputmc_a_prompt $name "_a" "eff_fit.root" $ptmin $ptmax $centmin $centmax $ymin $ymax
        ./tnp_converter_${tmp}.exe $inputmc_b_prompt $name "_b" "eff_fit.root" $ptmin $ptmax $centmin $centmax $ymin $ymax
    }
}
rootdir=rootfiles/$name$kinematic/

[[ $RUN_DRAW -eq 1 ]] && {
    ./draw_tnp_${tmp}.exe "$rootdir/tnp_a.root" $name "_a"
    ./draw_tnp_${tmp}.exe "$rootdir/tnp_b.root" $name "_b"
    [[ $RUN_PTWEIGHT -eq 1 ]] && {
        ./draw_tnp_${tmp}.exe "$rootdir/funs/fun-1/tnp_a.root" "$rootdir/funs/fun-2/tnp_a.root" "$rootdir/funs/fun-3/tnp_a.root" "$rootdir/funs/fun-4/tnp_a.root" "$rootdir/funs/fun-5/tnp_a.root" $name "_a"
        ./draw_tnp_${tmp}.exe "$rootdir/funs/fun-1/tnp_b.root" "$rootdir/funs/fun-2/tnp_b.root" "$rootdir/funs/fun-3/tnp_b.root" "$rootdir/funs/fun-4/tnp_b.root" "$rootdir/funs/fun-5/tnp_b.root" $name "_b"
    }
    ./draw_mu_${tmp}.exe "$rootdir/tnp_a.root" $name "_a"
    ./draw_mu_${tmp}.exe "$rootdir/tnp_b.root" $name "_b"
}

[[ $RUN_DRAWRATIO -eq 1 ]] && {
    ./draw_tnp_ratio_${tmp}.exe "$rootdir/drawtnp_a.root" "$rootdir/drawtnp_b.root" $name
}

rm draw_tnp_ratio_${tmp}.exe 2> /dev/null
rm draw_tnp_${tmp}.exe 2> /dev/null
rm draw_mu_${tmp}.exe 2> /dev/null
rm tnp_converter_${tmp}.exe 2> /dev/null

