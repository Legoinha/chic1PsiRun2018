#!/bin/bash

name=trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_finebin ; ptbins='float ptbins[] = {15., 20., 25., 30., 35., 40., 45., 50.};' ;
# name=trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_1bin ; ptbins='float ptbins[] = {15., 50.};' ;
idxi=(0 1)
input=(
    /export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
    /export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root
)

cp tnpcc.h tnpcc_tmp.h
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc.h

g++ tnp_converter.cc $(root-config --libs --cflags) -g -o tnp_converter.exe || exit 1
g++ draw_tnp.cc $(root-config --libs --cflags) -g -o draw_tnp.exe || { rm *.exe ; exit 1 ;  }
g++ draw_tnp_ratio.cc $(root-config --libs --cflags) -g -o draw_tnp_ratio.exe || { rm *.exe ; exit 1 ;  }

[[ ${1:-0} -eq 1 ]] && {
    for ii in ${idxi[@]}
    do
        echo ${input[ii]}
        ./tnp_converter.exe ${input[ii]} $name
    done
}

[[ ${2:-0} -eq 1 ]] && {
    for ii in `echo rootfiles/$name/tnp_*.root`
    do
        echo $ii
        ./draw_tnp.exe "$ii" $name
    done
}

inputforratio=
for ii in ${idxi[@]}
do
    inputname=${input[ii]}
    outputname=rootfiles/$name/draw_tnp_${inputname##*/}
    # echo $outputname
    inputforratio="$inputforratio $outputname"
done
# echo $inputforratio

[[ ${3:-0} -eq 1 ]] && {
    ./draw_tnp_ratio.exe $inputforratio $name
}

rm draw_tnp_ratio.exe
rm draw_tnp.exe
rm tnp_converter.exe

mv tnpcc_tmp.h tnpcc.h