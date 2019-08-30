#!/bin/bash

##
tcounts=(0)
types=(
    "pt"     # 0
    "absy"   # 1
)
#
counts=(0 1)
inputmc=(
    /raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
    /raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
)
inputtags=("a" "b")

##
RUN_TREE=${1:-0}
RUN_FIT=${2:-0}
RUN_ACC=${3:-0}

outputdir=rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_pt15-50_cent090_y0p0-1p6

set -x
[[ $RUN_TREE -eq 1 || $# == 0 ]] && { g++ acc_tree.cc $(root-config --libs --cflags) -g -o acc_tree.exe || { rm *.exe 2>/dev/null ; exit 1; } } 
[[ $RUN_FIT -eq 1 || $# == 0 ]] && { g++ acc_fit.cc $(root-config --libs --cflags) -g -o acc_fit.exe || { rm *.exe 2>/dev/null ; exit 1; } } 
[[ $RUN_ACC -eq 1 || $# == 0 ]] && { g++ acc_acc.cc $(root-config --libs --cflags) -g -o acc_acc.exe || { rm *.exe 2>/dev/null ; exit 1; } } 
set +x

[[ $RUN_TREE -eq 1 ]] && {
    for i in ${counts[@]}
    do
        output=${inputmc[i]%%.root}_skimacc.root
        [[ -f $output ]] && continue

        ./acc_tree.exe ${inputmc[i]} $output
    done
}

for i in ${tcounts[@]}
do
    output=${outputdir##*/}/${types[i]}
    [[ $RUN_FIT -eq 1 ]] && {
        ./acc_fit.exe $outputdir/${types[i]}/datamc_fithist.root $output ${types[i]}
    }
    [[ $RUN_ACC -eq 1 ]] && {
        for j in ${counts[@]}
        do
            inputskim=${inputmc[j]%%.root}_skimacc.root
            ./acc_acc.exe $inputskim $outputdir/${types[i]}/acc_fit.root $output ${types[i]} ${inputtags[j]}
        done
    }
done

rm acc_acc.exe 2>/dev/null
rm acc_fit.exe 2>/dev/null
rm acc_tree.exe 2>/dev/null
