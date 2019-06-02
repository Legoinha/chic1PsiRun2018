#!/bin/bash

ifs=(6) # input directory
jns=(1) # tree

## ifs
inputdirs=(
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327327/      # 0
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327400_327564_v2/   # 1
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326381_326577/             # 2
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326580_326855/             # 3
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326856_327078/             # 4
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327327/ # 5
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327400_327564/ # 6
)

########################################
## >>> do not change lines below >>>  ##
########################################

########################################
nts=( "ntKp" "ntmix" "ntphi" "ntKstar" "empty" "root")
#jns  0      1       2       3         4       5
########################################

tmp=$(date +%y%m%d%H%M%S)
cp merge.C merge_${tmp}.C
g++ merge_${tmp}.C $(root-config --libs --cflags) -g -o merge_${tmp}.exe || exit 1

for i in ${ifs[@]}
do
    inputdir=${inputdirs[i]}
    IFS='/'; subdir=($inputdir); unset IFS;
    request=${subdir[${#subdir[@]}-1]}
    primedir=${inputdir%%crab_*}

##
    filelist=filelist_${request}.txt
    [[ -f $filelist ]] && {
        echo "error: filelist $filelist exits. "
        echo "remove filelist? (y/n):"
        rewrite=
        while [[ $rewrite != 'y' && $rewrite != 'n' ]]
        do
            read rewrite
            if [[ $rewrite == 'y' ]] ; then { rm $filelist ; } ;
            elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; rm merge_${tmp}.exe ; exit 2 ; } ;
            else { echo "please input y/n" ; } ; fi ;
        done
    } 

    ls $inputdir/*.root -d > $filelist

    for j in ${jns[@]}
    do
        ntname=${nts[j]}
        set -x
        output=${primedir}/${request}_skimhltBsize_${ntname}.root
        set +x
        [[ -f $output ]] && {
            echo "error: output $output exits. "
            echo "remove output? (y/n):"
            rewrite=
            while [[ $rewrite != 'y' && $rewrite != 'n' ]]
            do
                read rewrite
                if [[ $rewrite == 'y' ]] ; then { echo "$output removed" ; rm $output ; } ;
                elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; rm merge_${tmp}.exe ; exit 2 ; } ;
                else { echo "please input y/n" ; } ; fi ;
            done
        }

        [[ ${1:-0} -eq 1 ]] && { ./merge_${tmp}.exe $output $filelist $ntname ; }
    done
done

rm merge_${tmp}.exe
rm merge_${tmp}.C

##

    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1_1031_NoJSON/ # 0
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON/ # 1
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON_Run327527_327564/ # 2
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON_ToComplete/ # 3
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1_1031_NoJSON_ToComplete/ # 4
