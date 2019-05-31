#!/bin/bash

jns=(0) # tree

inputfiles=(
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v2_1031_NoJSON_ToComplete_skimhltBsize_ntKstar.root
    HiForestAOD_971.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326580_326855/HiForestAOD_238.root
)
output=/export/d00/scratch/jwang/BntupleRun2018/HiForestAOD_971_238.root

## >>> do not change lines below >>>
echo $output

##
nts=("ntKp" "ntmix" "ntphi" "ntKstar")
#    (0)    (1)     (2)     (3) 
##
[[ -f $output ]] && {
    echo "error: output $output exits. "
    echo "remove output? (y/n):"
    rewrite=
    while [[ $rewrite != 'y' && $rewrite != 'n' ]]
    do
        read rewrite
        if [[ $rewrite == 'y' ]] ; then { echo "$output removed" ; rm $output ; } ;
        elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; exit 2 ; } ; 
        else { echo "please input y/n" ; } ; fi ;
    done
}

filelist=${output##*/}
filelist=filelist_${filelist%%.*}.txt

[[ -f $filelist ]] && {
    echo "error: filelist $filelist exits. "
    echo "remove filelist? (y/n):"
    rewrite=
    while [[ $rewrite != 'y' && $rewrite != 'n' ]]
    do
        read rewrite
        if [[ $rewrite == 'y' ]] ; then { rm $filelist ; touch $filelist ; } ;
        elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; exit 2 ; } ; 
        else { echo "please input y/n" ; } ; fi ;
    done
} || touch $filelist 

for i in ${inputfiles[@]} ; do echo $i >> $filelist ; done ;

thisone=$(date +%y%m%d%H%M%S)
cp merge.C merge_${thisone}.C

g++ merge_${thisone}.C $(root-config --libs --cflags) -g -o merge_${thisone}.exe || exit 1
for j in ${jns[@]}
do
    ntname=${nts[j]}
    [[ ${1:-0} -eq 1 ]] && { ./merge_${thisone}.exe $output $filelist $ntname ; }
done
rm merge_${thisone}.exe
rm merge_${thisone}.C

##
