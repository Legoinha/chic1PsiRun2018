#!/bin/bash

skim=0
ifs=({7..31}) # input directory
jns=(1) # tree

## ifs
inputdirs=(
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327327_skimBpt10/      # 0
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327400_327564_v2_skimBpt10/   # 1
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326381_326577_skimBpt10/             # 2
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326580_326855_skimBpt10/             # 3
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_20190513_HIDoubleMuon_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_326856_327078_skimBpt10/             # 4
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327327_skimBpt10/ # 5
    /export/d00/scratch/jwang/BntupleRun2018/ntmix_20190711_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327400_327564_skimBpt10/ # 6
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-10_1033p1_official_pt6tkpt0p9dls0/ # 7
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-15_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-30_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-50_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-5_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-10_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-15_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-30_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-50_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-5_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-10_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-15_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-30_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-50_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-5_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-10_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-15_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-30_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-50_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_Pthat-5_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-10_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-15_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-30_1033p1_official_pt6tkpt0p9dls0/
    # /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-50_1033p1_official_pt6tkpt0p9dls0/
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_Pthat-5_1033p1_official_pt6tkpt0p9dls0/ # 31
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
    [[ $i -lt ${#inputdirs[@]} ]] || break
    inputdir=${inputdirs[i]}
    IFS='/'; subdir=($inputdir); unset IFS;
    request=${subdir[${#subdir[@]}-1]}
    primedir=${inputdir%%${request}*}

    [[ ! -d $inputdir ]] && continue

    ## ======================================== #

    filelist=filelist_${request}.txt
    [[ -f $filelist ]] && {
        # echo "error: filelist $filelist exits. "
        # echo "remove filelist? (y/n):"
        # rewrite=
        # while [[ $rewrite != 'y' && $rewrite != 'n' ]]
        # do
        #     read rewrite
        #     if [[ $rewrite == 'y' ]] ; then { rm $filelist ; } ;
        #     elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; rm merge_${tmp}.exe ; continue ; } ;
        #     else { echo "please input y/n" ; } ; fi ;
        # done
        rm $filelist
    } 

    ls $inputdir/*.root -d > $filelist

    for j in ${jns[@]}
    do
        ntname=${nts[j]}
        set -x
        output=${primedir}/${request}_skimhltBsize_${ntname}.root
        set +x
        willrun=1
        [[ -f $output ]] && {
            echo "error: output $output exits. "
            echo "remove output? (y/n):"
            rewrite=
            while [[ $rewrite != 'y' && $rewrite != 'n' ]]
            do
                read rewrite
                if [[ $rewrite == 'y' ]] ; then { echo "$output removed" ; rm $output ; } ;
                elif [[ $rewrite == 'n' ]] ; then { echo "please change output file name" ; willrun=0 ; } ;
                else { echo "please input y/n" ; } ; fi ;
            done
        }

        [[ $willrun -eq 0 ]] && continue
        [[ ${1:-0} -eq 1 ]] && { ./merge_${tmp}.exe $output $filelist $ntname $skim ; }
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
