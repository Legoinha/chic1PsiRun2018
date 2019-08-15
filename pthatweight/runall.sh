#!/bin/bash

ichan=(0)
###
filelists=(
    "/export/d00/scratch/jwang/BntupleRun2018/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_Pthat-*_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix.root" # 0: nonprompt psi'
    "/export/d00/scratch/jwang/BntupleRun2018/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_Pthat-*_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix.root"  # 1: nonprompt X
    "/export/d00/scratch/jwang/BntupleRun2018/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_Pthat-*_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix.root"    # 2: prompt psi'
    "/export/d00/scratch/jwang/BntupleRun2018/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXPiPi_Pthat-*_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix.root"    # 3: prompt X pipi
    "/export/d00/scratch/jwang/BntupleRun2018/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_Pthat-*_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_ntmix.root"     # 4: prompt X rho
)
##
crosssec=(
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={1.230e+07, 2.868e+06, 8.902e+05, 7.930e+04, 9.631e+03, 0.}; int genSignal[1]={7};' # 0: nonprompt psi'
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={7.622e+06, 1.776e+06, 5.638e+05, 4.973e+04, 6.273e+03, 0.}; int genSignal[1]={7};' # 1: nonprompt X
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={1.296e+05, 1.954e+04, 4.672e+03, 3.522e+02, 4.040e+01, 0.}; int genSignal[1]={7};' # 2: prompt psi'
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={2.617e+04, 1.987e+04, 6.830e+03, 3.076e+02, 2.775e+01, 0.}; int genSignal[1]={7};' # 3: prompt X pipi
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={2.653e+04, 1.979e+04, 6.670e+03, 3.053e+02, 2.768e+01, 0.}; int genSignal[1]={7};' # 4: prompt X rho
)

##
for ii in ${ichan[@]}
do
##
    tmp=$(date +%y%m%d%H%M%S)
    sed '1i'"${crosssec[$ii]}" weighPurePthat.C > weighPurePthat_${tmp}.C

    g++ addbranch.C $(root-config --cflags --libs) -g -o addbranch_${tmp}.exe || { rm weighPurePthat_${tmp}.C ; exit 1 ; }
    g++ weighPurePthat_${tmp}.C $(root-config --cflags --libs) -g -o weighPurePthat_${tmp}.exe || { rm weighPurePthat_${tmp}.C ; rm addbranch_${tmp}.exe ; exit 1 ; }

    filelist=${filelists[ii]}
    echo "=========== add sample pthat cut value >>>>"
    mergelist=
    for ifile in `echo $filelist`
    do
        ifilecp=${ifile%%.root}
        ifilecp=${ifilecp}_addSamplePthat.root
        pthatcut=${ifile##*Pthat-} ; pthatcut=${pthatcut%%_*.root}
        echo "----------"
        echo "input:  $ifile"
        echo "output: $ifilecp"
        echo "pthatcut value: $pthatcut"
        [[ $ifile == $ifilecp ]] && { echo "invalid input for ./addbranch.exe" ; continue ; }
        [[ ${1:-0} -eq 1 ]] && {
            rsync --progress $ifile $ifilecp
            set -x
            yes y | ./addbranch_${tmp}.exe $ifilecp hiEvtAnalyzer/HiTree sample $pthatcut float
            set +x
        }
        mergelist="$mergelist "$ifilecp
    done

    echo

##
    echo "=========== merge >>>>"
    mergesuffix=${filelist##*'Pthat-*'} ; mergesuffix=${mergesuffix%%.root}
    mergepreffix=${filelist%%_Pthat*}
    mergeoutput=$mergepreffix${mergesuffix}_addSamplePthat_noweight.root
    [[ ${2:-0} -eq 1 ]] && {
        set -x
        hadd $mergeoutput $mergelist
        set +x
    }

    echo

##
    echo "=========== weight >>>>"
    weightoutput=$mergepreffix${mergesuffix}_addSamplePthat_pthatweight.root
    [[ ${3:-0} -eq 1 ]] && { 
        rsync --progress "$mergeoutput" "$weightoutput"
        set -x
        ./weighPurePthat_${tmp}.exe "$mergeoutput" "$weightoutput" 
        set +x
    }

    rm weighPurePthat_${tmp}.exe
    rm weighPurePthat_${tmp}.C
    rm addbranch_${tmp}.exe
done

