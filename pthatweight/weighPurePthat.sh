#!/bin/bash

###
inputfiles=(
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_noweight.root"
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_noweight.root"
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_noweight.root"
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_noweight.root"
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_noweight.root"
)
# outputfile="/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_pthatweight.root"
# chan=3 # ([0: Bs]; [1: prompt psi']; [2: prompt X->jpsi + rho]; [3: Bu]; [4: prompt psi']; [5: prompt X->jpsi + rho])

##
crosssec=(
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={1.230e+07, 2.868e+06, 8.902e+05, 7.930e+04, 9.631e+03, 0.}; int genSignal[1]={7};' # nonprompt psi'
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={7.622e+06, 1.776e+06, 5.638e+05, 4.973e+04, 6.273e+03, 0.}; int genSignal[1]={7};' # nonprompt X
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={1.296e+05, 1.954e+04, 4.672e+03, 3.522e+02, 4.040e+01, 0.}; int genSignal[1]={7};' # prompt psi'
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={2.617e+04, 1.987e+04, 6.830e+03, 3.076e+02, 2.775e+01, 0.}; int genSignal[1]={7};' # prompt X pipi
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={2.653e+04, 1.979e+04, 6.670e+03, 3.053e+02, 2.768e+01, 0.}; int genSignal[1]={7};' # prompt X rho
)

i=0
for ii in ${inputfiles[@]}
do
    tmp=$(date +%y%m%d%H%M%S)
    sed '1i'"${crosssec[$i]}" weighPurePthat.C > weighPurePthat_${tmp}.C

    g++ weighPurePthat_${tmp}.C $(root-config --cflags --libs) -g -o weighPurePthat_${tmp}.exe || { rm weighPurePthat_${tmp}.C ; exit 1 ; }

    [[ ${1:-0} -eq 1 ]] && { 
        outputfile=${ii%_noweight.root}_pthatweight.root
        rsync --progress "$ii" "$outputfile"
        echo $ii $outputfile
        set -x
        ./weighPurePthat_${tmp}.exe "$ii" "$outputfile" 
        set +x
    }
    rm weighPurePthat_${tmp}.exe
    rm weighPurePthat_${tmp}.C
    i=$((i+1))
done

