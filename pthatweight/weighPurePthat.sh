#!/bin/bash

###
inputfile="/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_noweight.root"
outputfile="/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_pthatweight.root"
chan=3 # ([0: Bs]; [1: prompt psi']; [2: prompt X->jpsi + rho]; [3: Bu]; [4: prompt psi']; [5: prompt X->jpsi + rho])

##
crosssec=(
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={6.123e+06, 2.864e+06, 9.652e+05, 1.299e+05, 1.895e+04, 0.}; int genSignal[1]={6};' # 0: Bs
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={1.296e+05, 1.954e+04, 4.672e+03, 3.522e+02, 4.040e+01, 0.}; int genSignal[1]={7};' # 1: prompt psi'
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={2.653e+04, 1.979e+04, 6.670e+03, 3.053e+02, 2.768e+01, 0.}; int genSignal[1]={7};' # 2: prompt X->jpsi + rho
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={5.790e+07, 1.734e+07, 5.364e+06, 5.793e+05, 7.805e+04, 0.}; int genSignal[1]={1};' # 3: Bu
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={1.230e+07, 2.868e+06, 8.902e+05, 7.930e+04, 9.631e+03, 0.}; int genSignal[1]={7};' # 4: nonprompt psi'
    'const int nBins=5; float pthatBin[nBins]={5, 10, 15, 30, 50}; float crosssec[nBins+1]={7.622e+06, 1.776e+06, 5.638e+05, 4.973e+04, 6.273e+03, 0.}; int genSignal[1]={7};' # 5: nonprompt X->jpsi + rho
)

tmp=$(date +%y%m%d%H%M%S)
sed '1i'"${crosssec[$chan]}" weighPurePthat.C > weighPurePthat_${tmp}.C

##

g++ addbranch.C $(root-config --cflags --libs) -g -o addbranch_${tmp}.exe || exit 1
g++ weighPurePthat_${tmp}.C $(root-config --cflags --libs) -g -o weighPurePthat_${tmp}.exe || { rm weighPurePthat_${tmp}.C ; rm addbranch_${tmp}.exe ; exit 1 ; }



[[ ${1:-0} -eq 1 ]] && { 
    rsync --progress "$inputfile" "$outputfile"
    set -x
    ./weighPurePthat_${tmp}.exe "$inputfile" "$outputfile" 
    set +x
}
rm weighPurePthat_${tmp}.exe
rm weighPurePthat_${tmp}.C

