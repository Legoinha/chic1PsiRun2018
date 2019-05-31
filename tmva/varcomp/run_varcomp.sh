#!/bin/bash

ptmin=10 ; ptmax=-1 ;
inputs="/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_X3872ToJpsiRho_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root"
inputb="/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntmix_Bpt${ptmin}.root"
output="rootfiles/root_PsiX_${ptmin}_${ptmax}"

##
nt="Bfinder/ntmix"
leg="Background Signal(#Psi') Signal(X)"

##
mkdir -p rootfiles plots
tmp=$(date +%y%m%d%H%M%S)

echo 'Usage: ./run_varcomp.sh [fill hist] [draw hist]'

g++ varcomp.C $(root-config --libs --cflags) -g -o varcomp_${tmp}.exe || exit 1
g++ drawvarcomp.C $(root-config --libs --cflags) -g -o drawvarcomp_${tmp}.exe || exit 1

[[ ${1:-0} -eq 1 ]] && { ./varcomp_${tmp}.exe $nt $output $ptmin $ptmax $inputb $inputs ; }
[[ ${2:-0} -eq 1 ]] && {
    [[ ! -f ${output}.root ]] && { rm varcomp_${tmp}.exe ; rm *_${tmp}.exe ; exit 1 ; }
    [[ -f ../../merge/showdir ]] && { ../../merge/showdir ${output}.root ; }
    ./drawvarcomp_${tmp}.exe $output $leg
}

rm drawvarcomp_${tmp}.exe
rm varcomp_${tmp}.exe

