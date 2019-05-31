#!/bin/bash

inputpars=(
    "/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_Pthat5_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat.root hiEvtAnalyzer/HiTree sample 5."
)

g++ addbranch.C $(root-config --libs --cflags) -g -o addbranch.exe || exit 1

for ii in ${inputpars[@]}
do
    pars=($ii)
    cpfile=${pars[0]%%_addSamplePthat.root}.root
    set -x
    cp $cpfile ${pars[0]}

    ./addbranch.exe $ii
done

rm addbranch.exe
set +x
