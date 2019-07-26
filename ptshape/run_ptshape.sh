#!/bin/bash

inputdirname=HEPData-ins1495026-v1-csv

idxi=(0 1 2 3)
input=(
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root promptPsi"
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root promptX"
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root nonpromptPsi"
    "/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9.root nonpromptX"
)

g++ fillptshape.cc $(root-config --libs --cflags) -g -o fillptshape.exe || exit 1
g++ fitptshape.cc $(root-config --libs --cflags) -g -o fitptshape.exe || exit 1

for ii in "${idxi[@]}"
do
    iis=(${input[ii]})
    outputname="rootfiles/$inputdirname/${iis[1]}.root"
    [[ ${1:-0} -eq 1 ]] && {
        set -x
        ./fillptshape.exe ${input[ii]} $outputname $inputdirname &
        set +x
    }
    wait
    [[ ${2:-0} -eq 1 ]] && {
        set -x
        ./fitptshape.exe $outputname ${iis[1]} $inputdirname
        set +x
    }
done

rm fitptshape.exe
rm fillptshape.exe

