#!/bin/bash

input="../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt15-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt15-50_cent090_y0p0-1p6"
lumiscale=5.88235

##
RUN_FITHIST=${1:-0}

tmp=$(date +%y%m%d%H%M%S)

##
echo -e "\e[32;1mcompiling...\e[0m"
[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "projection.cc
    g++ projection.cc -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o projection_${tmp}.exe || exit 1
}

inputtwo=($input)
inputfile=${inputtwo[0]}
outputdir=${inputtwo[1]}
echo -e "----------------------------------------"
echo -e "==> Input file: \e[4m$inputfile\e[0m"
echo -e "==> Output dir: \e[4m$outputdir\e[0m"
echo -e "----------------------------------------"

[[ $RUN_FITHIST -eq 1 ]] && {
    ./projection_${tmp}.exe "$inputfile" $outputdir $lumiscale
}

[[ -f projection_${tmp}.exe ]] && rm projection_${tmp}.exe

