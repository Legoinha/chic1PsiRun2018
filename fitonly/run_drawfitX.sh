#!/usr/local/bin/bash

subdir=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt15-50_cent090_y0p0-1p6

g++ drawfitX.cc -I"../includes" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o drawfitX.exe || exit 1

[[ ${1:-0} -eq 1 ]] && { ./drawfitX.exe rootfiles/$subdir/fitX_savehist.root $subdir ; }

rm drawfitX.exe

