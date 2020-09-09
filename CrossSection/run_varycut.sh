#!/bin/bash

g++ fitX_drawvary.C -I../includes/ $(root-config --libs --cflags) -g -o fitX_drawvary.exe || exit 1

[[ ${1:-0} -eq 1 ]] && ./fitX_drawvary.exe rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3_pt15-50_cent090_y0p0-1p6/fitX_fithist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3_pt15-50_cent090_y0p0-1p6

rm fitX_drawvary.exe
