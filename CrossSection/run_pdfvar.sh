#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

tags="BDTQvalue_PVz15_newL2L3"

fitopt=(
    "poly4"
    "poly3"
    "floatwidth"
    "3gaus"
)

RUN_FITHIST=${1:-0}
RUN_DRAWHIST=${2:-0}

tmp=$(date +%y%m%d%H%M%S)

##
g++ getfname.cc -I"../includes/" $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe

[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitX_fithist.C
    g++ fitX_fithist.C -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o fitX_fithist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}
[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "pdf_drawhist.C
    g++ pdf_drawhist.C -I"../includes/" $(root-config --libs --cflags) -g -o pdf_drawhist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

##
for ff in ${fitopt[@]}
do

    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${tags}
    echo -e "----------------------------------------"
    echo -e "==> File directory: \e[4m$name/$ff\e[0m"
    echo -e "----------------------------------------"

    rootdir=rootfiles/$name$kinematic/

    [[ $RUN_FITHIST -eq 1 ]] && {
        ./fitX_fithist_${tmp}.exe "$rootdir/fitX_savehist.root" "$name" $ff
    }
done

[[ $RUN_DRAWHIST -eq 1 ]] && ./pdf_drawhist_${tmp}.exe "$name$kinematic"

rm *.exe 2> /dev/null

