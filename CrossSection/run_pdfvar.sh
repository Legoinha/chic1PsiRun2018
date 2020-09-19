#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

fitopt=(
    "poly4 trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3 BDTQvalue_PVz15_newL2L3"
    "poly3 trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3 BDTQvalue_PVz15_newL2L3"
    "floatwidth trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3 BDTQvalue_PVz15_newL2L3"
    "3gaus trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3 BDTQvalue_PVz15_newL2L3"
    "range trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3_range BDTQvalue_PVz15_newL2L3"
    "poly2 trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3 BDTQvalue_PVz15_newL2L3"
    "poly1 trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3 BDTQvalue_PVz15_newL2L3"
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
for fff in "${fitopt[@]}"
do
    IFS=' ' read -r -a sff <<< "$fff"
    ff=${sff[0]}
    input=${sff[1]}
    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${sff[2]}
    rootdir=rootfiles/$input$kinematic/

    echo -e "----------------------------------------"
    echo -e "==> Output directory: \e[4m$name/$ff\e[0m"
    echo -e "==> Input directory: \e[4m$rootdir\e[0m"
    echo -e "----------------------------------------"

    [[ $RUN_FITHIST -eq 1 ]] && {
        ./fitX_fithist_${tmp}.exe "$rootdir/fitX_savehist.root" "$name" $ff
    }
done

[[ $RUN_DRAWHIST -eq 1 ]] && ./pdf_drawhist_${tmp}.exe "$name$kinematic"

rm *.exe 2> /dev/null

