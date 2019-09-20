#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

RUN_DRAWHIST=${1:-0}

tmp=$(date +%y%m%d%H%M%S)

##
g++ getfname.cc -I"../includes/" $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe

[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "algo_drawhist.C
    g++ algo_drawhist.C -I"../includes/" $(root-config --libs --cflags) -g -o algo_drawhist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

##
name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9__ALGO_

echo -e "----------------------------------------"
echo -e "==> File directory: \e[4m$name\e[0m"
echo -e "----------------------------------------"

# rootdir=rootfiles/$name$kinematic/
[[ $RUN_DRAWHIST -eq 1 ]] && ./algo_drawhist_${tmp}.exe "$name$kinematic"

rm *.exe 2> /dev/null

