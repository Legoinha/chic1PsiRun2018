#!/bin/bash
# #!/usr/local/bin/bash

subdir=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3_pt15-50_cent090_y0p0-1p6
subdir_range=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_PVz15_newL2L3_range_pt15-50_cent090_y0p0-1p6

[[ ${1:-0} == 0 ]] || { g++ drawfitX.cc -I"../includes" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o drawfitX.exe || exit 1 ; }
[[ ${2:-0} == 0 ]] || { g++ readfitX.cc -I"../includes" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o readfitX.exe || exit 1 ; }
[[ ${3:-0} == 0 ]] || { g++ printfitX.cc -I"../includes" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o printfitX.exe || exit 1 ; }

[[ ${1:-0} -eq 1 ]] && { ./drawfitX.exe ../CrossSection/rootfiles/$subdir/fitX_savehist.root $subdir ; }

[[ ${1:-0} -eq 2 ]] && { 
    ./drawfitX.exe ../CrossSection/rootfiles/$subdir/fitX_savehist.root $subdir ; 
    ./drawfitX.exe ../CrossSection/rootfiles/$subdir/fitX_savehist.root $subdir "pol4" ;
    ./drawfitX.exe ../CrossSection/rootfiles/$subdir/fitX_savehist.root $subdir "pol3" ;
    ./drawfitX.exe ../CrossSection/rootfiles/$subdir/fitX_savehist.root $subdir "3gaus" ;
    ./drawfitX.exe ../CrossSection/rootfiles/$subdir/fitX_savehist.root $subdir "floatwidth" ;
    ./drawfitX.exe ../CrossSection/rootfiles/$subdir_range/fitX_savehist.root $subdir "range" ;
}

[[ ${2:-0} -eq 1 ]] && { ./readfitX.exe $subdir "pol4,pol3,3gaus,floatwidth,range" ; }
# [[ ${2:-0} -eq 1 ]] && { ./readfitX.exe $subdir "pol4" ; }

[[ ${3:-0} -eq 1 ]] && { ./printfitX.exe $subdir "pol4,pol3,3gaus,floatwidth,range" ; }

[[ -f printfitX.exe ]] && rm printfitX.exe
[[ -f readfitX.exe ]] && rm readfitX.exe
[[ -f drawfitX.exe ]] && rm drawfitX.exe

