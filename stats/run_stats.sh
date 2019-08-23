#!/bin/bash
counts=(6)
inputs=(
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTDQvalue_pt15-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTDQvalue_pt15-50_cent090_y0p0-1p6" # 0
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTDQvalue_pt20-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTDQvalue_pt20-50_cent090_y0p0-1p6"
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTFQvalue_pt15-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTFQvalue_pt15-50_cent090_y0p0-1p6" # 2
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTFQvalue_pt20-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTFQvalue_pt20-50_cent090_y0p0-1p6"
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTGQvalue_pt15-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTGQvalue_pt15-50_cent090_y0p0-1p6" # 4
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTGQvalue_pt20-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTGQvalue_pt20-50_cent090_y0p0-1p6"
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt15-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt15-50_cent090_y0p0-1p6" # 6
    "../CrossSection/rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt20-50_cent090_y0p0-1p6/fitX_savehist.root trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue_pt20-50_cent090_y0p0-1p6"
)

##
RUN_FITHIST=${1:-0}
RUN_DRAWHIST=${2:-0}

tmp=$(date +%y%m%d%H%M%S)

##
echo -e "\e[32;1mcompiling...\e[0m"
[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "toymc.cc
    g++ toymc.cc $(root-config --libs --cflags) -lRooFit -lRooFitCore -g -o toymc_${tmp}.exe || exit 1
}
[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "drawtoymc.cc
    g++ drawtoymc.cc $(root-config --libs --cflags) -lRooFit -lRooFitCore -g -o drawtoymc_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
}

echo
for count in ${counts[@]}
do
    input=(${inputs[count]})
    inputfile=${input[0]}
    outputdir=${input[1]}
    echo -e "----------------------------------------"
    echo -e "==> Input file: \e[4m$inputfile\e[0m"
    echo -e "==> Output dir: \e[4m$outputdir\e[0m"
    echo -e "----------------------------------------"

    [[ $RUN_FITHIST -eq 1 ]] && {
        ./toymc_${tmp}.exe "$inputfile" $outputdir
    }
    [[ $RUN_DRAWHIST -eq 1 ]] && ./drawtoymc_${tmp}.exe rootfiles/$outputdir/toymc.root $outputdir

done

[[ -f drawtoymc_${tmp}.exe ]] && rm drawtoymc_${tmp}.exe
[[ -f toymc_${tmp}.exe ]] && rm toymc_${tmp}.exe

