#!/bin/bash

inputfiles=(
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXPiPi_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight.root
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190712_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight.root
)

g++ fillpthatweight.C $(root-config --libs --cflags) -g -o fillpthatweight.exe || exit
g++ drawpthatweight.C $(root-config --libs --cflags) -g -o drawpthatweight.exe || exit

[[ ${1:-0} -eq 1 ]] && {
    for ii in ${inputfiles[@]}
    do
        echo -e "\e[33;1m${ii}\e[0m"
        ./fillpthatweight.exe $ii
    done
}

[[ ${2:-0} -eq 1 ]] && {
    for ii in `echo rootfiles/*.root`
    do
        echo -e "\e[33;1m${ii}\e[0m"
        tag=${ii%%"_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight.root"}
        tag=${tag##"rootfiles/pthatweight_crab_Bfinder_20190712_Hydjet_Pythia8_"}
        echo $tag
        legend=
        if [[ $tag == NonPrompt* ]] ; then 
            legend="(b#rightarrow)"
        elif [[ $tag == Prompt* ]] ; then
            legend="Prompt"
        fi

        if [[ $tag == *Psi2S ]] ; then
            legend="$legend #psi(2S) #rightarrow J/#psi#pi#pi"
        elif [[ $tag == *XRho ]] ; then 
            legend="$legend X(3872) #rightarrow J/#psi#rho"
        elif [[ $tag == *XPiPi ]] ; then
            legend="$legend X(3872) #rightarrow J/#psi#pi#pi"
        fi
        echo $legend
        ./drawpthatweight.exe $ii "$legend"
        echo
    done
}

rm drawpthatweight.exe
rm fillpthatweight.exe