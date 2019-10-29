#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

counts=({0..9})
types=(
    "BDT"      # 0
    "Qvalue"   # 1
    "pt"       # 2
    "absy"     # 3
    "chi2cl"   # 4
    "trk1pt"   # 5
    "trk2pt"   # 6
    "dRtrk1"   # 7
    "dRtrk2"   # 8
    "trkptimb" # 9
    "BDTnoQ"   # 10
    "lxy"      # 11
)
precuts=(
    "&& (Bmass-3.096916-Btktkmass)<0.15"             # 0
    "&& BDT>0.07"                                    # 1
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 2
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 3
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 4 
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 5
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 6
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 7
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 8
    "&& (Bmass-3.096916-Btktkmass)<0.12 && BDT>0.07" # 9
    "&& BDT>-1"                                      # 10
    "&& BDT>0.04"                                    # 11
)

input=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190806_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_skim_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root

RUN_SAVEHIST=${1:-0}
RUN_FITHIST=${2:-0}
RUN_DRAWHIST=${3:-0}

tmp=$(date +%y%m%d%H%M%S)

##
g++ getfname.cc -I"../includes/" $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe
echo -e "\e[32;1mcompiling...\e[0m"

[[ $RUN_SAVEHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "datamc.cc
    g++ datamc.cc -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o datamc_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitdatamc.cc
    g++ fitdatamc.cc -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o fitdatamc_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "drawdatamc.cc
    g++ drawdatamc.cc -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o drawdatamc_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

echo
for count in ${counts[@]}
do
    cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
    cut="$cut && mvapref ${precuts[count]}"
    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9
    rootdir=rootfiles/$name$kinematic/${types[count]}/
    echo -e "----------------------------------------"
    echo -e "==> File directory: \e[4m$rootdir\e[0m"
    echo -e "----------------------------------------"
    
    [[ $RUN_SAVEHIST -eq 1 ]] && {
        ./datamc_${tmp}.exe $input $inputmc_a_prompt $inputmc_b_prompt "$cut" "${types[count]}" $name $ptmin $ptmax $centmin $centmax $ymin $ymax
    }
    [[ $RUN_FITHIST -eq 1 ]] && {
        ./fitdatamc_${tmp}.exe "${rootdir}/datamc_savehist.root" $name "${types[count]}"
    }
    [[ $RUN_DRAWHIST -eq 1 ]] && {
        ./drawdatamc_${tmp}.exe "${rootdir}/datamc_fithist.root" $name "${types[count]}"
    }
done

rm drawdatamc_${tmp}.exe 2> /dev/null
rm fitdatamc_${tmp}.exe 2> /dev/null
rm datamc_${tmp}.exe 2> /dev/null

