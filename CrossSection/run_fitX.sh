#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

counts=(1)
optcuts=(
    "&& BDTF > 0.3 && (Bmass-3.096916-Btktkmass) < 0.12"
    "&& BDT > 0.07 && (Bmass-3.096916-Btktkmass) < 0.12"
    "&& BDTD > 0.12 && (Bmass-3.096916-Btktkmass) < 0.12"
    "&& BDTG > 0.70 && (Bmass-3.096916-Btktkmass) < 0.12"
)
optcutntuples=(
    "ntp->BDTF[j] > 0.3 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.12"
    "ntp->BDT[j] > 0.07 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.12"
    "ntp->BDTD[j] > 0.12 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.12"
    "ntp->BDTG[j] > 0.70 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.12"
)
tags=(
    "BDTFQvalue"
    "BDTQvalue"
    "BDTDQvalue"
    "BDTGQvalue"
)
input=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190806_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_skim_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
input_ss=/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root

##

RUN_SAVEHIST=${1:-0}
RUN_CONVERTTNP=${2:-0}
RUN_FITHIST=${3:-0}
RUN_DRAWHIST=${4:-0}

##
g++ getfname.cc $(root-config --libs --cflags) -g -o getfname.exe || exit 1
kinematic=$(./getfname.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname.exe
echo -e "\e[32;1mcompiling...\e[0m"
cp ../tnp/tnpcc.h tnpcc_tmp.h
ptbins='float ptbins[] = {'$ptmin', '$ptmax'};'
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc_tmp.h
[[ $RUN_SAVEHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitX_flatten.C
    g++ fitX_flatten.C $(root-config --libs --cflags) -g -o fitX_flatten.exe || exit 1 
    echo " -- "fitX_savehist.C
    g++ fitX_savehist.C $(root-config --libs --cflags) -lRooFit -lRooFitCore -g -o fitX_savehist.exe || exit 1 
}
[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "draw_tnp.cc
    g++ draw_tnp.cc $(root-config --libs --cflags) -g -o draw_tnp.exe || exit 1 
    echo " -- "fitX_fithist.C
    g++ fitX_fithist.C $(root-config --libs --cflags) -lRooFit -lRooFitCore -g -o fitX_fithist.exe || exit 1
}
[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitX_drawhist.C
    g++ fitX_drawhist.C $(root-config --libs --cflags) -g -o fitX_drawhist.exe || exit 1
}

##
[[ $RUN_SAVEHIST -eq 1 ]] && {
    flatinput=${input%%.root}_flatten.root
    [[ -f $flatinput ]] || ./fitX_flatten.exe $input $flatinput
    flatinputmc_a_prompt=${inputmc_a_prompt%%.root}_flatten.root
    [[ -f $flatinputmc_a_prompt ]] || ./fitX_flatten.exe $inputmc_a_prompt $flatinputmc_a_prompt
    flatinputmc_b_prompt=${inputmc_b_prompt%%.root}_flatten.root
    [[ -f $flatinputmc_b_prompt ]] || ./fitX_flatten.exe $inputmc_b_prompt $flatinputmc_b_prompt
}

##
echo
for count in ${counts[@]}
do
    cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
    cut="$cut && mvapref ${optcuts[count]}"
    cutgen="1"
    cp ../tnp/tnp_converter.cc tnp_converter_tmp.cc
    sed -i "s/__CUTINPUT__/${optcutntuples[$count]}/g" tnp_converter_tmp.cc
    [[ $RUN_CONVERTTNP -eq 1 || $# -eq 0 ]] && {
        g++ tnp_converter_tmp.cc $(root-config --libs --cflags) -g -o tnp_converter.exe || continue
    }

    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${tags[count]}
    rootdir=rootfiles/$name$kinematic/
    echo -e "----------------------------------------"
    echo -e "==> File directory: \e[4m$rootdir\e[0m"
    echo -e "----------------------------------------"

    [[ $RUN_SAVEHIST -eq 1 ]] && {
        ./fitX_savehist.exe $input $inputmc_a_prompt $inputmc_b_prompt $inputmc_a_nonprompt $inputmc_b_nonprompt "$cut" "$cutgen" "$name" $ptmin $ptmax $centmin $centmax $ymin $ymax
    }
    [[ $RUN_CONVERTTNP -eq 1 ]] && {
        ./tnp_converter.exe $inputmc_a_prompt $name "_a" $ptmin $ptmax $centmin $centmax $ymin $ymax
        ./tnp_converter.exe $inputmc_b_prompt $name "_b" $ptmin $ptmax $centmin $centmax $ymin $ymax
    }
    [[ $RUN_FITHIST -eq 1 ]] && {
        ./draw_tnp.exe "$rootdir/tnp_a.root" $name "_a"
        ./draw_tnp.exe "$rootdir/tnp_b.root" $name "_b"
        ./fitX_fithist.exe "$rootdir/fitX_savehist.root" $name
    }
    [[ $RUN_DRAWHIST -eq 1 ]] && ./fitX_drawhist.exe "$rootdir/fitX_fithist.root" $name

    # syst
    # echo fitX_pdfvar.C
    # g++ fitX_pdfvar.C $(root-config --libs --cflags) -lRooFit -lRooFitCore -g -o fitX_pdfvar.exe || exit 1
    # [[ ${5:-0} -eq 1 ]] && ./fitX_pdfvar.exe "$outputdir/fitX_savehist" $name

    [[ -f tnp_converter.exe ]] && rm tnp_converter.exe
    rm tnp_converter_tmp.cc
done

[[ -f fitX_drawhist.exe ]] && rm fitX_drawhist.exe
[[ -f fitX_fithist.exe ]] && rm fitX_fithist.exe
[[ -f fitX_savehist.exe ]] && rm fitX_savehist.exe
[[ -f fitX_flatten.exe ]] && rm fitX_flatten.exe
[[ -f draw_tnp.exe ]] && rm draw_tnp.exe
rm tnpcc_tmp.h

# rm fitX_pdfvar.exe
