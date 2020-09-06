#!/bin/bash

ptmin=15
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

counts=(0)
optcuts=(
    "&& BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13" # 0: nominal
    "&& BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13" # 1
    "&& BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13" # 2
    "&& BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13" # 3
    "&& BDT > 0.06 && (Bmass-3.096916-Btktkmass) < 0.13" # 4
    "&& BDTF > 0.3 && (Bmass-3.096916-Btktkmass) < 0.13" # 5
    "&& BDTD > 0.11 && (Bmass-3.096916-Btktkmass) < 0.13" # 6
    "&& BDTG > 0.68 && (Bmass-3.096916-Btktkmass) < 0.13" # 7
    "&& BDT > 0.04"
    ""
)
optcutntuples=(
    "ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 0: nominal
    "ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 1
    "ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 2
    "ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 3
    "ntp->BDT[j] > 0.06 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 4
    "ntp->BDTF[j] > 0.3 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 5
    "ntp->BDTD[j] > 0.11 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 6
    "ntp->BDTG[j] > 0.68 \&\& (ntp->Bmass[j]-3.096916-ntp->Btktkmass[j]) < 0.13" # 7
    "ntp->BDT[j] > 0.04"
    "true"
)
# --> tags including "L2L3" to use correct tnp leg matching <-- see: tnp_converter.cc
# --> tags including "PVz15" to apply |PVz| < 15
# --> tags including "newL2L3" to apply correct tnp header
tags=(
    "BDTQvalue_PVz15_newL2L3"  # 0: nominal
    "BDTQvalue_PVz15_L2L3"     # 1
    "BDTQvalue_PVz15"          # 2
    "BDTQvalue"                # 3
    "BDTQvalue_wmissingfile"   # 4
    "BDTFQvalue_PVz15_newL2L3" # 5
    "BDTDQvalue_PVz15_newL2L3" # 6
    "BDTGQvalue_PVz15_newL2L3" # 7
    "BDTLoose2S"
    "noBDT"
)

input=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190806_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_skim_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
input_ss=/export/d00/scratch/jwang/BntupleRun2018/mva_output/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_mutrg_20190808_20200830rmevt_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_pthatweight_trainX_20190808ptdep.root
# inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_mutrg_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_pthatweight_trainX_20190808ptdep.root
# inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190808_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_pt6tkpt0p9dls0_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root

##

RUN_SAVEHIST=${1:-0}
RUN_CONVERTTNP=${2:-0}
RUN_FITHIST=${3:-0}
RUN_DRAWHIST=${4:-0}

tmp=$(date +%y%m%d%H%M%S)

##
g++ getfname.cc -I"../includes/" $(root-config --libs --cflags) -g -o getfname_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
kinematic=$(./getfname_${tmp}.exe $ptmin $ptmax $centmin $centmax $ymin $ymax)
rm getfname_${tmp}.exe
echo -e "\e[32;1mcompiling...\e[0m"
cp ../tnp/tnpcc.h tnpcc_tmp.h
ptbins='float ptbins[] = {'$ptmin', '$ptmax'};'
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc_tmp.h
[[ $RUN_SAVEHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitX_flatten.C
    g++ fitX_flatten.C -I"../includes/" $(root-config --libs --cflags) -g -o fitX_flatten_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; } 
    echo " -- "fitX_savehist.C
    g++ fitX_savehist.C -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o fitX_savehist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; } 
}
[[ $RUN_FITHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "draw_tnp.cc
    g++ draw_tnp.cc -I"../includes/" $(root-config --libs --cflags) -g -o draw_tnp_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; } 
    echo " -- "fitX_fithist.C
    g++ fitX_fithist.C -I"../includes/" $(root-config --libs --cflags) -lRooFit -lRooFitCore -lRooStats -g -o fitX_fithist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}
[[ $RUN_DRAWHIST -eq 1 || $# -eq 0 ]] && {
    echo " -- "fitX_drawhist.C
    g++ fitX_drawhist.C -I"../includes/" $(root-config --libs --cflags) -g -o fitX_drawhist_${tmp}.exe || { rm *_${tmp}.exe 2> /dev/null ; exit 1 ; }
}

##
[[ $RUN_SAVEHIST -eq 1 ]] && {
    flatinput=${input%%.root}_flatten.root
    [[ -f $flatinput ]] || ./fitX_flatten_${tmp}.exe $input $flatinput
    flatinputmc_a_prompt=${inputmc_a_prompt%%.root}_flatten.root
    [[ -f $flatinputmc_a_prompt ]] || ./fitX_flatten_${tmp}.exe $inputmc_a_prompt $flatinputmc_a_prompt
    flatinputmc_b_prompt=${inputmc_b_prompt%%.root}_flatten.root
    [[ -f $flatinputmc_b_prompt ]] || ./fitX_flatten_${tmp}.exe $inputmc_b_prompt $flatinputmc_b_prompt
}

##
echo
for count in ${counts[@]}
do
    precut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && mvapref"
    [[ "${tags[count]}" == *"PVz15"* ]] && { precut="$precut && fabs(PVz) < 15" ; }
    cut="$precut ${optcuts[count]}"
    cutgen="1"
    cutacc="(fabs(Gmu1eta) < 2.4) && ((fabs(Gmu1eta) < 1.2 && Gmu1pt >= 3.5) || (fabs(Gmu1eta) >= 1.2 && fabs(Gmu1eta) < 2.1 && Gmu1pt >= 5.47-1.89*fabs(Gmu1eta)) || (fabs(Gmu1eta) >= 2.1 && Gmu1pt >= 1.5)) && (fabs(Gmu2eta) < 2.4) && ((fabs(Gmu2eta) < 1.2 && Gmu2pt >= 3.5) || (fabs(Gmu2eta) >= 1.2 && fabs(Gmu2eta) < 2.1 && Gmu2pt >= 5.47-1.89*fabs(Gmu2eta)) || (fabs(Gmu2eta) >= 2.1 && Gmu2pt >= 1.5)) && Gtk1pt > 0.9 && Gtk2pt > 0.9 && fabs(Gtk1eta) < 2.4 && fabs(Gtk2eta) < 2.4"
    cp ../tnp/tnp_converter.cc tnp_converter_tmp_${tmp}.cc
    sed -i "s/__CUTINPUT__/${optcutntuples[$count]}/g" tnp_converter_tmp_${tmp}.cc
    [[ "${tags[count]}" != *"newL2L3"* ]] && { sed -i "s/tnp_weight_lowptPbPb/obsolete\/tnp_weight_lowptPbPb/g" tnp_converter_tmp_${tmp}.cc ; }
    [[ $RUN_CONVERTTNP -eq 1 || $# -eq 0 ]] && {
        g++ tnp_converter_tmp_${tmp}.cc -I"../includes/" $(root-config --libs --cflags) -g -o tnp_converter_${tmp}.exe || { rm tnp_converter_tmp_${tmp}.cc ; continue ; }
    }
    rm tnp_converter_tmp_${tmp}.cc

    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${tags[count]}
    echo -e "----------------------------------------"
    echo -e "==> File directory: \e[4m$name\e[0m"
    echo -e "----------------------------------------"

    [[ $RUN_SAVEHIST -eq 1 ]] && {
        ./fitX_savehist_${tmp}.exe $input $inputmc_a_prompt $inputmc_b_prompt $inputmc_a_nonprompt $inputmc_b_nonprompt "$cut" "$cutgen" "$cutacc" "$precut" "$name" $ptmin $ptmax $centmin $centmax $ymin $ymax
    }

    rootdir=rootfiles/$name$kinematic/

    [[ $RUN_CONVERTTNP -eq 1 ]] && {
        ./tnp_converter_${tmp}.exe $inputmc_a_prompt $name "_a" "noweight" $ptmin $ptmax $centmin $centmax $ymin $ymax
        ./tnp_converter_${tmp}.exe $inputmc_b_prompt $name "_b" "noweight" $ptmin $ptmax $centmin $centmax $ymin $ymax
        rm tnp_converter_${tmp}.exe
    }

    [[ $RUN_FITHIST -eq 1 ]] && {
        ./draw_tnp_${tmp}.exe "$rootdir/tnp_a.root" $name "_a"
        ./draw_tnp_${tmp}.exe "$rootdir/tnp_b.root" $name "_b"
        ./fitX_fithist_${tmp}.exe "$rootdir/fitX_savehist.root" $name
    }
    [[ $RUN_DRAWHIST -eq 1 ]] && ./fitX_drawhist_${tmp}.exe "$rootdir/fitX_fithist.root" "$name$kinematic"

done

rm *_${tmp}.exe 2>/dev/null
rm tnpcc_tmp.h 2>/dev/null
