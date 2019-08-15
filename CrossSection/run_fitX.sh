#!/bin/bash

ptmin=20
ptmax=50
centmin=0
centmax=90
ymin=0
ymax=1.6

counts=(0)
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
input=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/ntmix_20190806_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_HF_and_MuonJSON_skimBpt10_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
input_ss=/export/d00/scratch/jwang/BntupleRun2018/mva_output_old/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
inputmc_a_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/crab_Bfinder_20190712_Hydjet_Pythia8_PromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_prompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/crab_Bfinder_20190730_Hydjet_Pythia8_PromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_a_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptPsi2S_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root
inputmc_b_nonprompt=/raid5/data/wangj/BntupleRun2018/mva_output_20190808ptdep/crab_Bfinder_20190712_Hydjet_Pythia8_NonPromptXRho_1033p1_official_pt6tkpt0p9dls0_skimhltBsize_pthatweight_trainX_20190808ptdep_sideband_tktk0p2_BDT_BDTD_BDTG_BDTF_LD_15p0_50p0_0-10-1-2-9_1bin.root

cp ../tnp/tnpcc.h tnpcc_tmp.h
ptbins='float ptbins[] = {'$ptmin', '$ptmax'};'
sed -i "s/__PTBIN_INPUT__/$ptbins/g" tnpcc_tmp.h
g++ draw_tnp.cc $(root-config --libs --cflags) -g -o draw_tnp.exe || exit 1
g++ fitX_savehist.C $(root-config --libs --cflags) -g -o fitX_savehist.exe || exit 1
g++ fitX_fithist.C $(root-config --libs --cflags) -g -o fitX_fithist.exe || exit 1
g++ fitX_drawhist.C $(root-config --libs --cflags) -g -o fitX_drawhist.exe || exit 1

# syst
g++ fitX_pdfvar.C $(root-config --libs --cflags) -g -o fitX_pdfvar.exe || exit 1

for count in ${counts[@]}
do
    cut="HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
    cut="$cut && mvapref ${optcuts[count]}"
    cutgen="1"
    cp ../tnp/tnp_converter.cc tnp_converter_tmp.cc
    sed -i "s/__CUTINPUT__/${optcutntuples[$count]}/g" tnp_converter_tmp.cc
    g++ tnp_converter_tmp.cc $(root-config --libs --cflags) -g -o tnp_converter.exe || continue

    name=trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_${tags[count]}

    [[ ${1:-0} -eq 1 ]] && ./fitX_savehist.exe $input $inputmc_a_prompt $inputmc_b_prompt $inputmc_a_nonprompt $inputmc_b_nonprompt "$cut" "$cutgen" "$name" $ptmin $ptmax $centmin $centmax $ymin $ymax
    [[ ${2:-0} -eq 1 ]] && {
        ./tnp_converter.exe $inputmc_a_prompt $name "_a" $ptmin $ptmax $centmin $centmax $ymin $ymax ;
        ./tnp_converter.exe $inputmc_b_prompt $name "_b" $ptmin $ptmax $centmin $centmax $ymin $ymax ;
        ./draw_tnp.exe $name "_a" ;
        ./draw_tnp.exe $name "_b" ;
    }
    [[ ${3:-0} -eq 1 ]] && ./fitX_fithist.exe "$name"
    [[ ${4:-0} -eq 1 ]] && ./fitX_drawhist.exe "$name"

    # [[ ${5:-0} -eq 1 ]] && ./fitX_pdfvar.exe "$outputdir/fitX_savehist" "$name"

    rm tnp_converter.exe
    rm tnp_converter_tmp.cc
done

rm fitX_drawhist.exe
rm fitX_fithist.exe
rm fitX_savehist.exe
rm draw_tnp.exe
rm tnpcc_tmp.h

rm fitX_pdfvar.exe

# input=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
# input_ss=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
# inputmc_a_prompt=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
# inputmc_b_prompt=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
# inputmc_a_nonprompt=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
# inputmc_b_nonprompt=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight_trainX_sideband_tktk0p2_BDT_BDTG_CutsGA_CutsSA_LD_10p0_inf_0-10-1-2-9_oldLH.root
