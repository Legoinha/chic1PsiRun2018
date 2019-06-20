#!/bin/bash

##
ptmin=10 ; ptmax=-1 ; output=
outmvas=(0)

##
# -- signal sample
inputs=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root ; output=rootfiles/TMVA_trainX ;
# -- background sample
inputb=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntmix_Xpt10.root ; bkgstrategy=sideband ;
# inputb=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10.root ; bkgstrategy=samesign ;
output=${output}_$bkgstrategy
# -- data sample
# (((( ioutmva ))))
inputms=(
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntmix_Xpt10.root                 # 0
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_samesign_20190513_HIDoubleMuonPsi_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntmix_Xpt10.root # 1
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root                 # 2
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root                  # 3
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root              # 4
    /export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_X3872ToJpsiRho_nonprompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root               # 5
)
outputmvas=(
    /export/d00/scratch/jwang/BntupleRun2018/output_mva_data         # 0
    /export/d00/scratch/jwang/BntupleRun2018/output_mva_samesign     # 1
    /export/d00/scratch/jwang/BntupleRun2018/output_mva_promptPsi    # 2
    /export/d00/scratch/jwang/BntupleRun2018/output_mva_promptX      # 3
    /export/d00/scratch/jwang/BntupleRun2018/output_mva_nonpromptPsi # 4
    /export/d00/scratch/jwang/BntupleRun2018/output_mva_nonpromptX   # 5
)

##
# -- event filter
cut="pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter"
# -- jpsi
cut=$cut" && Bpt > 10 && TMath::Abs(Bmumumass-3.096916) < 0.05 && TMath::Abs(Bujeta) < 2.4"
# -- HLT
cut=$cut" && HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1"
# -- muon
cut=$cut" && (Bmu1SoftMuID && Bmu2SoftMuID && Bmu1isAcc && Bmu2isAcc && Bmu1isTriggered && Bmu2isTriggered)"
# -- track kinematics
cut=$cut" && Btrk1Pt > 0.9 && Btrk2Pt > 0.9 && TMath::Abs(Btrk1Eta) < 2.4 && TMath::Abs(Btrk2Eta) < 2.4"
# -- track qualirty
cut=$cut" && Btrk1highPurity && Btrk2highPurity && (Btrk1PixelHit+Btrk1StripHit) >= 11 && (Btrk2PixelHit+Btrk2StripHit) >= 11 && TMath::Abs(Btrk1PtErr/Btrk1Pt) < 0.1 && TMath::Abs(Btrk2PtErr/Btrk2Pt) < 0.1 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer)) < 0.18 && (Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer)) < 0.18"
# -- X prefilter
cut=$cut" && TMath::Abs(By) < 2.4 && Bchi2cl > 0.1"
# -- tricky selection
cut=$cut" && BsvpvDisErr>1.e-5 && BsvpvDisErr_2D>1.e-5"
# -- additional selection
cut=$cut" && (Bmass-3.096916-Btktkmass) < 0.2" ; output=${output}_tktk0p2 ;

##
algo="BDT,BDTG,CutsGA,CutsSA,LD"
# algo="BDT,BDTG"

# stages="0,10,1,2,9,8,4,5,6,7,15,16" ; sequence=1 ; # see definition below #
# stages="0,10,17,18,20" ; sequence=0 ; # see definition below #
stages="0,10,1,2,9" ; sequence=0 ; # see definition below #
# stages="0,10,1,2" ; sequence=0 ; # see definition below #
# stages="0,10,1,2,9,8" ; sequence=0 ; # see definition below #
# stages="0,10,1,2,9,8,4,6" ; sequence=0 ; # see definition below #

## ===== do not change the lines below =====
varstrategy=("Single set" "Sequence")
varlist=(
    '#  0  :  ("Bchi2cl"    , "Bchi2cl"                                                                                        , "FMax")  #' 
    '#  1  :  ("dRtrk1"     , "dRtrk1 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk1Phi)),2) + pow(Bujeta-Btrk1Eta,2))", "FMin")  #' 
    '#  2  :  ("dRtrk2"     , "dRtrk2 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk2Phi)),2) + pow(Bujeta-Btrk2Eta,2))", "FMin")  #' 
    '#  3  :  ("Qvalue"     , "Qvalue := (Bmass-3.096916-Btktkmass)"                                                           , "FMin")  #' 
    '#  4  :  ("Balpha"     , "Balpha"                                                                                         , "FMin")  #' 
    '#  5  :  ("costheta"   , "costheta := TMath::Cos(Bdtheta)"                                                                , "FMax")  #' 
    '#  6  :  ("dls3D"      , "dls3D := TMath::Abs(BsvpvDistance/BsvpvDisErr)"                                                 , "FMax")  #' 
    '#  7  :  ("dls2D"      , "dls2D := TMath::Abs(BsvpvDistance_2D/BsvpvDisErr_2D)"                                           , "FMax")  #' 
    '#  8  :  ("Btrk1Pt"    , "Btrk1Pt"                                                                                        , "FMax")  #' 
    '#  9  :  ("Btrk2Pt"    , "Btrk2Pt"                                                                                        , "FMax")  #' 
    '#  10 :  ("trkptimba"  , "trkptimba := TMath::Abs((Btrk1Pt-Btrk2Pt) / (Btrk1Pt+Btrk2Pt))"                                 , "FMax")  #' 
    '#  11 :  ("By"         , "By"                                                                                             , ""    )  #'
    '#  12 :  ("Bmass"      , "Bmass"                                                                                          , ""    )  #'
    '#  13 :  ("Btrk1Eta"   , "Btrk1Eta"                                                                                       , ""    )  #'
    '#  14 :  ("Btrk2Eta"   , "Btrk2Eta"                                                                                       , ""    )  #'
    '#  15 :  ("Btrk1DxySig", "Btrk1DxySig := TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)"                                            , ""    )  #'
    '#  16 :  ("Btrk2DxySig", "Btrk2DxySig := TMath::Abs(Btrk2Dxy1/Btrk2DxyError1)"                                            , ""    )  #'
    '#  17 :  ("dRtrkH"     , "dRtrkH := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-BtrkHPhi)),2) + pow(Bujeta-BtrkHEta,2))", "FMin")  #' 
    '#  18 :  ("dRtrkL"     , "dRtrkL := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-BtrkLPhi)),2) + pow(Bujeta-BtrkLEta,2))", "FMin")  #' 
    '#  19 :  ("BtrkHPt"    , "BtrkHPt"                                                                                        , "FMax")  #' 
    '#  20 :  ("BtrkLPt"    , "BtrkLPt"                                                                                        , "FMax")  #' 
    '#  21 :  ("BtrkHEta"   , "BtrkHEta"                                                                                       , ""    )  #'
    '#  22 :  ("BtrkLEta"   , "BtrkLEta"                                                                                       , ""    )  #'
    '#  23 :  ("BtrkHDxySig", "BtrkHDxySig := TMath::Abs(BtrkHDxy1/BtrkHDxyError1)"                                            , ""    )  #'
    '#  24 :  ("BtrkLDxySig", "BtrkLDxySig := TMath::Abs(BtrkLDxy1/BtrkLDxyError1)"                                            , ""    )  #'
)


cuts="${cut} && Bgen==23333 && BgencollisionId==0"
[[ $bkgstrategy == "sideband" ]] && cutb="${cut} && TMath::Abs(Bmass-3.8719) > 0.07 && TMath::Abs(Bmass-3.8719) < 0.128" # sideband
[[ $bkgstrategy == "samesign" ]] && cutb="${cut} && TMath::Abs(Bmass-3.8719) < 0.04" # samesign

## ===== do not do not do not change the lines below =====
IFS=','; allstages=($stages); unset IFS;
echo -e '

######################################################################################################################################
#                                                 Variables \e[1;32m(To be used)\e[0m                                                             #
######################################################################################################################################
#                                                                                                                                    #'
vv=0
while(( $vv < ${#varlist[@]})) ; do
    for ss in ${allstages[@]} ; do [[ $ss == $vv ]] && { echo -en "\e[1;32m" ; break ; } ; done ;
    echo -e "${varlist[vv]}\e[0m"
    vv=$((vv+1))
done
echo '#                                                                                                                                    #
######################################################################################################################################
'

##
echo -e "
######################################################################################################################################
#                                                   Training Configurations                                                          #
######################################################################################################################################

>>>>>> Variables training strategy
  >>>> \e[32m[${varstrategy[sequence]}]\e[0m

>>>>>> Background strategy
  >>>> \e[32m[$bkgstrategy]\e[0m

>>>>>> Selections
  >>>> Prefilters
    >> \e[32m\"$cut\"\e[0m
  >>>> Signal cut
    >> \e[32m\"${cuts##$cut}\"\e[0m
  >>>> Background cut
    >> \e[32m\"${cutb##$cut}\"\e[0m
"

##
mkdir -p $output ; rm -r $output ;
tmp=$(date +%y%m%d%H%M%S)

##
[[ $# -eq 0 ]] && echo "usage: ./run_TMVAClassification.sh [train] [draw curves] [create BDT tree]"

g++ TMVAClassification.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o TMVAClassification_${tmp}.exe || { exit 1 ; }
g++ guivariables.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guivariables_${tmp}.exe             || { rm *_${tmp}.exe ; exit 1 ; }
g++ guiefficiencies.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guiefficiencies_${tmp}.exe       || { rm *_${tmp}.exe ; exit 1 ; }
g++ guieffvar.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guieffvar_${tmp}.exe                   || { rm *_${tmp}.exe ; exit 1 ; }
g++ mvaprod.C $(root-config --libs --cflags) -lTMVA -lXMLIO -lstdc++fs -g -o mvaprod_${tmp}.exe              || { rm *_${tmp}.exe ; exit 1 ; }

[[ ${1:-0} -eq 1 ]] && {
    conf=
    echo -e "\e[2m==> Do you really want to run\e[0m \e[1mTMVAClassification.C\e[0m \e[2m(it might be very slow)?\e[0m [y/n]"
    read conf
    while [[ $conf != 'y' && $conf != 'n' ]] ; do { echo "warning: input [y/n]" ; read conf ; } ; done ;
    [[ $conf == 'n' ]] && { rm *_${tmp}.exe ; exit ; }
}

stage=$stages
while [[ $stage == *,* ]]
do
# train
    [[ ${1:-0} -eq 1 ]] && { ./TMVAClassification_${tmp}.exe $inputs $inputb "$cuts" "$cutb" $output $ptmin $ptmax "$algo" "$stage"; } &
    [[ $sequence -eq 0 ]] && break;
    while [[ $stage != *, ]] ; do stage=${stage%%[0-9]} ; done ;
    stage=${stage%%,}
done
wait

# draw curves
[[ ${2:-0} -eq 1 ]] && { 
    ./guivariables_${tmp}.exe $output $ptmin $ptmax "$algo" "$stages"
    ./guiefficiencies_${tmp}.exe $output $ptmin $ptmax "$algo" "$stages"
}

# draw curve vs. var
[[ ${2:-0} -eq 1 && $sequence -eq 1 ]] && ./guieffvar_${tmp}.exe $output $ptmin $ptmax "$algo" "$stages"

# produce mva values
for ioutmva in ${outmvas[@]}
{
    inputm=${inputms[ioutmva]} ; outputmva=${outputmvas[ioutmva]} ;
    [[ ${3:-0} -eq 1 ]] && ./mvaprod_${tmp}.exe $inputm "Bfinder/ntmix" $output $outputmva $ptmin $ptmax "$algo" "${stages}"
}
##
rm *_${tmp}.exe


# inputs=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190520_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_1033p1_pt6tkpt0p7dls0_v3_addSamplePthat_pthatweight.root ; output=rootfiles/TMVA_trainPsi2S ;
