#!/bin/bash

##
ptmin=10 ; ptmax=-1 ;
# inputs=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_Psi2SX3872_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root ; output=rootfiles/TMVA_XPsi2S_tktkmass2 ;
# inputs=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root ; output=rootfiles/TMVA_Psi2S_tktkmass2 ;
inputs=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_X3872ToJpsiRho_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root ; output=rootfiles/TMVA_X_tktk0p2 ;
inputb=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntmix_Bpt${ptmin}.root

inputm=$inputb ; outputmva=/export/d00/scratch/jwang/BntupleRun2018/output_mva ;
# inputm=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root ; outputmva=/export/d00/scratch/jwang/BntupleRun2018/output_mva_mcPsi ;
# inputm=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_X3872ToJpsiRho_prompt_20181231_pt5tkpt0p7dls0_pthatweight.root ; outputmva=/export/d00/scratch/jwang/BntupleRun2018/output_mva_mcX ;

cut="Bmu1TMOneStationTight && Bmu1InPixelLayer > 1 && (Bmu1InPixelLayer+Bmu1InStripLayer) > 6 && Bmu1dxyPV < 0.3 && Bmu1dzPV < 20 && Bmu1isGlobalMuon && TMath::Abs(Bmu1eta)<2.0 && Bmu1pt > 1.5 && Bmu2TMOneStationTight && Bmu2InPixelLayer > 1 && (Bmu2InPixelLayer+Bmu2InStripLayer) > 6 && Bmu2dxyPV < 0.3 && Bmu2dzPV < 20 && Bmu2isGlobalMuon && TMath::Abs(Bmu2eta)<2.0 && Bmu2pt > 1.5 && TMath::Abs(Bmumumass-3.096916) < 0.05 && TMath::Abs(Bujeta) < 2.0 && Btrk1highPurity &&  TMath::Abs(Btrk1Eta) < 2 && Btrk1Pt > 0.9 && (Btrk1PixelHit+Btrk1StripHit) > 11 && TMath::Abs(Btrk1PtErr/Btrk1Pt) < 0.1 && Btrk2highPurity &&  TMath::Abs(Btrk2Eta) < 2 && Btrk2Pt > 0.9 && (Btrk2PixelHit+Btrk2StripHit) > 11 && TMath::Abs(Btrk2PtErr/Btrk2Pt) < 0.1 && TMath::Abs(By) < 2.0 && Bchi2cl > 0.1 && BsvpvDisErr>1.e-5 && BsvpvDisErr_2D>1.e-5&&(Bmass-3.096916-Btktkmass)<0.2"
algo="BDT,BDTG,CutsGA,CutsSA,LD"
# algo="BDT,BDTG"

# stages="0,10,1,2,9,8,4,5,6,7" ; sequence=1 ; # see definition below #
stages="0,10,1,2,9" ; sequence=0 ; # see definition below #

## ===== do not change the lines below =====

varlist=(
    '#  0  :  ("Bchi2cl"  , "Bchi2cl"                                                                                        , "FMax")  #' 
    '#  1  :  ("dRtrk1"   , "dRtrk1 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk1Phi)),2) + pow(Bujeta-Btrk1Eta,2))", "FMin")  #' 
    '#  2  :  ("dRtrk2"   , "dRtrk2 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk2Phi)),2) + pow(Bujeta-Btrk2Eta,2))", "FMin")  #' 
    '#  3  :  ("Qvalue"   , "Qvalue := (Bmass-3.096916-Btktkmass)"                                                           , "FMin")  #' 
    '#  4  :  ("Balpha"   , "Balpha"                                                                                         , "FMin")  #' 
    '#  5  :  ("costheta" , "costheta := TMath::Cos(Bdtheta)"                                                                , "FMax")  #' 
    '#  6  :  ("dls3D"    , "dls3D := TMath::Abs(BsvpvDistance/BsvpvDisErr)"                                                 , "FMax")  #' 
    '#  7  :  ("dls2D"    , "dls2D := TMath::Abs(BsvpvDistance_2D/BsvpvDisErr_2D)"                                           , "FMax")  #' 
    '#  8  :  ("Btrk1Pt"  , "Btrk1Pt"                                                                                        , "FMax")  #' 
    '#  9  :  ("Btrk2Pt"  , "Btrk2Pt"                                                                                        , "FMax")  #' 
    '#  10 :  ("trkptimba", "trkptimba := TMath::Abs((Btrk1Pt-Btrk2Pt) / (Btrk1Pt+Btrk2Pt))"                                 , "FMax")  #' 
    '#  11 :  ("By"       , "By"                                                                                             , ""    )  #'
    '#  12 :  ("Bmass"    , "Bmass"                                                                                          , ""    )  #'
    '#  13 :  ("Btrk1Eta" , "Btrk1Eta"                                                                                       , ""    )  #'
    '#  14 :  ("Btrk2Eta" , "Btrk2Eta"                                                                                       , ""    )  #'
)


cuts="${cut}&&Bgen==23333"
# cutb="${cut}&&((Bmass>3.74&&Bmass<3.83) || (Bmass>3.60&&Bmass<3.65) || (Bmass>3.93&&Bmass<4))"
cutb="${cut}&&((Bmass>3.74&&Bmass<3.83) || (Bmass>3.93&&Bmass<4))"
# rootfiles=rootfiles

## ===== do not do not do not change the lines below =====
IFS=','; allstages=($stages); unset IFS;
echo -e '

####################################################################################################################################
#                                                Variables \e[1;32m(To be used)\e[0m                                                            #
####################################################################################################################################
#                                                                                                                                  #'
vv=0
while(( $vv < ${#varlist[@]})) ; do
    for ss in ${allstages[@]} ; do [[ $ss == $vv ]] && { echo -en "\e[1;32m" ; break ; } ; done ;
    echo -e "${varlist[vv]}\e[0m"
    vv=$((vv+1))
done
echo '#                                                                                                                                  #
####################################################################################################################################

'
##
mkdir -p $output ; rm -r $output ;
tmp=$(date +%y%m%d%H%M%S)

##
[[ $# -eq 0 ]] && echo "usage: ./run_TMVAClassification.sh [train] [draw curves] [create BDT tree]"

g++ TMVAClassification.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o TMVAClassification_${tmp}.exe || exit 1
g++ guivariables.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guivariables_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
g++ guiefficiencies.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guiefficiencies_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
g++ guieffvar.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guieffvar_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
g++ mvaprod.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o mvaprod_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }

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
[[ ${3:-0} -eq 1 ]] && ./mvaprod_${tmp}.exe $inputm "Bfinder/ntmix" $output $outputmva $ptmin $ptmax "$algo" "${stages}"

##
rm *_${tmp}.exe