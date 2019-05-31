#!/bin/bash

# input=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntmix_Bpt10_TMVA_X_10p0_inf_0-10-1-2-9.root
# input=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntmix_Bpt15_X_tktkmass2_15p0_inf_0-10-1-2-9.root
input=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntmix_Bpt10_X_tktk0p2_10p0_inf_0-10-1-2-9.root
inputmc_a=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_Psi2SToJpsiPiPi_prompt_20181231_pt5tkpt0p7dls0_pthatweight_X_tktk0p2_10p0_inf_0-10-1-2-9.root
inputmc_b=/export/d00/scratch/jwang/BntupleRun2018/crab_Bfinder_20190115_Hydjet_Pythia8_X3872ToJpsiRho_prompt_20181231_pt5tkpt0p7dls0_pthatweight_X_tktk0p2_10p0_inf_0-10-1-2-9.root

# cut="Bpt>15&&mvapref&&BDTG>0.78&&BsvpvDistance/BsvpvDisErr>0.8&&(Bmass-3.096916-Btktkmass)<0.4&&Btktkmass>0.45" ; name="tighttktkmass" ;
# cut="Bpt>15&&mvapref&&BDTG>0.85&&BsvpvDistance/BsvpvDisErr>0.8&&(Bmass-3.096916-Btktkmass)<0.2" ; name="rltktkmass" ;
# cut="Bpt>15&&mvapref&&(Bmass-3.096916-Btktkmass)<0.2&&BDTG>0.65&&BsvpvDistance/BsvpvDisErr>0.8" ; name="tktkmass2" ;
cut="Bpt>15&&mvapref&&(Bmass-3.096916-Btktkmass)<0.2&&BDTG>0.76&&BsvpvDistance/BsvpvDisErr>0.8" ; name="tktk0p2" ;
cutgen="Gpt>15&&TMath::Abs(Gy)<2.0&&GisSignal==7"

mkdir -p plots

# g++ fitX.C $(root-config --libs --cflags) -g -o fitX.exe || exit 1
# ./fitX.exe $input "$cut"
# rm fitX.exe

g++ fitXmc.C $(root-config --libs --cflags) -g -o fitXmc.exe || exit 1
g++ MCefficiency.C $(root-config --libs --cflags) -g -o MCefficiency.exe || exit 1


[[ ${1:-0} -eq 1 ]] && ./fitXmc.exe $input $inputmc_a $inputmc_b "${cut}" $name
[[ ${2:-0} -eq 1 ]] && ./MCefficiency.exe $inputmc_a $inputmc_b "${cut}" "$cutgen" $name


rm MCefficiency.exe
rm fitXmc.exe

