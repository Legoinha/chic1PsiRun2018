#ifndef _TMVACLASSIFICATION_H_
#define _TMVACLASSIFICATION_H_

#include <string>
#include <vector>
#include <TString.h>

#include "xjjcuti.h"
#include "ntuple.h"

namespace mytmva
{
  struct tmvavar
  {
    std::string varname;
    std::string vartex;
    std::string var;
    std::string cutsign;
    float varmin;
    float varmax;
    tmvavar(const std::string& varname_, const std::string& var_, const std::string& cutsign_, const std::string& vartex_, const float& varmin_, const float& varmax_) 
      : varname(varname_), var(var_), cutsign(cutsign_), vartex(vartex_), varmin(varmin_), varmax(varmax_) { ; }
  };

  const std::vector<mytmva::tmvavar> varlist = {
    mytmva::tmvavar("Bchi2cl"  , "Bchi2cl"                                                                                        , "FMax", "vertex #chi^{2} prob"                 , 0   , 1)  , // 0
    mytmva::tmvavar("dRtrk1"   , "dRtrk1 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk1Phi)),2) + pow(Bujeta-Btrk1Eta,2))", "FMin", "#DeltaR(#pi_{1},J/#psi)"              , 0   , 0.5), // 1
    mytmva::tmvavar("dRtrk2"   , "dRtrk2 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk2Phi)),2) + pow(Bujeta-Btrk2Eta,2))", "FMin", "#DeltaR(#pi_{2},J/#psi)"              , 0   , 0.5), // 2
    mytmva::tmvavar("Qvalue"   , "Qvalue := (Bmass-3.096916-Btktkmass)"                                                           , "FMin", "Q (GeV/c^{2})"                        , 0   , 0.5), // 3
    mytmva::tmvavar("Balpha"   , "Balpha"                                                                                         , "FMin", "#alpha"                               , 0   , 3.2), // 4
    mytmva::tmvavar("costheta" , "costheta := TMath::Cos(Bdtheta)"                                                                , "FMax", "cos(#theta)"                          , -1  , 1)  , // 5
    mytmva::tmvavar("dls3D"    , "dls3D := TMath::Abs(BsvpvDistance/BsvpvDisErr)"                                                 , "FMax", "l_{xyz}/#sigma(l_{xyz})"              , 0   , 10) , // 6
    mytmva::tmvavar("dls2D"    , "dls2D := TMath::Abs(BsvpvDistance_2D/BsvpvDisErr_2D)"                                           , "FMax", "l_{xy}/#sigma(l_{xy})"                , 0   , 10) , // 7
    mytmva::tmvavar("Btrk1Pt"  , "Btrk1Pt"                                                                                        , "FMax", "#pi_{1} p_{T} (GeV/c)"                , 0   , 10) , // 8
    mytmva::tmvavar("Btrk2Pt"  , "Btrk2Pt"                                                                                        , "FMax", "#pi_{2} p_{T} (GeV/c)"                , 0   , 10) , // 9
    mytmva::tmvavar("trkptimba", "trkptimba := TMath::Abs((Btrk1Pt-Btrk2Pt) / (Btrk1Pt+Btrk2Pt))"                                 , "FMax", "p_{T}(1)-p_{T}(2) / p_{T}(1)+p_{T}(2)", 0   , 1)  , // 10
    mytmva::tmvavar("By"       , "By"                                                                                             , ""    , "y"                                    , -2.4, 2.4), // 11
    mytmva::tmvavar("Bmass"    , "Bmass"                                                                                          , ""    , "m_{#mu#mu#pi#pi}"                     , 3.6 , 4.0), // 12
    mytmva::tmvavar("Btrk1Eta" , "Btrk1Eta"                                                                                       , ""    , "#pi_{1} #eta"                         , -2  , 2)  , // 13
    mytmva::tmvavar("Btrk2Eta" , "Btrk2Eta"                                                                                       , ""    , "#pi_{2} #eta"                         , -2  , 2)  , // 14
  };

  class varval
  {
  public:
    varval(mytmva::ntuple* nt) : fnt(nt), fvalid(true) { fval.clear(); fvalid = checkvarlist(); }
    varval(TTree* nttree) : fnt(0), fvalid(true) { fnt = new mytmva::ntuple(nttree); fval.clear(); fvalid = checkvarlist(); }
    float getval(std::string varname, int j) { refreshval(j); if(fval.find(varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<varname<<std::endl; return 0; } ; return fval[varname]; }
    mytmva::ntuple* getnt() { return fnt; }
    bool isvalid() { return fvalid; }

  private:
    bool fvalid;
    std::map<std::string, float> fval;
    mytmva::ntuple* fnt; //~
    void refreshval(int j)
    {
      fval["Bmass"]     = j<0?0:(float)fnt->Bmass[j];
      fval["dRtrk1"]    = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk1Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk1Eta[j], 2));
      fval["dRtrk2"]    = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk2Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk2Eta[j], 2));
      fval["Qvalue"]    = j<0?0:(float)(fnt->Bmass[j]-3.096916-fnt->Btktkmass[j]);
      fval["costheta"]  = j<0?0:(float)TMath::Cos(fnt->Bdtheta[j]);
      fval["dls3D"]     = j<0?0:(float)TMath::Abs(fnt->BsvpvDistance[j]/fnt->BsvpvDisErr[j]);
      fval["Bchi2cl"]   = j<0?0:(float)fnt->Bchi2cl[j];
      fval["Btrk1Pt"]   = j<0?0:(float)fnt->Btrk1Pt[j];
      fval["Btrk2Pt"]   = j<0?0:(float)fnt->Btrk2Pt[j];
      fval["Balpha"]    = j<0?0:(float)fnt->Balpha[j];
      fval["Btrk1dca"]  = j<0?0:(float)TMath::Abs(fnt->Btrk1Dxy[j]/fnt->Btrk1D0Err[j]);
      fval["Btrk2dca"]  = j<0?0:(float)TMath::Abs(fnt->Btrk2Dxy[j]/fnt->Btrk2D0Err[j]);
      fval["dls2D"]     = j<0?0:(float)TMath::Abs(fnt->BsvpvDistance_2D[j]/fnt->BsvpvDisErr_2D[j]);
      fval["trkptimba"] = j<0?0:(float)TMath::Abs((fnt->Btrk1Pt[j]-fnt->Btrk2Pt[j]) / (fnt->Btrk1Pt[j]+fnt->Btrk2Pt[j]));
      fval["By"]        = j<0?0:(float)fnt->By[j];
      fval["Btrk1Eta"]  = j<0?0:(float)fnt->Btrk1Eta[j];
      fval["Btrk2Eta"]  = j<0?0:(float)fnt->Btrk2Eta[j];
    }
    bool checkvarlist() 
    {  
      refreshval(-1);
      for(auto& vn : varlist) 
        { if(fval.find(vn.varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<vn.varname<<std::endl; return false; } }
      return true;
    }
  };

}

#endif
