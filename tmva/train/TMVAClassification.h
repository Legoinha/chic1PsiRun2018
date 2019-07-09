#ifndef _TMVACLASSIFICATION_H_
#define _TMVACLASSIFICATION_H_

#include <string>
#include <vector>
#include <TString.h>
#include <TRandom.h>

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
    /*0: */ mytmva::tmvavar("Bchi2cl"    , "Bchi2cl"                                                                                        , "FMax", "vertex #chi^{2} prob"                 , 0   , 1)  ,
    /*1: */ mytmva::tmvavar("dRtrk1"     , "dRtrk1 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk1Phi)),2) + pow(Bujeta-Btrk1Eta,2))", "FMin", "#DeltaR(#pi_{1},J/#psi)"              , 0   , 0.5),
    /*2: */ mytmva::tmvavar("dRtrk2"     , "dRtrk2 := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-Btrk2Phi)),2) + pow(Bujeta-Btrk2Eta,2))", "FMin", "#DeltaR(#pi_{2},J/#psi)"              , 0   , 0.5),
    /*3: */ mytmva::tmvavar("Qvalue"     , "Qvalue := (Bmass-3.096916-Btktkmass)"                                                           , "FMin", "Q (GeV/c^{2})"                        , 0   , 0.5),
    /*4: */ mytmva::tmvavar("Balpha"     , "Balpha"                                                                                         , "FMin", "#alpha"                               , 0   , 3.2),
    /*5: */ mytmva::tmvavar("costheta"   , "costheta := TMath::Cos(Bdtheta)"                                                                , "FMax", "cos(#theta)"                          , -1  , 1)  ,
    /*6: */ mytmva::tmvavar("dls3D"      , "dls3D := TMath::Abs(BsvpvDistance/BsvpvDisErr)"                                                 , "FMax", "l_{xyz}/#sigma(l_{xyz})"              , 0   , 10) ,
    /*7: */ mytmva::tmvavar("dls2D"      , "dls2D := TMath::Abs(BsvpvDistance_2D/BsvpvDisErr_2D)"                                           , "FMax", "l_{xy}/#sigma(l_{xy})"                , 0   , 10) ,
    /*8: */ mytmva::tmvavar("Btrk1Pt"    , "Btrk1Pt"                                                                                        , "FMax", "#pi_{1} p_{T} (GeV/c)"                , 0   , 10) ,
    /*9: */ mytmva::tmvavar("Btrk2Pt"    , "Btrk2Pt"                                                                                        , "FMax", "#pi_{2} p_{T} (GeV/c)"                , 0   , 10) ,
    /*10:*/ mytmva::tmvavar("trkptimba"  , "trkptimba := TMath::Abs((Btrk1Pt-Btrk2Pt) / (Btrk1Pt+Btrk2Pt))"                                 , "FMax", "p_{T}(1)-p_{T}(2) / p_{T}(1)+p_{T}(2)", 0   , 1)  ,
    /*11:*/ mytmva::tmvavar("By"         , "By"                                                                                             , ""    , "y"                                    , -2.4, 2.4),
    /*12:*/ mytmva::tmvavar("Bmass"      , "Bmass"                                                                                          , ""    , "m_{#mu#mu#pi#pi}"                     , 3.6 , 4.0),
    /*13:*/ mytmva::tmvavar("Btrk1Eta"   , "Btrk1Eta"                                                                                       , ""    , "#pi_{1} #eta"                         , -2  , 2)  ,
    /*14:*/ mytmva::tmvavar("Btrk2Eta"   , "Btrk2Eta"                                                                                       , ""    , "#pi_{2} #eta"                         , -2  , 2)  ,
    /*15:*/ mytmva::tmvavar("Btrk1DxySig", "Btrk1DxySig := TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)"                                            , ""    , "#pi_{1} |D_{xy}/#sigma(D_{xy})|"      , 0   , 4)  ,
    /*16:*/ mytmva::tmvavar("Btrk2DxySig", "Btrk2DxySig := TMath::Abs(Btrk2Dxy1/Btrk2DxyError1)"                                            , ""    , "#pi_{2} |D_{xy}/#sigma(D_{xy})|"      , 0   , 4)  ,
    /*17:*/ mytmva::tmvavar("dRtrkH"     , "dRtrkH := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-BtrkHPhi)),2) + pow(Bujeta-BtrkHEta,2))", "FMin", "#DeltaR(#pi_{lead},J/#psi)"           , 0   , 0.5),
    /*18:*/ mytmva::tmvavar("dRtrkL"     , "dRtrkL := TMath::Sqrt(pow(TMath::ACos(TMath::Cos(Bujphi-BtrkLPhi)),2) + pow(Bujeta-BtrkLEta,2))", "FMin", "#DeltaR(#pi_{sub},J/#psi)"            , 0   , 0.5),
    /*19:*/ mytmva::tmvavar("BtrkHPt"    , "BtrkHPt"                                                                                        , "FMax", "#pi_{lead} p_{T} (GeV/c)"             , 0   , 10) ,
    /*20:*/ mytmva::tmvavar("BtrkLPt"    , "BtrkLPt"                                                                                        , "FMax", "#pi_{sub} p_{T} (GeV/c)"              , 0   , 10) ,
    /*21:*/ mytmva::tmvavar("BtrkHEta"   , "BtrkHEta"                                                                                       , ""    , "#pi_{lead} #eta"                      , -2  , 2)  ,
    /*22:*/ mytmva::tmvavar("BtrkLEta"   , "BtrkLEta"                                                                                       , ""    , "#pi_{sub} #eta"                       , -2  , 2)  ,
    /*23:*/ mytmva::tmvavar("BtrkHDxySig", "BtrkHDxySig := TMath::Abs(BtrkHDxy1/BtrkHDxyError1)"                                            , ""    , "#pi_{lead} |D_{xy}/#sigma(D_{xy})|"   , 0   , 4)  ,
    /*24:*/ mytmva::tmvavar("BtrkLDxySig", "BtrkLDxySig := TMath::Abs(BtrkLDxy1/BtrkLDxyError1)"                                            , ""    , "#pi_{sub} |D_{xy}/#sigma(D_{xy})|"    , 0   , 4)  ,
  };
  const mytmva::tmvavar* findvar(std::string varlabel);

  class varval
  {
  public:
    varval(mytmva::ntuple* nt) : fnt(nt), fvalid(true), rr(new TRandom()) { fval.clear(); fvalid = checkvarlist(); }
    varval(TTree* nttree) : fnt(0), fvalid(true), rr(new TRandom()) { fnt = new mytmva::ntuple(nttree); fval.clear(); fvalid = checkvarlist(); }
    float getval(std::string varname, int j) { refreshval(j); if(fval.find(varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<varname<<std::endl; return 0; } ; return fval[varname]; }
    mytmva::ntuple* getnt() { return fnt; }
    bool isvalid() { return fvalid; }

  private:
    bool fvalid;
    std::map<std::string, float> fval;
    mytmva::ntuple* fnt; //~
    TRandom* rr;

    void refreshval(int j)
    {
      bool ll = rr->Integer(2) < 0.5;
      fval["Bmass"]       = j<0?0:(float)fnt->Bmass[j];
      fval["dRtrkH"]      = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->BtrkHPhi[j])), 2) + pow(fnt->Bujeta[j] - fnt->BtrkHEta[j], 2));
      fval["dRtrkL"]      = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->BtrkLPhi[j])), 2) + pow(fnt->Bujeta[j] - fnt->BtrkLEta[j], 2));
      fval["Qvalue"]      = j<0?0:(float)(fnt->Bmass[j]-3.096916-fnt->Btktkmass[j]);
      fval["costheta"]    = j<0?0:(float)TMath::Cos(fnt->Bdtheta[j]);
      fval["dls3D"]       = j<0?0:(float)TMath::Abs(fnt->BsvpvDistance[j]/fnt->BsvpvDisErr[j]);
      fval["Bchi2cl"]     = j<0?0:(float)fnt->Bchi2cl[j];
      fval["BtrkHPt"]     = j<0?0:(float)fnt->BtrkHPt[j];
      fval["BtrkLPt"]     = j<0?0:(float)fnt->BtrkLPt[j];
      fval["Balpha"]      = j<0?0:(float)fnt->Balpha[j];
      fval["dls2D"]       = j<0?0:(float)TMath::Abs(fnt->BsvpvDistance_2D[j]/fnt->BsvpvDisErr_2D[j]);
      fval["trkptimba"]   = j<0?0:(float)TMath::Abs((fnt->Btrk1Pt[j]-fnt->Btrk2Pt[j]) / (fnt->Btrk1Pt[j]+fnt->Btrk2Pt[j]));
      fval["By"]          = j<0?0:(float)fnt->By[j];
      fval["BtrkHEta"]    = j<0?0:(float)fnt->BtrkHEta[j];
      fval["BtrkLEta"]    = j<0?0:(float)fnt->BtrkLEta[j];
      fval["BtrkHDxySig"] = j<0?0:(float)TMath::Abs(fnt->BtrkHDxy1[j]/fnt->BtrkHDxyError1[j]);
      fval["BtrkLDxySig"] = j<0?0:(float)TMath::Abs(fnt->BtrkLDxy1[j]/fnt->BtrkLDxyError1[j]);
      if(ll)
        {
          fval["dRtrk1"]      = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk1Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk1Eta[j], 2));
          fval["dRtrk2"]      = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk2Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk2Eta[j], 2));
          fval["Btrk1Pt"]     = j<0?0:(float)fnt->Btrk1Pt[j];
          fval["Btrk2Pt"]     = j<0?0:(float)fnt->Btrk2Pt[j];
          fval["Btrk1Eta"]    = j<0?0:(float)fnt->Btrk1Eta[j];
          fval["Btrk2Eta"]    = j<0?0:(float)fnt->Btrk2Eta[j];
          fval["Btrk1DxySig"] = j<0?0:(float)TMath::Abs(fnt->Btrk1Dxy1[j]/fnt->Btrk1DxyError1[j]);
          fval["Btrk2DxySig"] = j<0?0:(float)TMath::Abs(fnt->Btrk2Dxy1[j]/fnt->Btrk2DxyError1[j]);          
        }
      else
        {
          fval["dRtrk2"]      = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk1Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk1Eta[j], 2));
          fval["dRtrk1"]      = j<0?0:(float)TMath::Sqrt(pow(TMath::ACos(TMath::Cos(fnt->Bujphi[j] - fnt->Btrk2Phi[j])), 2) + pow(fnt->Bujeta[j] - fnt->Btrk2Eta[j], 2));
          fval["Btrk2Pt"]     = j<0?0:(float)fnt->Btrk1Pt[j];
          fval["Btrk1Pt"]     = j<0?0:(float)fnt->Btrk2Pt[j];
          fval["Btrk2Eta"]    = j<0?0:(float)fnt->Btrk1Eta[j];
          fval["Btrk1Eta"]    = j<0?0:(float)fnt->Btrk2Eta[j];
          fval["Btrk2DxySig"] = j<0?0:(float)TMath::Abs(fnt->Btrk1Dxy1[j]/fnt->Btrk1DxyError1[j]);
          fval["Btrk1DxySig"] = j<0?0:(float)TMath::Abs(fnt->Btrk2Dxy1[j]/fnt->Btrk2DxyError1[j]);
        }
    }
    bool checkvarlist() 
    {  
      refreshval(-1);
      for(auto& vn : varlist) 
        { if(fval.find(vn.varname) == fval.end()) { std::cout<<"==> "<<__FUNCTION__<<": invalid varname key "<<vn.varname<<std::endl; return false; } }
      return true;
    }
  };

  std::vector<std::string> argmethods; std::vector<int> argstages;
  std::string mkname(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage,
                     std::vector<std::string> &methods = argmethods, std::vector<int> &stages = argstages);
}

const mytmva::tmvavar* mytmva::findvar(std::string varlabel)
{
  for(auto& vv : varlist)
    {
      if(vv.varname == varlabel) return &vv;
    }
  return 0;
}

std::string mytmva::mkname(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage, 
                           std::vector<std::string> &methods, std::vector<int> &stages)
{
  mymethod = xjjc::str_replaceall(mymethod, " ", "");
  stage = xjjc::str_replaceall(stage, " ", "");
  methods = xjjc::str_divide(mymethod, ",");
  for(auto& ss : xjjc::str_divide(stage, ",")) { stages.push_back(atoi(ss.c_str())); }
  std::string outfname(Form("%s_%s_%s_%s_%s.root", outputname.c_str(),xjjc::str_replaceallspecial(mymethod).c_str(),
                            xjjc::number_to_string(ptmin).c_str(), (ptmax<0?"inf":xjjc::number_to_string(ptmax).c_str()),
                            xjjc::str_replaceall(stage, ",", "-").c_str()));
  return outfname;
}


#endif
