#ifndef _VARCOMP_NTUPLE_H_
#define _VARCOMP_NTUPLE_H_

#include <TTree.h>
#include <TMath.h>

#ifndef MAX_XB
#define MAX_XB      20000
#endif

#ifndef MAX_GEN
#define MAX_GEN     6000
#endif

#define PDGID_JPSI  443
#define MASS_JPSI   3.096916
#define PDGID_X     20443
#define MASS_X      3.87169
#define PDGID_PSI2S 100443
#define MASS_PSI2S  3.686097

namespace mytmva
{
  class ntuple
  {
  public:
    ntuple(TTree* nt, TTree* gnt) : fnt(nt), fgnt(gnt) { setbranchaddress(); }
    ntuple(TTree* nt) : fnt(nt) { fgnt = 0; setbranchaddress(); }
    ~ntuple() { fnt = 0; fgnt = 0; }
    bool passedpre(int j);
    bool passedsig(int j) { return (Bgen[j]==23333 && BgencollisionId[j]==0); }
    bool passedevtfil() { return (HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 && pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && hiBin >=0 && hiBin < 180); }
    // bool passedbkg(int j) { return ((Bmass[j]>3.74 && Bmass[j]<3.82) || (Bmass[j]>3.60 && Bmass[j]<3.65) || (Bmass[j]>3.92 && Bmass[j]<4.00)); }
    // bool passedbkg(int j) { return ((Bmass[j]>3.74 && Bmass[j]<3.83) || (Bmass[j]>3.93 && Bmass[j]<4.00)); }
    bool passedbkg(int j) { return (TMath::Abs(Bmass[j]-3.8719) > 0.07 && TMath::Abs(Bmass[j]-3.8719) < 0.128); }
    bool isweight() { return fweight; }
    bool ishlt() { return fhlt; }
    bool isskim() { return fskim; }

    // >>>
    bool signalregionl(int j) { return (TMath::Abs(Bmass[j]-MASS_PSI2S) < 0.018); }
    bool signalregionh(int j) { return (TMath::Abs(Bmass[j]-MASS_X) < 0.018); }
    bool sidebandl(int j) { return (TMath::Abs(Bmass[j]-MASS_PSI2S) > 0.04 && TMath::Abs(Bmass[j]-MASS_PSI2S) < 0.12); }
    bool sidebandh(int j) { return (TMath::Abs(Bmass[j]-MASS_X) > 0.02 && TMath::Abs(Bmass[j]-MASS_X) < 0.10); }

    TTree* getnt() { return fnt; }
    TTree* getgnt() { return fgnt; }

    float  pthatweight;
    float  Ncoll;
    int    hiBin;
    int    HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;
    int    pprimaryVertexFilter;
    int    phfCoincFilter2Th4;
    int    pclusterCompatibilityFilter;
    float  BDT[MAX_XB];
    float  BDTG[MAX_XB];
    bool   mvapref[MAX_XB];

    int    Bsize;
    float  Bgen[MAX_XB];
    int    BgencollisionId[MAX_XB];
    float  Bpt[MAX_XB];
    float  Balpha[MAX_XB];
    float  Blxy[MAX_XB];
    float  Btrk1Pt[MAX_XB];
    float  Btrk2Pt[MAX_XB];
    float  Bmass[MAX_XB];
    float  Btktkmass[MAX_XB];
    float  Bujpt[MAX_XB];
    float  Bchi2cl[MAX_XB];
    float  Bdtheta[MAX_XB];
    float  Bujphi[MAX_XB];
    float  Btrk1Phi[MAX_XB];
    float  Btrk2Phi[MAX_XB];
    float  BsvpvDistance[MAX_XB];
    float  BsvpvDisErr[MAX_XB];
    float  BsvpvDistance_2D[MAX_XB];
    float  BsvpvDisErr_2D[MAX_XB];
    float  Bujeta[MAX_XB];
    float  Btrk1Eta[MAX_XB];
    float  Btrk2Eta[MAX_XB];
    float  Btrk1Dxy1[MAX_XB];
    float  Btrk2Dxy1[MAX_XB];
    float  Btrk1DxyError1[MAX_XB];
    float  Btrk2DxyError1[MAX_XB];
    float  Btrk1Chi2ndf[MAX_XB];
    float  Btrk2Chi2ndf[MAX_XB];
    float  Btrk1nStripLayer[MAX_XB];
    float  Btrk2nStripLayer[MAX_XB];
    float  Btrk1nPixelLayer[MAX_XB];
    float  Btrk2nPixelLayer[MAX_XB];
    float  By[MAX_XB];

    //
    bool   Bmu1TMOneStationTight[MAX_XB];
    int    Bmu1InPixelLayer[MAX_XB];
    int    Bmu1InStripLayer[MAX_XB];
    float  Bmu1dxyPV[MAX_XB];
    float  Bmu1dzPV[MAX_XB];
    bool   Bmu1isGlobalMuon[MAX_XB];
    float  Bmu1eta[MAX_XB];
    float  Bmu1pt[MAX_XB];
    bool   Bmu2TMOneStationTight[MAX_XB];
    int    Bmu2InPixelLayer[MAX_XB];
    int    Bmu2InStripLayer[MAX_XB];
    float  Bmu2dxyPV[MAX_XB];
    float  Bmu2dzPV[MAX_XB];
    bool   Bmu2isGlobalMuon[MAX_XB];
    float  Bmu2eta[MAX_XB];
    float  Bmu2pt[MAX_XB];
    float  Bmumumass[MAX_XB];
    bool   Btrk1highPurity[MAX_XB];
    float  Btrk1PixelHit[MAX_XB];
    float  Btrk1StripHit[MAX_XB];
    float  Btrk1PtErr[MAX_XB];
    bool   Btrk2highPurity[MAX_XB];
    float  Btrk2PixelHit[MAX_XB];
    float  Btrk2StripHit[MAX_XB];
    float  Btrk2PtErr[MAX_XB];

    bool   Bmu1isTriggered[MAX_XB];
    bool   Bmu2isTriggered[MAX_XB];
    bool   Bmu1SoftMuID[MAX_XB];
    bool   Bmu2SoftMuID[MAX_XB];
    bool   Bmu1isAcc[MAX_XB];
    bool   Bmu2isAcc[MAX_XB];

    int    Gsize;
    float  Gy[MAX_GEN];
    float  Gpt[MAX_GEN];
    int    GpdgId[MAX_GEN];
    int    GcollisionId[MAX_GEN];
    int    GisSignal[MAX_GEN];

  private:
    TTree* fnt;
    TTree* fgnt;
    bool   fweight;
    bool   fskim;
    bool   fhlt;
    bool   fmva;
    bool   fhi;
    void setbranchaddress();

  };
}

bool mytmva::ntuple::passedpre(int j)
{
  if(
     Bpt[j] > 10 && TMath::Abs(Bmumumass[j]-3.096916) < 0.05 && TMath::Abs(Bujeta[j]) < 2.4 &&
     (Bmu1SoftMuID[j] && Bmu2SoftMuID[j] && Bmu1isAcc[j] && Bmu2isAcc[j] && Bmu1isTriggered[j] && Bmu2isTriggered[j]) &&
     Btrk1Pt[j] > 0.9 && Btrk2Pt[j] > 0.9 && TMath::Abs(Btrk1Eta[j]) < 2.4 && TMath::Abs(Btrk2Eta[j]) < 2.4 &&
     Btrk1highPurity[j] && Btrk2highPurity[j] && (Btrk1PixelHit[j]+Btrk1StripHit[j]) >= 11 && (Btrk2PixelHit[j]+Btrk2StripHit[j]) >= 11 && TMath::Abs(Btrk1PtErr[j]/Btrk1Pt[j]) < 0.1 && TMath::Abs(Btrk2PtErr[j]/Btrk2Pt[j]) < 0.1 && (Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j])) < 0.18 && (Btrk2Chi2ndf[j]/(Btrk2nStripLayer[j]+Btrk2nPixelLayer[j])) < 0.18 && 
     // TMath::Abs(Btrk1Dxy1[j]/Btrk1Dxy1[j]) < 3 && TMath::Abs(Btrk2Dxy1[j]/Btrk2DxyError1[j]) < 3 &&
     TMath::Abs(By[j]) < 2.4 && Bchi2cl[j] > 0.1 &&
     BsvpvDisErr[j] > 1.e-5 && BsvpvDisErr_2D[j] > 1.e-5 && 
     (Bmass[j]-3.096916-Btktkmass[j]) < 0.2
     ) return true;
  // if( Bmu1TMOneStationTight[j] && Bmu1InPixelLayer[j] > 1 && (Bmu1InPixelLayer[j]+Bmu1InStripLayer[j]) > 6 && Bmu1dxyPV[j] < 0.3 && Bmu1dzPV[j] < 20 && Bmu1isGlobalMuon[j] && TMath::Abs(Bmu1eta[j]) < 2 && Bmu1pt[j] > 1.5 && // mu1
  //     Bmu2TMOneStationTight[j] && Bmu2InPixelLayer[j] > 1 && (Bmu2InPixelLayer[j]+Bmu2InStripLayer[j]) > 6 && Bmu2dxyPV[j] < 0.3 && Bmu2dzPV[j] < 20 && Bmu2isGlobalMuon[j] && TMath::Abs(Bmu2eta[j]) < 2 && Bmu2pt[j] > 1.5 && // mu2
  //     TMath::Abs(Bmumumass[j]-3.096916) < 0.05 && TMath::Abs(Bujeta[j]) < 2.0 && // jpsi
  //     Btrk1highPurity[j] &&  TMath::Abs(Btrk1Eta[j]) < 2 && Btrk1Pt[j] > 0.9 && (Btrk1PixelHit[j]+Btrk1StripHit[j]) > 11 && TMath::Abs(Btrk1PtErr[j]/Btrk1Pt[j]) < 0.1 && // trk1
  //     Btrk2highPurity[j] &&  TMath::Abs(Btrk2Eta[j]) < 2 && Btrk2Pt[j] > 0.9 && (Btrk2PixelHit[j]+Btrk2StripHit[j]) > 11 && TMath::Abs(Btrk2PtErr[j]/Btrk2Pt[j]) < 0.1 && // trk2
  //     // (Bmass[j]-3.096916-Btktkmass[j]) < 0.2 && // Btktkmass[j] > 0.47 && // ! tktkmass
  //     TMath::Abs(By[j]) < 2.0 && Bchi2cl[j] > 0.1 // B
  //     ) return true;
  return false;
  /*
   */
}

void mytmva::ntuple::setbranchaddress()
{
  fweight = fnt->FindBranch("pthatweight") && fnt->FindBranch("Ncoll");
  // if(!fweight) { std::cout<<__FUNCTION__<<": warning: no branch \e[4mpthatweight\e[0m or \e[4mNcoll\e[0m."<<std::endl; }
  fhlt = fnt->FindBranch("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1");
  if(!fhlt) { std::cout<<__FUNCTION__<<": \e[31mwarning:\e[0m no branch \e[4m"<<"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1"<<"\e[0m"
                       <<"."<<std::endl; }
  fskim = fnt->FindBranch("pprimaryVertexFilter") && fnt->FindBranch("phfCoincFilter2Th4") && fnt->FindBranch("pclusterCompatibilityFilter");
  if(!fskim) { std::cout<<__FUNCTION__<<": \e[31mwarning:\e[0m no branch \e[4m"<<"pprimaryVertexFilter"<<"\e[0m"
                        <<" or \e[4m"<<"phfCoincFilter2Th4"<<"\e[0m"
                        <<" or \e[4m"<<"pclusterCompatibilityFilter"<<"\e[0m"
                        <<"."<<std::endl; }
  fmva = fnt->FindBranch("mvapref");
  if(!fmva) { std::cout<<__FUNCTION__<<": \e[31mwarning:\e[0m no branch \e[4m"<<"mvapref"<<"\e[0m"
                       <<"."<<std::endl; }
  fhi = fnt->FindBranch("hiBin");
  // public:
  if(fweight) 
    { 
      fnt->SetBranchAddress("pthatweight", &pthatweight); 
      fnt->SetBranchAddress("Ncoll", &Ncoll); 
    }
  if(fhlt) 
    { 
      fnt->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1", &HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1); 
    }
  if(fskim)
    {
      fnt->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter); 
      fnt->SetBranchAddress("phfCoincFilter2Th4", &phfCoincFilter2Th4); 
      fnt->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);       
    }
  if(fmva)
    {
      fnt->SetBranchAddress("mvapref", mvapref);
      fnt->SetBranchAddress("BDT", BDT);
      fnt->SetBranchAddress("BDTG", BDTG);
    }
  if(fhi)
    {
      fnt->SetBranchAddress("hiBin", &hiBin);
    }

  fnt->SetBranchAddress("Bsize", &Bsize);
  fnt->SetBranchAddress("Bgen", Bgen);
  fnt->SetBranchAddress("BgencollisionId", BgencollisionId);
  fnt->SetBranchAddress("Bpt", Bpt);
  fnt->SetBranchAddress("Balpha", Balpha);
  fnt->SetBranchAddress("Blxy", Blxy);
  fnt->SetBranchAddress("Btrk1Pt", Btrk1Pt);
  fnt->SetBranchAddress("Btrk2Pt", Btrk2Pt);
  fnt->SetBranchAddress("Bmass", Bmass);
  fnt->SetBranchAddress("Btktkmass", Btktkmass);
  fnt->SetBranchAddress("Bujpt", Bujpt);
  fnt->SetBranchAddress("Bchi2cl", Bchi2cl);
  fnt->SetBranchAddress("Bdtheta", Bdtheta);
  fnt->SetBranchAddress("Bujphi", Bujphi);
  fnt->SetBranchAddress("Btrk1Phi", Btrk1Phi);
  fnt->SetBranchAddress("Btrk2Phi", Btrk2Phi);
  fnt->SetBranchAddress("BsvpvDistance", BsvpvDistance);
  fnt->SetBranchAddress("BsvpvDisErr", BsvpvDisErr);
  fnt->SetBranchAddress("BsvpvDistance_2D", BsvpvDistance_2D);
  fnt->SetBranchAddress("BsvpvDisErr_2D", BsvpvDisErr_2D);
  fnt->SetBranchAddress("Bujeta", Bujeta);
  fnt->SetBranchAddress("Btrk1Eta", Btrk1Eta);
  fnt->SetBranchAddress("Btrk2Eta", Btrk2Eta);
  fnt->SetBranchAddress("Btrk1Dxy1", Btrk1Dxy1);
  fnt->SetBranchAddress("Btrk2Dxy1", Btrk2Dxy1);
  fnt->SetBranchAddress("Btrk1DxyError1", Btrk1DxyError1);
  fnt->SetBranchAddress("Btrk2DxyError1", Btrk2DxyError1);
  fnt->SetBranchAddress("Btrk1Chi2ndf", Btrk1Chi2ndf);
  fnt->SetBranchAddress("Btrk2Chi2ndf", Btrk2Chi2ndf);
  fnt->SetBranchAddress("Btrk1nStripLayer", Btrk1nStripLayer);
  fnt->SetBranchAddress("Btrk2nStripLayer", Btrk2nStripLayer);
  fnt->SetBranchAddress("Btrk1nPixelLayer", Btrk1nPixelLayer);
  fnt->SetBranchAddress("Btrk2nPixelLayer", Btrk2nPixelLayer);
  // private
  fnt->SetBranchAddress("Bmu1TMOneStationTight", Bmu1TMOneStationTight);
  fnt->SetBranchAddress("Bmu1SoftMuID", Bmu1SoftMuID);
  fnt->SetBranchAddress("Bmu1isAcc", Bmu1isAcc);
  fnt->SetBranchAddress("Bmu1isTriggered", Bmu1isTriggered);
  fnt->SetBranchAddress("Bmu1InPixelLayer", Bmu1InPixelLayer);
  fnt->SetBranchAddress("Bmu1InStripLayer", Bmu1InStripLayer);
  fnt->SetBranchAddress("Bmu1dxyPV", Bmu1dxyPV);
  fnt->SetBranchAddress("Bmu1dzPV", Bmu1dzPV);
  fnt->SetBranchAddress("Bmu1isGlobalMuon", Bmu1isGlobalMuon);
  fnt->SetBranchAddress("Bmu1eta", Bmu1eta);
  fnt->SetBranchAddress("Bmu1pt", Bmu1pt);
  fnt->SetBranchAddress("Bmu2TMOneStationTight", Bmu2TMOneStationTight);
  fnt->SetBranchAddress("Bmu2SoftMuID", Bmu2SoftMuID);
  fnt->SetBranchAddress("Bmu2isAcc", Bmu2isAcc);
  fnt->SetBranchAddress("Bmu2isTriggered", Bmu2isTriggered);
  fnt->SetBranchAddress("Bmu2InPixelLayer", Bmu2InPixelLayer);
  fnt->SetBranchAddress("Bmu2InStripLayer", Bmu2InStripLayer);
  fnt->SetBranchAddress("Bmu2dxyPV", Bmu2dxyPV);
  fnt->SetBranchAddress("Bmu2dzPV", Bmu2dzPV);
  fnt->SetBranchAddress("Bmu2isGlobalMuon", Bmu2isGlobalMuon);
  fnt->SetBranchAddress("Bmu2eta", Bmu2eta);
  fnt->SetBranchAddress("Bmu2pt", Bmu2pt);
  fnt->SetBranchAddress("Bmumumass", Bmumumass);
  fnt->SetBranchAddress("Btrk1highPurity", Btrk1highPurity);
  fnt->SetBranchAddress("Btrk1PixelHit", Btrk1PixelHit);
  fnt->SetBranchAddress("Btrk1StripHit", Btrk1StripHit);
  fnt->SetBranchAddress("Btrk1PtErr", Btrk1PtErr);
  fnt->SetBranchAddress("Btrk2highPurity", Btrk2highPurity);
  fnt->SetBranchAddress("Btrk2PixelHit", Btrk2PixelHit);
  fnt->SetBranchAddress("Btrk2StripHit", Btrk2StripHit);
  fnt->SetBranchAddress("Btrk2PtErr", Btrk2PtErr);
  fnt->SetBranchAddress("By", By);
  fnt->SetBranchAddress("Btktkmass", Btktkmass);

  if(fgnt)
    {
      fgnt->SetBranchAddress("Gsize", &Gsize);
      fgnt->SetBranchAddress("Gy", Gy);
      fgnt->SetBranchAddress("Gpt", Gpt);
      fgnt->SetBranchAddress("GpdgId", GpdgId);
      fgnt->SetBranchAddress("GcollisionId", GcollisionId);
      fgnt->SetBranchAddress("GisSignal", GisSignal);
    }
}

#endif
