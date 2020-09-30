#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#define MAX_XB       20000
#define MAX_MUON     10000
#define MAX_TRIGGER  30

int skimmu(std::string inputname, std::string inputnamemuon, std::string outfile)
{
  TFile* inf = TFile::Open(inputname.c_str());
  TFile* infmuon = TFile::Open(inputnamemuon.c_str());

  std::cout<<inputname<<std::endl;
  TTree* nt = (TTree*)inf->Get("Bfinder/ntmix");
  int Bsize; nt->SetBranchAddress("Bsize", &Bsize);
  float Bmu1pt[MAX_XB]; nt->SetBranchAddress("Bmu1pt", Bmu1pt);
  float Bmu1eta[MAX_XB]; nt->SetBranchAddress("Bmu1eta", Bmu1eta);
  float Bmu1phi[MAX_XB]; nt->SetBranchAddress("Bmu1phi", Bmu1phi);
  float Bmu2pt[MAX_XB]; nt->SetBranchAddress("Bmu2pt", Bmu2pt);
  float Bmu2eta[MAX_XB]; nt->SetBranchAddress("Bmu2eta", Bmu2eta);
  float Bmu2phi[MAX_XB]; nt->SetBranchAddress("Bmu2phi", Bmu2phi);
  
  std::cout<<inputnamemuon<<std::endl;
  TTree* muinfo = (TTree*)infmuon->Get("muoninfo/root");
  int MuonInfo_size; muinfo->SetBranchAddress("MuonInfo.size", &MuonInfo_size);
  float MuonInfo_pt[MAX_MUON]; muinfo->SetBranchAddress("MuonInfo.pt", MuonInfo_pt);
  float MuonInfo_eta[MAX_MUON]; muinfo->SetBranchAddress("MuonInfo.eta", MuonInfo_eta);
  float MuonInfo_phi[MAX_MUON]; muinfo->SetBranchAddress("MuonInfo.phi", MuonInfo_phi);
  int MuonInfo_MuTrgMatchFilterSize; muinfo->SetBranchAddress("MuonInfo.MuTrgMatchFilterSize", &MuonInfo_MuTrgMatchFilterSize);
  float MuonInfo_MuTrgMatchFilterTrgObjE[MAX_MUON*MAX_TRIGGER]; muinfo->SetBranchAddress("MuonInfo.MuTrgMatchFilterTrgObjE", MuonInfo_MuTrgMatchFilterTrgObjE);

  std::cout<<outfile<<std::endl;
  TFile* outf = new TFile(outfile.c_str(), "recreate");
  TDirectory* dmuinfo = outf->mkdir("muinfo", "");
  dmuinfo->cd();
  TTree* ntmu = new TTree("mu", "");;
  int size; ntmu->Branch("size", &size, "size/I");
  int mu1match[MAX_XB]; ntmu->Branch("mu1match", mu1match, "mu1match[size]/I");
  float mu1pt[MAX_XB]; ntmu->Branch("mu1pt", mu1pt, "mu1pt[size]/F");
  float mu1eta[MAX_XB]; ntmu->Branch("mu1eta", mu1eta, "mu1eta[size]/F");
  float mu1phi[MAX_XB]; ntmu->Branch("mu1phi", mu1phi, "mu1phi[size]/F");
  int mu2match[MAX_XB]; ntmu->Branch("mu2match", mu2match, "m2match[size]/I");
  float mu2pt[MAX_XB]; ntmu->Branch("mu2pt", mu2pt, "mu2pt[size]/F");
  float mu2eta[MAX_XB]; ntmu->Branch("mu2eta", mu2eta, "mu2eta[size]/F");
  float mu2phi[MAX_XB]; ntmu->Branch("mu2phi", mu2phi, "mu2phi[size]/F");
  float mu1hltL2L3lastflt[MAX_XB]; ntmu->Branch("mu1hltL2L3lastflt", mu1hltL2L3lastflt, "mu1hltL2L3lastflt[size]/F");
  float mu1hltL2flt[MAX_XB]; ntmu->Branch("mu1hltL2flt", mu1hltL2flt, "mu1hltL2flt[size]/F");
  float mu1hltL3flt[MAX_XB]; ntmu->Branch("mu1hltL3flt", mu1hltL3flt, "mu1hltL3flt[size]/F");
  float mu2hltL2L3lastflt[MAX_XB]; ntmu->Branch("mu2hltL2L3lastflt", mu2hltL2L3lastflt, "mu2hltL2L3lastflt[size]/F");
  float mu2hltL2flt[MAX_XB]; ntmu->Branch("mu2hltL2flt", mu2hltL2flt, "mu2hltL2flt[size]/F");
  float mu2hltL3flt[MAX_XB]; ntmu->Branch("mu2hltL3flt", mu2hltL3flt, "mu2hltL3flt[size]/F");

  int nentries = nt->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      if(i%10000==0 || i==nentries) xjjc::progressbar(i, nentries);
      nt->GetEntry(i);
      muinfo->GetEntry(i);

      size = Bsize;
      for(int j=0; j<Bsize; j++)
        {
          int k1 = 0;
          for(;k1<MuonInfo_size; k1++)
            {
              if(fabs(Bmu1pt[j]-MuonInfo_pt[k1]) < 0.0001 &&
                 fabs(Bmu1eta[j]-MuonInfo_eta[k1]) < 0.0001 && 
                 fabs(Bmu1phi[j]-MuonInfo_phi[k1]) < 0.0001) break;
            }
          if(k1 < MuonInfo_size)
            {
              mu1match[j] = 1;
              mu1pt[j] = MuonInfo_pt[k1];
              mu1eta[j] = MuonInfo_eta[k1];
              mu1phi[j] = MuonInfo_phi[k1];
              mu1hltL2L3lastflt[j] = MuonInfo_MuTrgMatchFilterTrgObjE[MuonInfo_MuTrgMatchFilterSize*k1 + 0];
              mu1hltL2flt[j] = MuonInfo_MuTrgMatchFilterTrgObjE[MuonInfo_MuTrgMatchFilterSize*k1 + 1];
              mu1hltL3flt[j] = MuonInfo_MuTrgMatchFilterTrgObjE[MuonInfo_MuTrgMatchFilterSize*k1 + 2];
            }
          else 
            {
              std::cout<<i<<std::endl;
              std::cout<<"Bmu1pt["<<j<<"] = "<<Bmu1pt[j]<<std::endl;
              for(int k=0; k<MuonInfo_size; k++) { std::cout<<MuonInfo_pt[k]<<", "; }
              std::cout<<std::endl;
              mu1match[j] = 0;
              mu1pt[j] = -999;
              mu1eta[j] = -999;
              mu1phi[j] = -999;
              mu1hltL2L3lastflt[j] = -999;
              mu1hltL2flt[j] = -999;
              mu1hltL3flt[j] = -999;
              return 2;
            }

          //
          int k2 = 0;
          for(;k2<MuonInfo_size; k2++)
            {
              if(fabs(Bmu2pt[j]-MuonInfo_pt[k2]) < 0.0001 &&
                 fabs(Bmu2eta[j]-MuonInfo_eta[k2]) < 0.0001 && 
                 fabs(Bmu2phi[j]-MuonInfo_phi[k2]) < 0.0001) break;
            }
          if(k2 < MuonInfo_size)
            {
              mu2match[j] = 1;
              mu2pt[j] = MuonInfo_pt[k2];
              mu2eta[j] = MuonInfo_eta[k2];
              mu2phi[j] = MuonInfo_phi[k2];
              mu2hltL2L3lastflt[j] = MuonInfo_MuTrgMatchFilterTrgObjE[MuonInfo_MuTrgMatchFilterSize*k2 + 0];
              mu2hltL2flt[j] = MuonInfo_MuTrgMatchFilterTrgObjE[MuonInfo_MuTrgMatchFilterSize*k2 + 1];
              mu2hltL3flt[j] = MuonInfo_MuTrgMatchFilterTrgObjE[MuonInfo_MuTrgMatchFilterSize*k2 + 2];
            }
          else 
            {
              mu2match[j] = 0;
              mu2pt[j] = -999;
              mu2eta[j] = -999;
              mu2phi[j] = -999;
              mu2hltL2L3lastflt[j] = -999;
              mu2hltL2flt[j] = -999;
              mu2hltL3flt[j] = -999;
            }
        }
      ntmu->Fill();
    }
  xjjc::progressbar_summary(nentries);

  ntmu->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==4) { return skimmu(argv[1], argv[2], argv[3]); }
  return 1;
}
