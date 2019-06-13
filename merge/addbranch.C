#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <string>
#include <regex>

#include "addbranch.h"

#define MAX_XB       20000

void xjjc::addbranch(std::string inputname, std::string outputname, std::string treename)
{
  char warning = 'x';
  std::cout<<std::endl<<" -- Warning..."<<std::endl;
  std::cout<<"*** This file will be changed: "<<std::endl<<"  "<<outputname<<std::endl;
  while(warning != 'n' && warning != 'y')
    {
      std::cout<<"*** Do you want to continue? [y/n]"<<std::endl;
      std::cin>>warning;
    }
  if(warning == 'n') { return; }

  std::cout<<std::endl<<" -- Checking input/output files"<<std::endl;
  std::cout<<"input file: "<<inputname<<std::endl;
  std::cout<<"output file: "<<outputname<<std::endl;
  std::cout<<"tree: "<<treename<<std::endl;

  std::cout<<std::endl<<" -- Building tree"<<std::endl;
  TFile* inf = TFile::Open(inputname.c_str());
  if(!inf->IsOpen()) { std::cout<<__FUNCTION__<<": error: invalid input."<<std::endl<<"  "<<inputname<<std::endl; return; }
  TTree* trread = (TTree*)inf->Get(treename.c_str());
  TFile* outf = TFile::Open(outputname.c_str(), "update");
  if(!outf->IsOpen()) { std::cout<<__FUNCTION__<<": error: invalid output."<<std::endl<<"  "<<outputname<<std::endl; return; }
  TTree* tr = (TTree*)outf->Get(treename.c_str());
  //
  std::vector<std::string> tdirs = xjjc::str_divide(treename, "/"); 
  std::string trname = tdirs[tdirs.size()-1];
  tdirs.pop_back();
  //
  outf->cd();
  if(!tr || !trread) { std::cout<<__FUNCTION__<<": error: tree does not exist."<<std::endl; return; }
  std::cout<<"TTree "<<treename<<" is opened."<<std::endl; 

  std::cout<<std::endl<<" -- Setting branch address"<<std::endl;
  Int_t Bsize; trread->SetBranchAddress("Bsize", &Bsize);
  Float_t Btrk1Pt[MAX_XB];        trread->SetBranchAddress("Btrk1Pt",        Btrk1Pt);
  Float_t Btrk1Eta[MAX_XB];       trread->SetBranchAddress("Btrk1Eta",       Btrk1Eta);
  Float_t Btrk1Phi[MAX_XB];       trread->SetBranchAddress("Btrk1Phi",       Btrk1Phi);
  Float_t Btrk1Dxy1[MAX_XB];      trread->SetBranchAddress("Btrk1Dxy1",      Btrk1Dxy1);
  Float_t Btrk1DxyError1[MAX_XB]; trread->SetBranchAddress("Btrk1DxyError1", Btrk1DxyError1);
  Float_t Btrk2Pt[MAX_XB];        trread->SetBranchAddress("Btrk2Pt",        Btrk2Pt);
  Float_t Btrk2Eta[MAX_XB];       trread->SetBranchAddress("Btrk2Eta",       Btrk2Eta);
  Float_t Btrk2Phi[MAX_XB];       trread->SetBranchAddress("Btrk2Phi",       Btrk2Phi);
  Float_t Btrk2Dxy1[MAX_XB];      trread->SetBranchAddress("Btrk2Dxy1",      Btrk2Dxy1);
  Float_t Btrk2DxyError1[MAX_XB]; trread->SetBranchAddress("Btrk2DxyError1", Btrk2DxyError1);

  std::cout<<std::endl<<" -- Building and filling branch"<<std::endl;
  Int_t   BsizeH;                 TBranch* br_BsizeH         = tr->Branch("BsizeH",         &BsizeH);
  Int_t   BtrkWhichH[MAX_XB];     TBranch* br_BtrkWhichH     = tr->Branch("BtrkWhichH",     BtrkWhichH,     "BtrkWhichH[BsizeH]/I");
  Float_t BtrkLPt[MAX_XB];        TBranch* br_BtrkLPt        = tr->Branch("BtrkLPt",        BtrkLPt,        "BtrkLPt[BsizeH]/F");
  Float_t BtrkLEta[MAX_XB];       TBranch* br_BtrkLEta       = tr->Branch("BtrkLEta",       BtrkLEta,       "BtrkLEta[BsizeH]/F");
  Float_t BtrkLPhi[MAX_XB];       TBranch* br_BtrkLPhi       = tr->Branch("BtrkLPhi",       BtrkLPhi,       "BtrkLPhi[BsizeH]/F");
  Float_t BtrkLDxy1[MAX_XB];      TBranch* br_BtrkLDxy1      = tr->Branch("BtrkLDxy1",      BtrkLDxy1,      "BtrkLDxy1[BsizeH]/F");
  Float_t BtrkLDxyError1[MAX_XB]; TBranch* br_BtrkLDxyError1 = tr->Branch("BtrkLDxyError1", BtrkLDxyError1, "BtrkLDxyError1[BsizeH]/F");
  Float_t BtrkHPt[MAX_XB];        TBranch* br_BtrkHPt        = tr->Branch("BtrkHPt",        BtrkHPt,        "BtrkHPt[BsizeH]/F");
  Float_t BtrkHEta[MAX_XB];       TBranch* br_BtrkHEta       = tr->Branch("BtrkHEta",       BtrkHEta,       "BtrkHEta[BsizeH]/F");
  Float_t BtrkHPhi[MAX_XB];       TBranch* br_BtrkHPhi       = tr->Branch("BtrkHPhi",       BtrkHPhi,       "BtrkHPhi[BsizeH]/F");
  Float_t BtrkHDxy1[MAX_XB];      TBranch* br_BtrkHDxy1      = tr->Branch("BtrkHDxy1",      BtrkHDxy1,      "BtrkHDxy1[BsizeH]/F");
  Float_t BtrkHDxyError1[MAX_XB]; TBranch* br_BtrkHDxyError1 = tr->Branch("BtrkHDxyError1", BtrkHDxyError1, "BtrkHDxyError1[BsizeH]/F");

  int nentries = trread->GetEntries();
  for(int i=0; i<nentries; i++) 
    {
      if(i%10000==0) xjjc::progressbar(i, nentries);
      trread->GetEntry(i);
      
      BsizeH = Bsize;
      for(int j=0; j<Bsize; j++)
        {
          BtrkWhichH[j] = 0;
          if(Btrk1Pt[j] > Btrk2Pt[j])
            {
              BtrkWhichH[j]     = 1;
              BtrkHPt[j]        = Btrk1Pt[j];
              BtrkHEta[j]       = Btrk1Eta[j];
              BtrkHPhi[j]       = Btrk1Phi[j];
              BtrkHDxy1[j]      = Btrk1Dxy1[j];
              BtrkHDxyError1[j] = Btrk1DxyError1[j];
              BtrkLPt[j]        = Btrk2Pt[j];
              BtrkLEta[j]       = Btrk2Eta[j];
              BtrkLPhi[j]       = Btrk2Phi[j];
              BtrkLDxy1[j]      = Btrk2Dxy1[j];
              BtrkLDxyError1[j] = Btrk2DxyError1[j];
            }
          else
            {
              BtrkWhichH[j]     = 2;
              BtrkHPt[j]        = Btrk2Pt[j];
              BtrkHEta[j]       = Btrk2Eta[j];
              BtrkHPhi[j]       = Btrk2Phi[j];
              BtrkHDxy1[j]      = Btrk2Dxy1[j];
              BtrkHDxyError1[j] = Btrk2DxyError1[j];
              BtrkLPt[j]        = Btrk1Pt[j];
              BtrkLEta[j]       = Btrk1Eta[j];
              BtrkLPhi[j]       = Btrk1Phi[j];
              BtrkLDxy1[j]      = Btrk1Dxy1[j];
              BtrkLDxyError1[j] = Btrk1DxyError1[j];
            }
        }
      br_BsizeH->Fill();
      br_BtrkWhichH->Fill();
      br_BtrkHPt->Fill();
      br_BtrkHEta->Fill();
      br_BtrkHPhi->Fill();
      br_BtrkHDxy1->Fill();
      br_BtrkHDxyError1->Fill();
      br_BtrkLPt->Fill();
      br_BtrkLEta->Fill();
      br_BtrkLPhi->Fill();
      br_BtrkLDxy1->Fill();
      br_BtrkLDxyError1->Fill();
    }
  xjjc::progressbar_summary(nentries);
  outf->cd();
  for(auto& dir : tdirs) { if(!gDirectory->cd(dir.c_str())) 
      { std::cout<<__FUNCTION__<<": error: did not manage to mkdir at beginning."<<std::endl; return; } }
  tr->Write("", TObject::kOverwrite);
  outf->Close();

  std::cout<<std::endl<<" -- End"<<std::endl<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc == 3)
    {
      std::string argvinput(argv[1]);
      std::string argvoutput = xjjc::str_replaceall(argvinput, ".root",
                                                    Form("__add_%s_%s.root", xjjc::str_replaceallspecial(argv[2]).c_str(), "BtrkWhichH"));
      std::cout<<std::endl<<" -- Create new file"<<std::endl;
      gSystem->Exec(Form("rsync --progress %s %s", argvinput.c_str(), argvoutput.c_str()));
      xjjc::addbranch(argvinput, argvoutput, argv[2]);
      return 0;
    }

  std::cout<<__FUNCTION__<<": ./addbranch.exe [inputfile] [treename]"<<std::endl;
  return 1;
}
