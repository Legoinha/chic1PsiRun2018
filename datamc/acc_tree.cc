#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

#include "ntuple.h"
#include "packtree.h"
#include "tree_flatten.h"

#define MAX_GEN      6000

void acc_tree(std::string input, std::string output="")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if(output=="") { output = xjjc::str_replaceall(input, ".root", "_skimacc.root"); }
  std::cout<<output<<std::endl;
  TFile* inf = TFile::Open(input.c_str());
  if(!inf->IsOpen()) return;
  xjjroot::packtree* ptr = new xjjroot::packtree(inf, "Bfinder/ntmix", "mc", "Bfinder/ntGen");
  mytmva::ntuple* nt = ptr->ntp();

  std::map<std::string, fitX::br_element<float, MAX_GEN>*> branchesf = {
    {"Gpt", new fitX::br_element<float, MAX_GEN>("Gpt")},
    {"Gy", new fitX::br_element<float, MAX_GEN>("Gy")},
    {"Gmu1pt", new fitX::br_element<float, MAX_GEN>("Gmu1pt")},
    {"Gmu2pt", new fitX::br_element<float, MAX_GEN>("Gmu2pt")},
    {"Gmu1eta", new fitX::br_element<float, MAX_GEN>("Gmu1eta")},
    {"Gmu2eta", new fitX::br_element<float, MAX_GEN>("Gmu2eta")},
    {"weight", new fitX::br_element<float, MAX_GEN>("weight")},
  };
  std::map<std::string, fitX::br_element<bool, MAX_GEN>*> brancheso = {
    {"isAcc", new fitX::br_element<bool, MAX_GEN>("isAcc")},
  };

  TFile* outf = new TFile(output.c_str(), "recreate");
  outf->cd();
  TTree* outnt = new TTree("gnt", "");
  for(auto& br : branchesf) { outnt->Branch(br.second->name.c_str(), &(br.second->value)); }
  for(auto& br : brancheso) { outnt->Branch(br.second->name.c_str(), &(br.second->value)); }

  int nentries = nt->gnt()->GetEntries();
  for(int i=0; i<nentries; i++) // 
    {
      if(i%1000==0) { xjjc::progressbar(i, nentries); }
      ptr->getentry(i);

      if(!(nt->hiBin >= fitX::centmincut && nt->hiBin < fitX::centmaxcut)) continue;
      if(!(nt->passedevtfil() && fabs(nt->PVz)<15)) continue;

      branchesf["weight"]->value = nt->pthatweight * nt->Ncoll;
      for(int j=0; j<nt->Gsize; j++)
        {
          if(!(nt->GisSignal[j]==7 && nt->GcollisionId[j]==0 
               && nt->Gpt[j] > fitX::ptmincut && nt->Gpt[j] < fitX::ptmaxcut 
               && TMath::Abs(nt->Gy[j]) >= fitX::ymincut && TMath::Abs(nt->Gy[j]) < fitX::ymaxcut
               )) continue;
          branchesf["Gpt"]->value = nt->Gpt[j];
          branchesf["Gy"]->value = nt->Gy[j];
          branchesf["Gmu1pt"]->value = nt->Gmu1pt[j];
          branchesf["Gmu2pt"]->value = nt->Gmu2pt[j];
          branchesf["Gmu1eta"]->value = nt->Gmu1eta[j];
          branchesf["Gmu2eta"]->value = nt->Gmu2eta[j];
          brancheso["isAcc"]->value = (
                                       (fabs(nt->Gmu1eta[j]) < 2.4) &&
                                       ((fabs(nt->Gmu1eta[j]) < 1.2 && nt->Gmu1pt[j] >= 3.5) ||
                                        (fabs(nt->Gmu1eta[j]) >= 1.2 && fabs(nt->Gmu1eta[j]) < 2.1 && nt->Gmu1pt[j] >= 5.47-1.89*fabs(nt->Gmu1eta[j])) ||
                                        (fabs(nt->Gmu1eta[j]) >= 2.1 && nt->Gmu1pt[j] >= 1.5)) &&
                                       (fabs(nt->Gmu2eta[j]) < 2.4) &&
                                       ((fabs(nt->Gmu2eta[j]) < 1.2 && nt->Gmu2pt[j] >= 3.5) ||
                                        (fabs(nt->Gmu2eta[j]) >= 1.2 && fabs(nt->Gmu2eta[j]) < 2.1 && nt->Gmu2pt[j] >= 5.47-1.89*fabs(nt->Gmu2eta[j])) ||
                                        (fabs(nt->Gmu2eta[j]) >= 2.1 && nt->Gmu2pt[j] >= 1.5)) &&
                                       nt->Gtk1pt[j] > 0.9 && nt->Gtk2pt[j] > 0.9 &&
                                       fabs(nt->Gtk1eta[j]) < 2.4 && fabs(nt->Gtk2eta[j]) < 2.4
                                       );
          outnt->Fill();
        }
    }
  xjjc::progressbar_summary(nentries);

  outnt->Write();
  outf->Close();
  std::cout<<"output file: "<<output<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==3) { acc_tree(argv[1], argv[2]); return 0; }
  if(argc==2) { acc_tree(argv[1]); return 0; }
  return 1;
}

