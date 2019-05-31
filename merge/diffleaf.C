#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TLeaf.h>
#include <TBranch.h>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

void diffleaf(std::string input_a, std::string input_b, std::string treename)
{
  TFile* inf_a = new TFile(input_a.c_str());
  TTree* nt_a = (TTree*)inf_a->Get(treename.c_str());
  TFile* inf_b = new TFile(input_b.c_str());
  TTree* nt_b = (TTree*)inf_b->Get(treename.c_str());

  auto branchlist_a = nt_a->GetListOfBranches(); // TObjArray*
  const int nbr_a = nt_a->GetNbranches();
  std::vector<std::string> brname_a;
  for(int i=0; i<nbr_a; i++) { brname_a.push_back(branchlist_a->At(i)->GetName()); 
    // std::cout<<branchlist_a->At(i)->GetName()<<std::endl; 
  }

  std::cout<<std::endl<<std::string(25, '#')<<std::endl;

  auto branchlist_b = nt_b->GetListOfBranches(); // TObjArray*
  const int nbr_b = nt_b->GetNbranches();
  std::vector<std::string> brname_b;
  for(int i=0; i<nbr_b; i++) { brname_b.push_back(branchlist_b->At(i)->GetName()); 
    // std::cout<<branchlist_b->At(i)->GetName()<<std::endl; 
  }

  std::vector<std::string> only_a, only_b, both_ab, brname_b_copy(brname_b);
  for(auto& ia : brname_a) 
    {
      auto pos = std::find(brname_b_copy.begin(), brname_b_copy.end(), ia);
      if(pos == brname_b_copy.end())
        {
          only_a.push_back(ia);
          brname_b_copy.erase(pos);
        }
      else { both_ab.push_back(ia); }
    }
  for(auto& ib : brname_b_copy)
    {
      auto pos = std::find(brname_a.begin(), brname_a.end(), ib);
      if(pos == brname_a.end()) { only_b.push_back(ib); }
    }

  std::cout<<"\e[41;97;1m"<<input_a<<"\e[0m"<<std::endl<<"\e[41;97;1monly_a >>>\e[31;40;1m"<<std::endl;
  for(auto& o : only_a) { std::cout<<o<<std::endl; }
  std::cout<<std::endl;
  std::cout<<"\e[42;97;1m"<<input_b<<"\e[0m"<<std::endl<<"\e[42;97;1monly_b >>>\e[32;40;1m"<<std::endl;
  for(auto& o : only_b) { std::cout<<o<<std::endl; }
  std::cout<<"\e[0m";

  //
  // for(auto& o : only_a)
  //   {
  //     TBranch* bro = nt_a->GetBranch(o.c_str());
  //     TLeaf* leo = nt_a->GetLeaf(o.c_str());
  //     nt_a->GetListOfBranches()->Remove(bro);
  //     nt_a->GetListOfLeaves()->Remove(leo);
  //     nt_a->Write();
  //   }
  // for(auto& o : only_b)
  //   {
  //     TBranch* bro = nt_b->GetBranch(o.c_str());
  //     TLeaf* leo = nt_b->GetLeaf(o.c_str());
  //     nt_b->GetListOfBranches()->Remove(bro);
  //     nt_b->GetListOfLeaves()->Remove(leo);
  //     nt_b->Write();
  //   }
}

int main(int argc, char* argv[])
{
  if(argc==4) { diffleaf(argv[1], argv[2], argv[3]); return 0; }
  std::cout<<__FUNCTION__<<": ./diffleaf [file1] [file2] [treename]"<<std::endl;
  return 1;
}
