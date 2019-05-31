#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TClass.h>
#include <TObject.h>
#include <TROOT.h>

#include <iostream>
#include <map>
#include <string>

namespace xjjc
{
  class copydir_recur
  {
  public:
    std::map<TTree*, TDirectory*> trees_output;
    std::map<TTree*, TDirectory*> trees_input;
    copydir_recur(TDirectory* source, TDirectory* des, int fnevt=-1);
    int getentries() { return fnentries; }

  private:
    int nevt; // default = -1
    void copydir(TDirectory* source);
    int fll;
    std::vector<bool> lastll;
    bool lastiter(const TIter& current); 
    int fnentries;
  };
}

xjjc::copydir_recur::copydir_recur(TDirectory* src, TDirectory* des, int fnevt/*=-1*/) : nevt(fnevt), fll(0), fnentries(0)
{
  trees_output.clear();
  trees_input.clear();
  lastll.clear();
  lastll.push_back(true);
  des->cd();
  std::cout<<std::endl<<"   \e[1m./\e[0m"<<std::endl;
  copydir(src);
  std::cout<<std::endl;
  des->cd();
}

//
void xjjc::copydir_recur::copydir(TDirectory* source)
{
  // source->ls();
  TDirectory* savdir = gDirectory;
  TDirectory* adir = savdir;
  if(source->GetMotherDir()) adir = savdir->mkdir(source->GetName());
  adir->cd();
  TKey* key;
  TIter nextkey(source->GetListOfKeys());
  lastll.push_back(true); fll = lastll.size()-1;
  while((key = (TKey*)nextkey()))
    {
      const std::string classname(key->GetClassName()), keyname(key->GetName());
      TClass* cl = gROOT->GetClass(classname.c_str());

      // >> print
      lastll[fll] = lastiter(nextkey);
      for(auto lll = lastll.begin();lll<lastll.end()-1;lll++) { std::cout<<(*lll?" ":"\u2502")<<"  "; }
      std::cout<<(lastll[fll]?"\u2514":"\u251C")<<"\u2500 \e[2m("<<classname<<") =>\e[0m \e[1m"<<keyname<<"\e[0m";
      if(cl->InheritsFrom(TTree::Class())) { std::cout<<" \e[2m("<<((TTree*)source->Get(keyname.c_str()))->GetEntries()<<")\e[0m"; }
      std::cout<<std::endl;
      // << print

      if(!cl) { std::cout<<"\e[31;1m warning: unknown class "<<classname<<"\e[0m"<<std::endl; continue; }
      if(cl->InheritsFrom(TDirectory::Class()))
        {
          source->cd(key->GetName());
          TDirectory* subdir = gDirectory;
          adir->cd();
          copydir(subdir);
          adir->cd();
        }
      else if(cl->InheritsFrom(TTree::Class()))
        {
          TTree* nt = (TTree*)source->Get(key->GetName());
          trees_input.insert(std::pair(nt, source));
          if(nt->GetEntries()<fnentries || fnentries==0) { fnentries = nt->GetEntries(); }
          adir->cd();
          TTree* nt_new = nt->CloneTree(nevt); // nevt
          trees_output.insert(std::pair(nt_new, adir));
        }
      else
        {
          source->cd();
          TObject* obj = key->ReadObj();
          adir->cd();
          obj->Write();
          delete obj;
        }
    }
  lastll.pop_back(); fll = lastll.size()-1;
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

bool xjjc::copydir_recur::lastiter(const TIter& current)
{
  TIter iterclone(current);
  return iterclone()==0;
}
