#ifndef __MVAPROD_H_
#define __MVAPROD_H_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <experimental/filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "TMVAClassification.h"
#include "xjjcuti.h"

#ifndef MAX_XB
#define MAX_XB       20000
#endif

namespace fs = std::experimental::filesystem;
namespace mytmva
{
  void createmva(TTree* nttree, TFile* outf, std::vector<std::string> xmlname, int nevt=-1);
  std::string titlecolor = "\e[34;3m", nocolor = "\e[0m", contentcolor = "\e[34m", errorcolor = "\e[31;1m";
  void mvaprob(std::string inputname, std::string treename, std::string outputname, std::string weightdir,
               int nevt=-1, std::string rootfname="");
}

void mytmva::mvaprob(std::string inputname, std::string treename, std::string outputname, std::string weightdir, int nevt, std::string rootfname)
{
  std::cout<<std::endl;
  if(weightdir.back() == '/') { weightdir.pop_back(); }
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": directory of weight files:"<<mytmva::nocolor<<std::endl;
  std::cout<<weightdir<<std::endl;

  // resolve methods
  std::vector<std::string> xmlname;
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": found weight files:"<<mytmva::nocolor<<std::endl;
  for (const auto & entry : fs::directory_iterator(weightdir))
    {
      std::string entrypath(entry.path());
      if(xjjc::str_contains(entrypath, ".weights.xml")) 
        {
          xmlname.push_back(entrypath);
          std::cout<<entrypath<<std::endl;
        }
    }

  // input/output file  
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": input file:"<<mytmva::nocolor<<std::endl<<inputname<<mytmva::nocolor<<std::endl;
  std::string weightlabel = xjjc::str_replaceall(xjjc::str_replaceallspecial(weightdir), "dataset_weights_rootfiles_TMVA_", "");
  std::string outfname(Form("%s_%s.root", outputname.c_str(), weightlabel.c_str())); outfname = xjjc::str_replaceall(outfname, "_root.root", ".root");
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": output file:"<<mytmva::nocolor<<std::endl<<outfname<<mytmva::nocolor<<std::endl;
  if(std::experimental::filesystem::exists(outfname)) { std::cout<<mytmva::errorcolor<<"==> "<<__FUNCTION__<<": warning: output file already exists."<<mytmva::nocolor<<std::endl; }
  std::cout<<"==> "<<__FUNCTION__<<": warning: application of mva values will take long time. would you want to continue? [y/n]"<<std::endl; char ans='x';
  while(ans!='y' && ans!='n') { std::cin>>ans; }
  if(ans=='n') return;
  gSystem->Exec(Form("rsync --progress %s %s", inputname.c_str(), outfname.c_str()));

  // training rootfile
  if(rootfname == "")
    {
      std::string reviserootf = xjjc::str_replaceall(weightdir, "dataset/weights/rootfiles_", "rootfiles/");
      reviserootf = xjjc::str_replaceall(reviserootf, "_root", ".root");
      rootfname = reviserootf;
    }
  bool findrootf = !gSystem->AccessPathName(rootfname.c_str());
  std::string cuts = "", cutb = "", varinfo = "";
  if(findrootf)
    {
      TString *cuts_ = 0, *cutb_ = 0; std::string *varinfo_ = 0;
      TFile* rootf = TFile::Open(rootfname.c_str());
      std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": opening file:"<<mytmva::nocolor<<std::endl<<rootfname<<mytmva::nocolor<<std::endl;
      if(!rootf) { std::cout<<mytmva::errorcolor<<"==> "<<__FUNCTION__<<": error: file is not opened."<<mytmva::nocolor<<std::endl; return; }
      TTree* rinfo = (TTree*)rootf->Get("dataset/tmvainfo");
      if(!rinfo) { std::cout<<mytmva::errorcolor<<"==> "<<__FUNCTION__<<": error: tree is not opened."<<mytmva::nocolor<<std::endl; return; }
      rinfo->SetBranchAddress("cuts", &cuts_);
      rinfo->SetBranchAddress("cutb", &cutb_);
      rinfo->SetBranchAddress("var", &varinfo_);
      // std::cout<<mytmva::titlecolor<<std::endl; rinfo->Show(0); std::cout<<mytmva::nocolor<<std::endl;
      std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": mva info:"<<mytmva::nocolor<<std::endl;
      rinfo->Show(0); std::cout<<std::endl;
      rinfo->GetEntry(0);
      cuts = *cuts_; cutb = *cutb_; varinfo = *varinfo_;
      rootf->Close();
    }
  else { std::cout<<"\e[33m"<<"==> "<<__FUNCTION__<<": warning: file:"<<rootfname.c_str()<<" doesn't exist. skipped."<<mytmva::nocolor<<std::endl; }

  // fill cut info
  TFile* inf = TFile::Open(inputname.c_str());
  TTree* nttree = (TTree*)inf->Get(treename.c_str());
  TFile* outf = TFile::Open(outfname.c_str(), "update");
  outf->mkdir("dataset");
  outf->cd("dataset");
  TTree* info = new TTree("tmvainfo", "TMVA info");
  info->Branch("cuts", &cuts);
  info->Branch("cutb", &cutb);
  info->Branch("var", &varinfo);
  info->Fill();
  info->Write("", TObject::kOverwrite);

  outf->cd();
  mytmva::createmva(nttree, outf, xmlname, nevt);
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": output file:"<<mytmva::nocolor<<std::endl<<outfname<<mytmva::nocolor<<std::endl;
}

void mytmva::createmva(TTree* nttree, TFile* outf, std::vector<std::string> xmlname, int nevt)
{
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Set branch address:"<<mytmva::nocolor<<std::endl;
  mytmva::varval* values = new mytmva::varval(nttree);
  if(!values->isvalid()) { return; }

  // read weight file
  std::vector<std::string> methods;
  std::map<std::string, std::vector<std::string>> varlabels;
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Found method:"<<mytmva::nocolor<<std::endl;
  for(auto& ww : xmlname)
    {
      const char* filename = ww.c_str();
      void *doc = TMVA::gTools().xmlengine().ParseFile(filename,TMVA::gTools().xmlenginebuffersize());
      void* rootnode = TMVA::gTools().xmlengine().DocGetRootElement(doc); //node "MethodSetup"
      // method
      std::string fullmethodname("");
      TMVA::gTools().ReadAttr(rootnode, "Method", fullmethodname);
      std::string method = fullmethodname;
      std::size_t found = fullmethodname.find("::");
      method.erase(0, found+2);
      methods.push_back(method);
      std::cout<<std::left<<std::setw(10)<<method<<" // "<<fullmethodname<<mytmva::nocolor<<std::endl;
      // variable
      void* variables = TMVA::gTools().GetChild(rootnode, "Variables");
      UInt_t NVar=0;
      TMVA::gTools().ReadAttr(variables, "NVar", NVar);
      void* var = TMVA::gTools().GetChild(variables, "Variable");
      for(unsigned int k=0;k<NVar;k++)
        {
          std::string varlabel("");
          TMVA::gTools().ReadAttr(var, "Label", varlabel);
          varlabels[method].push_back(varlabel);
          var = TMVA::gTools().GetNextChild(var);
        }
    }
  if(methods.size() <= 0) { std::cout<<__FUNCTION__<<": error: no method is registered."<<std::endl; return; }

  // 
  int nmeth = methods.size();
  std::vector<std::string> varnames = varlabels[methods[0]];
  int nvar = varnames.size();
  for(auto& mm : methods)
    {
      if(varlabels[mm].size() != nvar)
        { std::cout<<__FUNCTION__<<": error: inconsistent variable number among methods."<<std::endl; return; }
      for(int vv=0; vv<nvar; vv++)
        {
          if(varlabels[mm].at(vv) != varnames[vv]) 
            { std::cout<<__FUNCTION__<<": error: inconsistent variable among methods."<<std::endl; return; }
        }
    }
  //
  std::string varnote("");
  std::vector<float> __varval(nvar, 0);
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Add variable:"<<mytmva::nocolor<<std::endl;
  for(int vv=0; vv<nvar; vv++)
    {
      std::cout<<std::left<<std::setw(10)<<varnames[vv]<<" // "<<mytmva::findvar(varnames[vv])->var.c_str()<<std::endl;
      reader->AddVariable(mytmva::findvar(varnames[vv])->var.c_str(), &(__varval[vv])); 
      varnote += (";"+varnames[vv]);
    }
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Book method:"<<mytmva::nocolor<<std::endl;
  for(int mm=0; mm<nmeth; mm++)
    {
      std::string methodtag(methods[mm] + " method");
      reader->BookMVA( methodtag.c_str(), xmlname[mm].c_str() ); // ~
    }

  outf->cd();
  if(!outf->cd("dataset")) { outf->mkdir("dataset"); outf->cd("dataset"); }
  TTree* note = new TTree("weightnote", "");
  note->Branch("varnote", &varnote);
  note->Fill();

  int mvaBsize;
  std::vector<float[MAX_XB]> __mvaval(nmeth);
  outf->cd("dataset");
  TTree* mvatree = new TTree("mva", "");
  mvatree->Branch("mvaBsize", &mvaBsize);
  for(int mm=0; mm<nmeth; mm++)
    { mvatree->Branch(methods[mm].c_str(), __mvaval[mm], Form("%s[mvaBsize]/F", methods[mm].c_str())); }
  bool __mvapref[MAX_XB];
  mvatree->Branch("mvapref", __mvapref, "mvapref[mvaBsize]/O");
  
  std::cout<<mytmva::titlecolor<<"==> "<<__FUNCTION__<<": Filling mva values:"<<mytmva::nocolor<<std::endl;
  outf->cd();
  int nentries = nevt>0&&nevt<values->getnt()->getnt()->GetEntries()?nevt:values->getnt()->getnt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      if(i%100==0) xjjc::progressbar(i, nentries);
      values->getnt()->getnt()->GetEntry(i);

      mvaBsize = values->getnt()->Bsize;
      for(int j=0; j<values->getnt()->Bsize; j++)
        {
          for(int vv=0; vv<nvar; vv++)
            { __varval[vv] = values->getval(varnames[vv], j); }
          for(int mm=0; mm<nmeth; mm++)
            {
              std::string methodtag(methods[mm] + " method");
              __mvaval.at(mm)[j] = reader->EvaluateMVA(methodtag.c_str());
            }
          __mvapref[j] = values->getnt()->passedpre(j);
        }
      outf->cd("dataset"); 
      mvatree->Fill(); 
    }
  xjjc::progressbar_summary(nentries);

  outf->cd("dataset");
  mvatree->Write("", TObject::kOverwrite);
  // outf->Write();
  outf->cd();
  outf->Close();
}

#endif
