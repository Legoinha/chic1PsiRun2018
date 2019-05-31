#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "TMVA/Reader.h"

#include "TMVAClassification.h"
#include "xjjcuti.h"

#define MAX_XB       20000

namespace mytmva
{
  void createmva(TTree* nttree, TFile* outf, std::vector<std::string> methods, std::vector<std::string> xmlname, std::vector<int> stages);
}

const int nevt = -1;
void mvaprob(std::string inputname, std::string treename, std::string outputname, std::string outputfilename,
             float ptmin, float ptmax, std::string mymethod = "", std::string stage = "0,1,2,3,4,5,6,7,8,9,10")
{
  mymethod = xjjc::str_replaceall(mymethod, " ", "");
  std::vector<std::string> methods(xjjc::str_divide(mymethod, ","));
  stage = xjjc::str_replaceall(stage, " ", "");
  std::vector<int> stages;
  for(auto& ss : xjjc::str_divide(stage, ",")) { stages.push_back(atoi(ss.c_str())); }

  std::cout<<"\e[34;1m"<<"==> "<<__FUNCTION__<<": input file:"<<inputname<<"\e[0m"<<std::endl;
  std::string outfname(Form("%s_%s_%s_%s_%s.root", outputname.c_str(),xjjc::str_replaceallspecial(mymethod).c_str(),
                            xjjc::number_to_string(ptmin).c_str(), (ptmax<0?"inf":xjjc::number_to_string(ptmax).c_str()),
                            xjjc::str_replaceall(stage, ",", "-").c_str()));
  outfname = xjjc::str_replaceallspecial(outfname);

  std::string rootfname(Form("%s", xjjc::str_replaceall(xjjc::str_replaceall(outfname, "_root", ".root"), "rootfiles_", "rootfiles/").c_str()));
  bool findrootf = !gSystem->AccessPathName(rootfname.c_str());

  std::vector<std::string> xmlname, methodsclone(methods);
  for(int mm=0; mm<methodsclone.size(); mm++)
    { 
      std::string xml(Form("dataset/weights/%s/TMVAClassification_%s.weights.xml", outfname.c_str(), methodsclone[mm].c_str()));
      if(gSystem->AccessPathName(xml.c_str()))
        {
          std::cout<<"\e[31;1m"<<"==> "<<__FUNCTION__<<": skip "<<xml<<"\e[0m"<<std::endl;
          methods.erase(std::find(methods.begin(), methods.end(), methodsclone[mm]));
          continue;
        }
      xmlname.push_back(xml.c_str());
    }
  if(!methods.size()) return;
  for(auto& me : methods) { std::cout<<"\e[34;1m"<<"==> "<<__FUNCTION__<<": Registered method "<<me<<"\e[0m"<<std::endl; }

  std::string cuts = "", cutb = "", varinfo = "";
  if(findrootf)
    {
      TString *cuts_ = 0, *cutb_ = 0; std::string *varinfo_ = 0;
      TFile* rootf = TFile::Open(rootfname.c_str());
      std::cout<<"\e[34;1m"<<"==> "<<__FUNCTION__<<": Opening file:"<<rootfname<<"."<<"\e[0m"<<std::endl;
      if(!rootf) { std::cout<<"\e[31;1m"<<"==> "<<__FUNCTION__<<": error: file is not opened."<<"\e[0m"<<std::endl; return; }
      TTree* rinfo = (TTree*)rootf->Get("dataset/tmvainfo");
      if(!rinfo) { std::cout<<"\e[31;1m"<<"==> "<<__FUNCTION__<<": error: tree is not opened."<<"\e[0m"<<std::endl; return; }
      rinfo->SetBranchAddress("cuts", &cuts_);
      rinfo->SetBranchAddress("cutb", &cutb_);
      rinfo->SetBranchAddress("var", &varinfo_);
      std::cout<<"\e[34;1m"<<std::endl; rinfo->Show(0); std::cout<<"\e[0m"<<std::endl;
      rinfo->GetEntry(0);
      cuts = *cuts_; cutb = *cutb_; varinfo = *varinfo_;
      rootf->Close();
    }
  else { std::cout<<"\e[31;1m"<<"==> "<<__FUNCTION__<<": warning: file:"<<rootfname.c_str()<<" doesn't exist. skipped."<<"\e[0m"<<std::endl; }

  TFile* inf = TFile::Open(inputname.c_str());
  TTree* nttree = (TTree*)inf->Get(treename.c_str());
  outputfilename += ("_" + outfname + ".root");
  TFile* outf = TFile::Open(outputfilename.c_str(), "recreate");
  outf->mkdir("dataset");
  outf->cd("dataset");
  TTree* info = new TTree("tmvainfo", "TMVA info");
  info->Branch("cuts", &cuts);
  info->Branch("cutb", &cutb);
  info->Branch("var", &varinfo);
  info->Fill();

  outf->cd();
  mytmva::createmva(nttree, outf, methods, xmlname, stages);
  std::cout<<"\e[34;1m"<<"==> "<<__FUNCTION__<<": output file:"<<outputfilename<<"\e[0m"<<std::endl;
}

void mytmva::createmva(TTree* nttree, TFile* outf, std::vector<std::string> methods, std::vector<std::string> xmlname, std::vector<int> stages)
{
  //
  // mytmva::ntuple* nt = new mytmva::ntuple(nttree);
  mytmva::varval* values = new mytmva::varval(nttree);
  if(!values->isvalid()) { return; }

  //
  int nvar = stages.size();
  int nmeth = methods.size();

  std::string varnote("");
  std::vector<float> __varval(nvar, 0);
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 
  for(int vv=0; vv<nvar; vv++)
    { 
      reader->AddVariable(mytmva::varlist[stages[vv]].var.c_str(), &(__varval[vv])); 
      varnote += (";"+mytmva::varlist[stages[vv]].varname);
    }
  for(int mm=0; mm<nmeth; mm++)
    {
      std::string methodtag(methods[mm] + " method");
      reader->BookMVA( methodtag.c_str(), xmlname[mm].c_str() ); // ~
    }
  std::cout<<std::endl;

  outf->cd();
  // if no dataset, mkdir dataset !!
  outf->cd("dataset");
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
            { __varval[vv] = values->getval(mytmva::varlist[stages[vv]].varname, j); }
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

  outf->cd();
  outf->Write();
  outf->Close();
}

int main(int argc, char* argv[])  
{
  if(argc==9)
    {
      mvaprob(argv[1], argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]), argv[7], argv[8]);
      return 0;
    }
  return 1;
}
