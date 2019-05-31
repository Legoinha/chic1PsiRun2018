// https://root.cern.ch/doc/v610/TMVAGui_8cxx_source.html
// https://root.cern.ch/doc/v608/variables_8cxx_source.html
#include <iostream>
#include <string>

#include <TKey.h>
#include <TList.h>
#include <TObjString.h>
#include "TMVA/variables.h"

#include "xjjcuti.h"

namespace mytmva
{
  static TList* TMVAGui_keyContent;
  void guivariables(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage = "0,1,2,3");
  void variables(std::string outfname);
  TList* GetKeyList( const TString& pattern );
}

void mytmva::guivariables(std::string outputname, float ptmin, float ptmax, std::string mymethod, std::string stage/* = "0,1,2,3"*/)
{
  mymethod = xjjc::str_replaceall(mymethod, " ", "");
  std::string outfname(Form("%s_%s_%s_%s_%s.root", outputname.c_str(), xjjc::str_replaceallspecial(mymethod).c_str(), 
                            xjjc::number_to_string(ptmin).c_str(), (ptmax<0?"inf":xjjc::number_to_string(ptmax).c_str()), 
                            xjjc::str_replaceall(stage, ",", "-").c_str()));
  mytmva::variables(outfname);
}

void mytmva::variables(std::string outfname)
{
  TString dataset("");
  std::string outputstr = xjjc::str_replaceallspecial(outfname);
  // Set up dataset
  TFile* file = TFile::Open( outfname.c_str() );
  if(!file)
    { std::cout << "==> Abort " << __FUNCTION__ << ", please verify filename." << std::endl; return; }
  if(file->GetListOfKeys()->GetEntries()<=0)
    { std::cout << "==> Abort " << __FUNCTION__ << ", please verify if dataset exist." << std::endl; return; }
  // --->>
  if( (dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()==1))
    { dataset = ((TKey*)file->GetListOfKeys()->At(0))->GetName(); }
  // <<---
  else if((dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()>=1))
    {
      std::cout << "==> Warning " << __FUNCTION__ << ", more than 1 dataset." << std::endl;
      for(int i=0;i<file->GetListOfKeys()->GetEntries();i++)
        {
          dataset = ((TKey*)file->GetListOfKeys()->At(i))->GetName();
          std::cout << "    " << dataset << std::endl;
        }
      return;
    }
  else { return; }

  TMVAGui_keyContent = (TList*)file->GetDirectory(dataset.Data())->GetListOfKeys()->Clone(); //
  
  TList* keylist = GetKeyList("InputVariables");
  TListIter it( keylist );
  TObjString* str = 0;
  while((str = (TObjString*)it()))
    {
      TString tmp   = str->GetString();
      TString title = Form( "Input variables '%s'-transformed (training sample)", 
                            tmp.ReplaceAll("InputVariables_","").Data() );
      if(tmp.Contains( "Id" )) title = "Input variables (training sample)";
      // void TMVA::variables (TString dataset, TString fin="TMVA.root", TString dirName="InputVariables_Id", TString title="TMVA Input Variables", Bool_t isRegression=kFALSE, Bool_t useTMVAStyle=kTRUE)
      TMVA::variables(dataset.Data(), outfname.c_str(), str->GetString().Data(), title.Data());
    }

  gSystem->Exec(Form("rm %s/plots/*.png", dataset.Data()));
  gSystem->Exec(Form("mkdir -p %s/plots/%s", dataset.Data(), outputstr.c_str()));
  gSystem->Exec(Form("mv %s/plots/*.eps %s/plots/%s/", dataset.Data(), dataset.Data(), outputstr.c_str()));
}

TList* mytmva::GetKeyList( const TString& pattern )
{
  TList* list = new TList();
  TIter next( TMVAGui_keyContent );
  TKey* key(0);
  while((key = (TKey*)next())) {         
    if (TString(key->GetName()).Contains( pattern )) { list->Add( new TObjString( key->GetName() ) ); }
  }
  return list;
}

int main(int argc, char* argv[])
{
  if(argc==6)
    { mytmva::guivariables(argv[1], atof(argv[2]), atof(argv[3]), argv[4], argv[5]); return 0; }
  return 1;
}
