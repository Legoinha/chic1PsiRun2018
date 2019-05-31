#include "copydir_recur.h"

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#define MAX_XB       20000

namespace xjjc
{
  void progressbar(int index_, int total_, int morespace_=0);
  void progressbar_summary(int total_);
  template<typename T> std::string number_to_string(T param_);
  std::string str_replaceall(std::string strs, std::string sub, std::string newsub);
  bool str_contains(std::string str1, std::string str2) { return str1.find(str2)!=std::string::npos; }
}

void pickpart(std::string inputname, float compress, std::string outputname="")
{
  if(compress > 1) { std::cout<<__FUNCTION__<<": error: [compress] should be smaller than 1."<<std::endl; return; }
  if(outputname=="")
    { outputname = xjjc::str_replaceall(inputname, ".root", "_"+xjjc::number_to_string(compress*100)+"per.root"); }
  else if(!xjjc::str_contains(outputname, ".root") || inputname==outputname)
    { std::cout<<__FUNCTION__<<": error: invalid outputname."<<std::endl<<"  "<<outputname<<std::endl; return; }
  std::cout<<outputname<<std::endl;

  TFile* inf = new TFile(inputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  // >> 
  xjjc::copydir_recur dirs(inf, outf, 0);
  // <<
  int interval = (int)(1/compress);

  // **************************************************
  TTree* ntmix = (TTree*)inf->Get("Bfinder/ntmix");
  int Bsize; ntmix->SetBranchAddress("Bsize", &Bsize);
  float Bpt[MAX_XB]; ntmix->SetBranchAddress("Bpt", Bpt);
  // **************************************************

  int nentries = dirs.getentries();
  std::cout<<" -- Event reading"<<std::endl;
  for(int i=0;i<nentries;i+=interval)
    {
      if((i/interval)%100==0) { xjjc::progressbar(i, nentries); }

      for(auto& inp : dirs.trees_input) { inp.first->GetEntry(i); }

      // --> selections
      int pass = 0;
      for(int j=0; j<Bsize; j++) { if(Bpt[j]>10) { pass++; break; } }
      if(!pass) { continue; }
      // <-- selections

      for(auto& outp : dirs.trees_output)
        { 
          outp.second->cd();
          outp.first->Fill(); 
        }
    }
  xjjc::progressbar_summary(nentries);
      
  outf->Write();
  std::cout<<" -- Writing new trees done"<<std::endl;
  outf->Close();
      
}

/* ==================== */
void xjjc::progressbar(int index_, int total_, int morespace_/*=0*/)
{
  std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10+morespace_)<<index_<<"\033[0m"<<" / "<<std::setw(10+morespace_)<<total_<<" ] "<<"\033[1;36m"<<(int)(100.*index_/total_)<<"%\033[0m"<<"\r"<<std::flush;
}
void xjjc::progressbar_summary(int total_)
{
  std::cout<<std::endl<<"  Processed "<<"\033[1;31m"<<total_<<"\033[0m event(s)."<<std::endl;
}
template<typename T>
std::string xjjc::number_to_string(T param_)
{
  if(param_<0) return "";
  std::string str = std::to_string(param_);
  std::size_t found = str.find('.');
  if(found==std::string::npos)
    {
      str+="p0";
      return str;
    }
  str.replace(found,1,"p");
  while(*(str.end()-1)=='0' && *(str.end()-2)!='p' && !str.empty()) str.erase(str.end()-1);
  return str;
}
std::string xjjc::str_replaceall(std::string strs, std::string sub, std::string newsub)
{
  std::string result(strs);
  std::string str(strs);
  size_t pos = str.find(sub, 0);
  size_t thispos = 0;
  while(pos != std::string::npos)
    {
      result.replace(pos + thispos, sub.size(), newsub);
      thispos += (newsub.size() - sub.size());
      pos = str.find(sub, pos+1);
    }
  return result;
}

int main(int argc, char* argv[])
{
  if(argc==4) { pickpart(argv[1], atof(argv[2]), argv[3]); return 0; }
  if(argc==3) { pickpart(argv[1], atof(argv[2])); return 0; }
  std::cout<<__FUNCTION__<<": ./pickpart.exe [input] [compress frac] (optional)[output]"<<std::endl;
  return 1;
}


