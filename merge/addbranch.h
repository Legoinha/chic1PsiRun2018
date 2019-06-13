#include <iostream>
#include <iomanip>
#include <cxxabi.h>
#include <string>
#include <map>

#include <RtypesCore.h>

namespace xjjc
{
  template<typename T>
  class typeinfo
  {
  public:
    typeinfo(std::string n__) : name(n__) { status = 0; division = std::string(61, '-'); isroottype = roottype.find(name)!=roottype.end(); }
    T t;
    std::string name;
    char* gettype();
    void printinfo(bool title=false, bool end=false);
  private:
    int status;
    std::string division;
    bool isroottype;
    void printtitle();
    void printend() { std::cout << std::left << division << std::endl; }
    std::map<std::string, char> roottype = {
      std::pair<std::string, char>("char", 'B'),
      std::pair<std::string, char>("Char_t", 'B'),
      std::pair<std::string, char>("UChar_t", 'b'),
      std::pair<std::string, char>("Short_t", 'S'),
      std::pair<std::string, char>("UShort_t", 's'),
      std::pair<std::string, char>("int", 'I'),
      std::pair<std::string, char>("Int_t", 'I'),
      std::pair<std::string, char>("UInt_t", 'i'),
      std::pair<std::string, char>("float", 'F'),
      std::pair<std::string, char>("Float_t", 'F'),
      std::pair<std::string, char>("Float16_t", 'f'),
      std::pair<std::string, char>("double", 'D'),
      std::pair<std::string, char>("Double_t", 'D'),
      std::pair<std::string, char>("Double32_t", 'd'),
      std::pair<std::string, char>("Long64_t", 'L'),
      std::pair<std::string, char>("ULong64_t", 'l'),
      std::pair<std::string, char>("bool", 'O'),
      std::pair<std::string, char>("Bool_t", 'O'),
      std::pair<std::string, char>("TString", 'C'),
    };
  };
  void addbranch(std::string inputname, std::string outputname, std::string treename);
  bool str_isnumber(std::string strs) { return (std::regex_match(strs, std::regex("-?[0-9]+([.][0-9]*)?")) || std::regex_match(strs, std::regex("-?[0-9]*[.][0-9]+"))); }
  bool str_isinteger(std::string strs) { return std::regex_match(strs, std::regex("-?[0-9]+")); }
  std::vector<std::string> str_divide(std::string str, std::string div);
  std::string str_replaceall(std::string strs, std::string sub, std::string newsub);
  std::vector<std::string> speciallist = {" ", "/", "(", ")", "^", "#", "%", "$", ",", ".", "*"};
  std::string str_replaceallspecial(std::string strs);
  void progressbar(int index_, int total_, int morespace_=0);
  void progressbar_summary(int total_);
}

template<typename T>
char* xjjc::typeinfo<T>::gettype()
{
  const std::type_info &ti = typeid(t);
  char* humantypename = abi::__cxa_demangle(ti.name(), 0, 0, &status);
  return humantypename;
}

template<typename T>
void xjjc::typeinfo<T>::printinfo(bool title, bool end)
{
  if(title) printtitle();
  std::cout << std::left << division << std::endl << "| " 
            << std::setw(20) << name 
            << std::setw(25) << gettype() 
            << std::setw(8) << sizeof(t) 
            << std::setw(5) << (isroottype?roottype[name]:'-')
            << "|" <<std::endl;
  if(end) printend();
}

template<typename T>
void xjjc::typeinfo<T>::printtitle()
{ 
  std::cout << std::left << division << std::endl << "| " 
            << std::setw(20) << "Type" 
            << std::setw(25) << "xjjc::gettype" 
            << std::setw(8) << "sizeof" 
            << std::setw(5) << "ROOT"
            << "|" <<std::endl; 
}

std::vector<std::string> xjjc::str_divide(std::string str, std::string div)
{
  std::vector<std::string> token;
  size_t lastpos = 0;
  size_t pos = str.find(div, 0);
  while(pos != std::string::npos)
    {
      std::string thistoken = str.substr(lastpos, pos-lastpos);
      token.push_back(thistoken);
      lastpos = pos+1;
      pos = str.find(div, pos+1);
    }
  if(lastpos != 0) { token.push_back(str.substr(lastpos, str.size()-lastpos)); }
  return token;
}

void xjjc::progressbar(int index_, int total_, int morespace_/*=0*/)
{
  std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10+morespace_)<<index_<<"\033[0m"<<" / "<<std::setw(10+morespace_)<<total_<<" ] "<<"\033[1;36m"<<(int)(100.*index_/total_)<<"%\033[0m"<<"\r"<<std::flush;
}

void xjjc::progressbar_summary(int total_)
{
  std::cout<<std::endl<<"  Processed "<<"\033[1;31m"<<total_<<"\033[0m event(s)."<<std::endl;
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

std::string xjjc::str_replaceallspecial(std::string strs)
{
  std::string result(strs);
  for(auto& isp : xjjc::speciallist) { result = xjjc::str_replaceall(result, isp, "-"); }
  return result;
}
