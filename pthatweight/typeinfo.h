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
