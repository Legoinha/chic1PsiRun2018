#include <string>

#include "TMVAClassification.h"
#include "mvaprod.h"

void mvaprob_main(std::string inputname, std::string treename, std::string outputname, std::string outputfilename,
                  float ptmin, float ptmax, std::string mymethod = "", std::string stage = "0,1,2,3,4,5,6,7,8,9,10")
{
  std::string outfname = xjjc::str_replaceallspecial(mytmva::mkname(outputname, ptmin, ptmax, mymethod, stage));
  mytmva::mvaprob(inputname, "Bfinder/ntmix", outputfilename, Form("dataset/weights/%s", outfname.c_str()));
}

int main(int argc, char* argv[])  
{
  if(argc==9)
    {
      mvaprob_main(argv[1], argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]), argv[7], argv[8]);
      return 0;
    }
  return 1;
}
