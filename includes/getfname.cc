#include "fitX.h"
#include <iostream>

int main(int argc, char* argv[])
{
  if(argc==7)
    { fitX::init(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]));
      std::cout<<fitX::tagname()<<std::endl;
      return 0; }
  return 1;
}
