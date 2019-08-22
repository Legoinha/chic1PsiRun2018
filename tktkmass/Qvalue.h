#include <vector>
#include <string>
#include "xjjrootuti.h"

namespace Qvalue
{
  int NBIN = 40;
  float BIN_MIN = 0, BIN_MAX = 0.2;
  std::string Q = "Bmass-3.096916-Btktkmass";
  std::vector<std::string> types = {"a", "brho", "bpipi"};
  std::vector<std::string> legtrs = {"#psi' (#pi#pi)", "X (#rho resonance)", "X (#pi#pi)"};
  std::vector<int> colors = {xjjroot::mycolor_middle["red"], xjjroot::mycolor_middle["green"], xjjroot::mycolor_middle["azure"]};
  std::vector<Style_t> line = {7, 7, 7};
    
}

namespace ptimb
{
  int NBIN = 40;
  float BIN_MIN = 0, BIN_MAX = 1.;
}
