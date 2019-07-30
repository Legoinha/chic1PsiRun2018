#include <map>
#include <vector>
#include <string>

namespace tnpcc
{
  const int filterId = 1; // * filterId = 1: Jpsi L3 filter
  std::map<int, std::string> idxname = {
    std::pair<int, std::string>(-1, "syst_u"),
    std::pair<int, std::string>(-2, "syst_d"),
    std::pair<int, std::string>(+1, "stat_u"),
    std::pair<int, std::string>(+2, "stat_d"),
    std::pair<int, std::string>(0, "nominal"),
  };
  const std::vector<std::string> types = {"nominal", "trg", "trk", "muid", "total"};
  std::map<std::string, Color_t> typecolor = {
    std::pair<std::string, Color_t>("nominal", kBlack),
    std::pair<std::string, Color_t>("trg", xjjroot::mycolor_middle["azure"]),
    std::pair<std::string, Color_t>("trk", xjjroot::mycolor_middle["red"]),
    std::pair<std::string, Color_t>("muid", xjjroot::mycolor_middle["green"]),
    std::pair<std::string, Color_t>("total", xjjroot::mycolor_middle["yellow"]),
  };
  std::map<std::string, Style_t> idxstyle = {
    std::pair<std::string, Style_t>("syst_u", 2),
    std::pair<std::string, Style_t>("syst_d", 2),
    std::pair<std::string, Style_t>("stat_u", 4),
    std::pair<std::string, Style_t>("stat_d", 4),
    std::pair<std::string, Style_t>("nominal", 1),
  };
  std::vector<std::string> err({"stat", "syst"});

  __PTBIN_INPUT__
  int nptbins = sizeof(ptbins)/sizeof(ptbins[0]) - 1;
  float scalemin = 0.9, scalemax = 1.4;
}
