#include <vector>
#include <map>
#include "xjjrootuti.h"
#include "xjjcuti.h"

#include "fitX.h"

namespace variation
{
  std::vector<float> bdtgcut = {0, 0.2, 0.4, 0.6, 0.7, 0.8, 1.0};
  std::vector<float> bdtg = {0., 0.6, 0.8, 1.0};
  std::vector<float> chi2cl = {0.1, 0.5, 0.8, 1.0};
  std::vector<float> trk2pt = {0.9, 1.2, 2.0, 5.0};
  std::vector<float> trkptimba = {0, 0.2, 0.4, 0.6};
  std::map<std::string, std::vector<float>> varvectors =
    {
      std::pair<std::string, std::vector<float>>("BDTGcut", bdtgcut),
      std::pair<std::string, std::vector<float>>("BDTG", bdtg),
      std::pair<std::string, std::vector<float>>("#chi^{2} probability", chi2cl),
      std::pair<std::string, std::vector<float>>("|p_{T}^{trk1}-p_{T}^{trk2}|/(p_{T}^{trk1}+p_{T}^{trk2})", trkptimba),
      std::pair<std::string, std::vector<float>>("p_{T}^{trk}", trk2pt),
    };
  float massmin_a = 3.62, massmax_a = 3.80, nmass_a = 18;
  float massmin_b = 3.80, massmax_b = 3.98, nmass_b = 18;

  void drawalltext()
  {
    xjjroot::drawCMSleft();
    xjjroot::drawCMSright();
    xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
    xjjroot::drawtex(0.24, 0.84-0.042, "X(3872)", 0.038, 12, 62, fitX::color_b);
    xjjroot::drawtex(0.90, 0.84, Form("p_{T} > %.0f GeV/c", fitX::ptcut), 0.038, 32, 62);
    xjjroot::drawtex(0.90, 0.84-0.042, Form("|y| < %.1f", fitX::ycut), 0.038, 32, 62);
  }
  void drawalltext_simulation()
  {
    xjjroot::drawCMSleft("Simulation");
    xjjroot::drawCMSright();
    xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
    xjjroot::drawtex(0.24, 0.84-0.042, "X(3872)", 0.038, 12, 62, fitX::color_b);
    xjjroot::drawtex(0.90, 0.84, Form("p_{T} > %.0f GeV/c", fitX::ptcut), 0.038, 32, 62);
    xjjroot::drawtex(0.90, 0.84-0.042, Form("|y| < %.1f", fitX::ycut), 0.038, 32, 62);
  }
  void drawtext(std::string text)
  {
    xjjroot::drawCMSleft();
    xjjroot::drawCMSright();
    xjjroot::drawtex(0.24, 0.84, text.c_str(), 0.038, 12, 62, kBlack);
    xjjroot::drawtex(0.90, 0.84, Form("p_{T} > %.0f GeV/c", fitX::ptcut), 0.038, 32, 62);
    xjjroot::drawtex(0.90, 0.84-0.042, Form("|y| < %.1f", fitX::ycut), 0.038, 32, 62);
  }

  void drawtexlist(int iscutordis, std::vector<float>& vvector, int nv, std::string vartitle, float xoffset=0, float yoffset=0)
  {
    for(int l=0; l<nv; l++)
      {
        if(iscutordis) { xjjroot::drawtex(0.70+xoffset+0.08*(l%3), 0.85+yoffset-0.04*(l/3), Form("%s%s", xjjc::number_remove_zero(vvector[l]).c_str(), (l%3==(3-1)?"":", ")), 0.035, 12, 62, xjjroot::colorlist_dark[l]); }
        else { xjjroot::drawtex(0.70+xoffset+0.12*(l%2), 0.85+yoffset-0.04*(l/2), Form("(%s,%s)%s", xjjc::number_remove_zero(vvector[l]).c_str(), xjjc::number_remove_zero(vvector[l+1]).c_str(), (l%2==(2-1)?"":", ")), 0.035, 12, 62, xjjroot::colorlist_dark[l]); }
      }
    xjjroot::drawtex(0.70+xoffset-0.05, 0.85+yoffset, Form("%s %s",vartitle.c_str(),iscutordis?">":"#in"), 0.035, 32, 62, kBlack);
  }
}
