#include <vector>

// std::vector<float> bdtg = {-1, -0.9, -0.8, 0, 0.4, 0.6, 0.70, 0.76, 0.80}; int nbdtg = bdtg.size();
std::vector<float> bdtg = {-1   , -0.9 , -0.8, -0.5 , -0.4 , -0.3 , -0.2 , -0.1 , 0    , 0.1  , 0.2  , 0.3  , 0.4 , 
                           0.5  , 0.55 , 0.6 , 0.65 , 0.70 , 0.75 , 0.80 , 0.85 , 0.9  , 1};
std::vector<bool> pbdtg = {true , true , true, false, false, false, false, false, true , false, false, false, true, 
                           false, false, true, false, false, true , true , false, false, false};
 int nbdtg = bdtg.size();
std::vector<float> dls = {0, 0.8}; int ndls = dls.size();

void drawalltext()
{
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.042, "X(3872)", 0.038, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.90, 0.84, "p_{T} < 15 GeV/c", 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.042, "|y| < 1.5", 0.038, 32, 62);
}
