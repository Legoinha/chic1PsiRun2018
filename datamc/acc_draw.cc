#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <string>
#include <iostream>

#include "xjjrootuti.h"
#include "fitX.h"

namespace acc
{
  class elements
  {
  public:
    elements(std::string input, std::string name_, std::string tex_,  Color_t color_);
    TH1F* hacc() { return fhacc; }
    float nominal() { return fnominal; }
    std::string name() { return fname; }
    std::string tex() { return ftex; }
    double mean() { return fmean; }
    double rms() { return frms; }
    int entries() { return fentries; }
    bool valid() { return fvalid; }
    void draw(float x, float y);
  private:
    bool fvalid;
    TH1F* fhacc;
    float fnominal;
    double fmean;
    double frms;
    int fentries;
    std::string fname;
    std::string ftex;
    TLegend* fleg;
    Color_t fcolor;
  };
}

void acc_draw(std::string inputpt, std::string inputabsy, std::string output)
{
  bool aorb = !xjjc::str_contains(inputpt, "_b.root");
  std::string ptex = aorb?"#psi(2S)":"X(3872)";
  std::string taorb = aorb?"a":"b";
  acc::elements* lpt = new acc::elements(inputpt, "pt", "p_{T}", xjjroot::mycolor_middle["orange"]);
  acc::elements* labsy = new acc::elements(inputabsy, "absy", "y", xjjroot::mycolor_middle["cyan"]);
  float ymax = std::max(lpt->hacc()->GetMaximum(), labsy->hacc()->GetMaximum());
  lpt->hacc()->SetMinimum(0);
  lpt->hacc()->SetMaximum(ymax*1.3);
  labsy->hacc()->SetMinimum(0);
  labsy->hacc()->SetMaximum(ymax*1.3);

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  labsy->hacc()->Draw("hist");
  lpt->draw(0.24, 0.60);
  labsy->draw(0.64, 0.60);
  xjjroot::drawtex(0.24, 0.85, ptex.c_str(), 0.038, 12, 62);
  xjjroot::drawtex(0.24, 0.84-0.04, "Toy MC study", 0.038, 12, 52);
  xjjroot::drawtex(0.24, 0.84-0.04*2, "Nominal", 0.036, 12, 62, kGray+1);
  xjjroot::drawtex(0.24, 0.84-0.04*3, Form("#alpha = %.5f", lpt->nominal()), 0.036, 12, 42, kGray+1);
  fitX::drawkinematics();
  xjjroot::drawCMS();
  std::string outputname = "plots/"+output+"/acceptance/caccsyst_"+taorb+".pdf";
  xjjroot::mkdir(outputname);
  c->RedrawAxis();
  c->SaveAs(outputname.c_str());

  float per = sqrt((lpt->mean()-lpt->nominal())*(lpt->mean()-lpt->nominal()) + (labsy->mean()-labsy->nominal())*(labsy->mean()-labsy->nominal()));
  per /= lpt->nominal();
  std::cout<<"Acceptance & "<<Form("%.1f", per*1.e+2)<<"\\% \\\\"<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==4) { acc_draw(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}

acc::elements::elements(std::string input, std::string name_, std::string tex_, Color_t color_) : fname(name_), ftex(tex_), fcolor(color_)
{
  TFile* inf = TFile::Open(input.c_str());
  fhacc = (TH1F*)inf->Get("hacc"); 
  fhacc->SetName(Form("hacc_%s", fname.c_str()));
  xjjroot::printhist(fhacc);
  fhacc->GetXaxis()->SetNdivisions(505);
  fmean = fhacc->GetMean();
  frms = fhacc->GetRMS();
  fentries = fhacc->GetEntries();
  xjjroot::sethempty(fhacc, 0, 0);
  xjjroot::setthgrstyle(fhacc, fcolor, 21, 0.8, fcolor, 1, 4);
  TH1F* hnominal = (TH1F*)inf->Get("hnominal");
  fnominal = hnominal->GetBinContent(1);
  delete hnominal;
  fvalid = true;
  
}

void acc::elements::draw(float x, float y)
{
  xjjroot::drawline(fnominal, fhacc->GetMinimum(), fnominal, fhacc->GetMaximum(), kGray+1, 7, 3);
  fhacc->Draw("hist same");
  xjjroot::drawtex(x, y+0.001, Form("%s variation", ftex.c_str()), 0.036, 13, 62, fcolor);
  xjjroot::drawtex(x, y-0.043, "Entries", 0.036, 13, 42, fcolor); xjjroot::drawtex(x+0.24, y-0.043, Form("%d", fentries), 0.036, 33, 42, fcolor);
  xjjroot::drawtex(x, y-0.043*2, "Mean", 0.036, 13, 42, fcolor);  xjjroot::drawtex(x+0.24, y-0.043*2, Form("%.5f", fmean), 0.036, 33, 42, fcolor);
  xjjroot::drawtex(x, y-0.043*3, "RMS", 0.036, 13, 42, fcolor);   xjjroot::drawtex(x+0.24, y-0.043*3, Form("%.5f", frms), 0.036, 33, 42, fcolor);
}

