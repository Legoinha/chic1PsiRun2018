#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TSystem.h>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

void drawpthatweight(std::string inputname, std::string tag)
{
  std::string outputname = xjjc::str_replaceall(inputname, "rootfiles/", "");
  outputname = xjjc::str_replaceall(outputname, ".root", "");
  outputname = "plots/c"+outputname+".pdf";
  gSystem->Exec("mkdir -p plots");

  std::cout<<"==> "<<__FUNCTION__<<": "<<inputname<<std::endl;
  std::cout<<"<== "<<__FUNCTION__<<": "<<outputname<<std::endl;

  TFile* inf = TFile::Open(inputname.c_str());
  TH1F* hpthat = (TH1F*)inf->Get("hpthat");
  TH1F* hpthatweight = (TH1F*)inf->Get("hpthatweight");
  hpthat->Scale(1./hpthat->Integral());
  hpthatweight->Scale(1./hpthatweight->Integral());
  xjjroot::sethempty(hpthat);
  xjjroot::sethempty(hpthatweight);
  xjjroot::setthgrstyle(hpthat, xjjroot::colorlist_middle[3], 21, 0.5, xjjroot::colorlist_middle[3], 1, 2);
  xjjroot::setthgrstyle(hpthatweight, xjjroot::colorlist_middle[2], 21, 0.5, xjjroot::colorlist_middle[2], 1, 2);

  TLegend* leg = new TLegend(0.30, 0.25, 0.45, 0.25+2*0.045);
  xjjroot::setleg(leg, 0.040);
  leg->AddEntry(hpthat, "MC", "l");
  leg->AddEntry(hpthatweight, "MC weighted", "l");

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  gPad->SetLogy();
  hpthatweight->Draw("histe");
  hpthat->Draw("histe same");
  xjjroot::drawtex(0.88, 0.85, "PYTHIA (+ HYDJET)", 0.04, 33, 62);
  xjjroot::drawtex(0.88, 0.80, tag.c_str(), 0.04, 33);
  xjjroot::drawCMS("Simulation");
  leg->Draw();
  c->SaveAs(outputname.c_str());
}

int main(int argc, char* argv[])
{
  if(argc==3) { drawpthatweight(argv[1], argv[2]); return 0; }
  return 1;
}
