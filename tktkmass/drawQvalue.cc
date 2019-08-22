#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include <iostream>

#include "Qvalue.h"
#include "xjjrootuti.h"

void drawQvalue(std::string name)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::vector<TH1F*> hq(Qvalue::types.size(), 0);
  std::vector<TH1F*> hptimb(Qvalue::types.size(), 0);
  TLegend* leg = new TLegend(0.60, 0.85-0.04*4, 0.85, 0.85-0.04);
  xjjroot::setleg(leg, 0.038);
  for(int i=0; i<Qvalue::types.size(); i++)
    {
      auto ff = TFile::Open(Form("rootfiles/%s/Qvalue_%s.root", name.c_str(), Qvalue::types[i].c_str()));
      hq[i] = (TH1F*)ff->Get("hq");
      hq[i]->SetName(Form("hq_%s", Qvalue::types[i].c_str()));
      hq[i]->Scale(1./hq[i]->Integral());
      xjjroot::setthgrstyle(hq[i], Qvalue::colors[i], 21, 0.7, Qvalue::colors[i], Qvalue::line[i], 5);
      xjjroot::sethempty(hq[i], 0, 0);
      hq[i]->GetXaxis()->SetNdivisions(-505);

      hptimb[i] = (TH1F*)ff->Get("hptimb");
      hptimb[i]->SetName(Form("hptimb_%s", Qvalue::types[i].c_str()));
      hptimb[i]->Scale(1./hptimb[i]->Integral());
      xjjroot::setthgrstyle(hptimb[i], Qvalue::colors[i], 21, 0.7, Qvalue::colors[i], Qvalue::line[i], 5);
      xjjroot::sethempty(hptimb[i], 0, 0);
      hptimb[i]->GetXaxis()->SetNdivisions(-505);

      leg->AddEntry(hq[i], Qvalue::legtrs[i].c_str(), "l");
    }

  xjjroot::setgstyle(1);
  TCanvas* chq = new TCanvas("chq", "", 600, 600);
  hq[1]->SetMaximum(hq[1]->GetMaximum()*1.3);
  hq[1]->Draw("histe");
  hq[0]->Draw("histe same");
  hq[2]->Draw("histe same");
  leg->Draw();
  xjjroot::drawtex(0.90, 0.85, "PYTHIA + HYDJET", 0.038, 33, 62);
  xjjroot::drawtex(0.24, 0.85, "|y| < 1.6", 0.038, 13, 42);
  xjjroot::drawtex(0.24, 0.85-0.04, "Cent. 0-90%", 0.038, 13, 42);
  xjjroot::drawCMS("Simulation");
  std::string outputnamehq = Form("plots/%s/mcQvalue.pdf", name.c_str());
  xjjroot::mkdir(outputnamehq.c_str());
  chq->SaveAs(outputnamehq.c_str());

  TCanvas* chptimb = new TCanvas("chptimb", "", 600, 600);
  hptimb[1]->SetMaximum(hptimb[1]->GetMaximum()*1.3);
  hptimb[1]->Draw("histe");
  hptimb[0]->Draw("histe same");
  hptimb[2]->Draw("histe same");
  leg->Draw();
  xjjroot::drawtex(0.90, 0.85, "PYTHIA + HYDJET", 0.038, 33, 62);
  xjjroot::drawtex(0.24, 0.85, "|y| < 1.6", 0.038, 13, 42);
  xjjroot::drawtex(0.24, 0.85-0.04, "Cent. 0-90%", 0.038, 13, 42);
  xjjroot::drawCMS("Simulation");
  std::string outputnamehptimb = Form("plots/%s/mcPtimb.pdf", name.c_str());
  xjjroot::mkdir(outputnamehptimb.c_str());
  chptimb->SaveAs(outputnamehptimb.c_str());
}

int main(int argc, char* argv[])
{
  if(argc==2) { drawQvalue(argv[1]); return 0; }
  return 1;
}

