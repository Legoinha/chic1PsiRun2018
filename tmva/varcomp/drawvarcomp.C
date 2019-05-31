#include "varcomp.h"
#include "xjjcuti.h"
#include "xjjrootuti.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>

int ncolumn = 3, nrow = 2, npad = ncolumn*nrow;
void drawvarcomp(std::string inputdir, std::vector<std::string> legtitle)
{
  xjjroot::setgstyle(1);
  std::string outputname = xjjc::str_replaceallspecial(inputdir);

  TFile* inf = TFile::Open(Form("%s.root", inputdir.c_str()));
  TTree* note = (TTree*)inf->Get("note");
  std::string* texpt = 0; note->SetBranchAddress("texpt", &texpt); note->GetEntry(0);
  std::vector<varcomp::variable*> var;
  int n = 0;
  while(true)
    {
      varcomp::variable* ivar = new varcomp::variable(inf, Form("nt%d", n));
      if(ivar->isvalid()) { var.push_back(ivar); }
      else break;
      var[n]->sethist(n);
      n++;
    }
  if(!var.size()) return;
  int nvar = var[0]->nvar;

  TLegend* leg = new TLegend(0.55, 0.85-n*0.06, 0.88, 0.85);
  xjjroot::setleg(leg, 0.048);
  for(int k=0; k<n; k++) { leg->AddEntry(var[k]->hist[0], legtitle[k].c_str(), "f"); }

  int vv = 0, cc = 0;
  TCanvas* c;
  while(vv < nvar)
    {
      if(vv%npad == 0) 
        {
          nrow = std::min((nvar-vv)/ncolumn+1, nrow);
          c = new TCanvas("c", "", 400*ncolumn, 400*nrow);
          c->Divide(ncolumn, nrow);
        }
      float ymax = 0;
      for(auto& v : var)
        {
          if(v->hist[vv]->GetMaximum() > ymax)
            { ymax = v->hist[vv]->GetMaximum(); }
        }
      var[0]->hist[vv]->SetMaximum(ymax*1.2);
      var[0]->hist[vv]->SetMinimum(0.);
      xjjroot::sethempty(var[0]->hist[vv]);
      c->cd(vv%6 + 1);
      for(int k=n-1; k>=0; k--)
        {
          var[k]->hist[vv]->Draw(Form("histe %s", k==n-1?"same":""));
        }
      xjjroot::drawCMS();
      leg->Draw();
      xjjroot::drawtex(0.22, 0.82, (*texpt + " GeV/c").c_str(), 0.045, 12, 62);
      if((vv+1)%npad == 0 || vv == nvar-1)
        {
          c->SaveAs(Form("plots/cvar_%s_%d.pdf", outputname.c_str(), cc));
          cc++;
          delete c; 
        }
      vv++;
    }  
}

int main(int argc, char* argv[])
{
  if(argc > 2)
    {
      std::vector<std::string> legtitle;
      for(int a=2; a<argc; a++) { legtitle.push_back(argv[a]); }
      drawvarcomp(argv[1], legtitle);
      return 0; 
    }
  std::cout<<__FUNCTION__<<": invalid argument number."<<std::endl;
  return 1;
}


