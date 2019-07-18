#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

#include "fitX.h"
#include "project.h"

namespace tktkmass
{
  const int NBIN = 10;
  const float BIN_MIN = 0.6, BIN_MAX = 0.8;
  float BIN_WIDTH = (BIN_MAX-BIN_MIN)*1./NBIN;
}

void tktkmass_main(std::string input, std::string inputmcp, std::string cut, std::string output, std::string name)
{

  std::cout<<cut<<std::endl;
  
  std::string mcweight = "(pthatweight*Ncoll)";

  //
  TTree* ntmix = (TTree*)fitX::getnt(input, "ntmix"); if(!ntmix) { return; }
  TTree* ntmixmcp = fitX::getnt(inputmcp, "ntmix"); if(!ntmixmcp) { return; }

  TH1F* htktkmass_sigreg = new TH1F("htktkmass_sigreg", Form(";m_{#pi#pi} (GeV/c^{2});Entries / %.0f MeV", tktkmass::BIN_WIDTH*1.e+3), tktkmass::NBIN, tktkmass::BIN_MIN, tktkmass::BIN_MAX); 
  htktkmass_sigreg->Sumw2();
  TH1F* htktkmass_sidbnd = new TH1F("htktkmass_sidbnd", Form(";m_{#pi#pi} (GeV/c^{2});Entries / %.0f MeV", tktkmass::BIN_WIDTH*1.e+3), tktkmass::NBIN, tktkmass::BIN_MIN, tktkmass::BIN_MAX); 
  htktkmass_sidbnd->Sumw2();
  TH1F* htktkmass_mcp = new TH1F("htktkmass_mcp", Form(";m_{#pi#pi} (GeV/c^{2});Entries / %.0f MeV", tktkmass::BIN_WIDTH*1.e+3), tktkmass::NBIN, tktkmass::BIN_MIN, tktkmass::BIN_MAX); 
  htktkmass_mcp->Sumw2();
  
  std::string cutreco = Form("(%s) && Bpt>%f && TMath::Abs(By)<%f", cut.c_str(), fitX::ptcut, fitX::ycut);
  std::string cutmcreco = Form("%s && Bgen>=23333 && BgencollisionId==0", cutreco.c_str(), fitX::ptcut, fitX::ycut);

  std::cout<<" == data ==>"<<std::endl;
  ntmix->Project("htktkmass_sigreg", "Btktkmass", TCut(Form("%s && TMath::Abs(Bmass-3.8719)<0.01", cutreco.c_str())));
  std::cout<<htktkmass_sigreg->GetEntries()<<std::endl;
  ntmix->Project("htktkmass_sidbnd", "Btktkmass", TCut(Form("%s && TMath::Abs(Bmass-3.8719)>0.02 && TMath::Abs(Bmass-3.8719)<0.12", cutreco.c_str())));
  std::cout<<htktkmass_sidbnd->GetEntries()<<std::endl;

  std::cout<<" == mcp ==>"<<std::endl;
  ntmixmcp->Project("htktkmass_mcp", "Btktkmass", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  std::cout<<htktkmass_mcp->GetEntries()<<std::endl;

  htktkmass_sidbnd->Scale(1./10);
  TH1F* htktkmass_sigsub = (TH1F*)htktkmass_sigreg->Clone("htktkmass_sigsub");
  htktkmass_sigsub->Add(htktkmass_sidbnd, -1.);

  htktkmass_mcp->Scale(1./htktkmass_mcp->Integral());
  htktkmass_sigsub->Scale(1./htktkmass_sigsub->Integral());
  htktkmass_sigsub->SetMinimum(0);
  xjjroot::setthgrstyle(htktkmass_mcp, kRed+1, 21, 0.5, kRed+1, 2, 2);
  xjjroot::setthgrstyle(htktkmass_sigsub, kBlack, 20, 1.2, kBlack, 1, 2);
  xjjroot::sethempty(htktkmass_mcp);
  xjjroot::sethempty(htktkmass_sigsub);

  xjjroot::setgstyle(3);
  TCanvas* ctktkmass = new TCanvas("ctktkmass", "", 600, 600);
  htktkmass_sigsub->Draw("ple");
  htktkmass_mcp->Draw("hist e same");
  xjjroot::drawCMS();
  ctktkmass->SaveAs(Form("plots/%s/ctktkmass.pdf", name.c_str()));

  TFile* outf = new TFile(Form("%s.root", output.c_str()), "recreate");
  outf->cd();
  htktkmass_sigreg->Write();
  htktkmass_sidbnd->Write();
  htktkmass_sigsub->Write();
  htktkmass_mcp->Write();
  outf->cd();
  TTree* info = new TTree("info", "cut info");
  info->Branch("input", &input);
  info->Branch("inputmcp", &inputmcp);
  info->Branch("cutreco", &cutreco);
  info->Branch("cutmcreco", &cutmcreco);
  info->Fill();
  info->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==6) { tktkmass_main(argv[1], argv[2], argv[3], argv[4], argv[5]); return 0; }
  return 1;
}
