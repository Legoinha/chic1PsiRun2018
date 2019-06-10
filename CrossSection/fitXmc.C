#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

#include "xjjcuti.h"

float weight_ss = 1./0.681;

void fitXmc(std::string input, std::string input_ss, std::string inputmc_a, std::string inputmc_b, 
            std::string cut, std::string output)
{

  std::cout<<cut<<std::endl;

  bool is_ss = false;

  std::cout<<" -- data >>>>"<<std::endl;
  TFile* inf = TFile::Open(input.c_str());
  if(!inf->IsOpen()) return;
  TH1F* h = new TH1F("h", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h->Sumw2();
  TTree* ntmix = (TTree*)inf->Get("Bfinder/ntmix");
  ntmix->AddFriend("hiEvtAnalyzer/HiTree");
  ntmix->AddFriend("hltanalysis/HltTree");
  ntmix->AddFriend("skimanalysis/HltTree");
  ntmix->AddFriend("dataset/mva");
  ntmix->Project("h", "Bmass", TCut(Form("(%s)", cut.c_str())));
  std::cout<<h->GetEntries()<<std::endl;
  std::cout<<" -- same-sign >>>>"<<std::endl;
  TFile* inf_ss = TFile::Open(input_ss.c_str());
  TH1F* h_ss = new TH1F("h_ss", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3), fitX::NBIN, fitX::BIN_MIN, fitX::BIN_MAX); h_ss->Sumw2();
  if(inf_ss->IsOpen()) 
    {  
      is_ss = true;
      TTree* ntmix_ss = (TTree*)inf_ss->Get("Bfinder/ntmix");
      ntmix_ss->AddFriend("hiEvtAnalyzer/HiTree");
      ntmix_ss->AddFriend("hltanalysis/HltTree");
      ntmix_ss->AddFriend("skimanalysis/HltTree");
      ntmix_ss->AddFriend("dataset/mva");
      ntmix_ss->Project("h_ss", "Bmass", TCut(Form("(%s)", cut.c_str())));
      std::cout<<h_ss->GetEntries()<<std::endl;
      h_ss->Scale(weight_ss);
    }

  std::cout<<" -- mc_a >>>>"<<std::endl;
  TFile* infmc_a = TFile::Open(inputmc_a.c_str());
  if(!infmc_a->IsOpen()) return;
  TH1F* hmc_a = new TH1F("hmc_a", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3), fitX::NBIN_L, fitX::BIN_MIN_L, fitX::BIN_MAX_L); hmc_a->Sumw2();
  TTree* ntmixmc_a = (TTree*)infmc_a->Get("Bfinder/ntmix");
  ntmixmc_a->AddFriend("hiEvtAnalyzer/HiTree");
  ntmixmc_a->AddFriend("hltanalysis/HltTree");
  ntmixmc_a->AddFriend("skimanalysis/HltTree");
  ntmixmc_a->AddFriend("dataset/mva");
  ntmixmc_a->Project("hmc_a", "Bmass", TCut("1")*TCut(Form("(%s) && Bgen==23333 && BgencollisionId==0", cut.c_str())));
  std::cout<<hmc_a->GetEntries()<<std::endl;
  hmc_a->Scale(hmc_a->GetEntries()/hmc_a->Integral());

  std::cout<<" -- mc_b >>>>"<<std::endl;
  TFile* infmc_b = TFile::Open(inputmc_b.c_str());
  if(!infmc_b->IsOpen()) return;
  TH1F* hmc_b = new TH1F("hmc_b", Form(";m_{#mu#mu#pi#pi} (GeV/c^{2});Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3), fitX::NBIN_H, fitX::BIN_MIN_H, fitX::BIN_MAX_H); hmc_b->Sumw2();
  TTree* ntmixmc_b = (TTree*)infmc_b->Get("Bfinder/ntmix");
  ntmixmc_b->AddFriend("hiEvtAnalyzer/HiTree");
  ntmixmc_b->AddFriend("hltanalysis/HltTree");
  ntmixmc_b->AddFriend("skimanalysis/HltTree");
  ntmixmc_b->AddFriend("dataset/mva");
  ntmixmc_b->Project("hmc_b", "Bmass", TCut("1")*TCut(Form("(%s) && Bgen==23333 && BgencollisionId==0", (xjjc::str_contains(cut,"Blxy > 0.1")?xjjc::str_replaceall(cut, "&& Blxy > 0.1", ""):cut).c_str())));
  std::cout<<hmc_b->GetEntries()<<std::endl;
  hmc_b->Scale(hmc_b->GetEntries()/hmc_b->Integral());

  // fit + yield
  std::vector<TF1*> funs = fitX::fit(h, (is_ss?h_ss:0), hmc_a, hmc_b, output, false, true); // fix mean = false
  float ysig1 = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
  float ysig1err = funs[0]->GetParError(5)*ysig1/funs[0]->GetParameter(5);
  float ysig2 = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
  float ysig2err = funs[0]->GetParError(10)*ysig2/funs[0]->GetParameter(10);
  float ybkg1 = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH;
  float ybkg2 = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH;

  // yield
  TH1F* hyield_a = new TH1F("hyield_a", "", 5, 0, 5); hyield_a->Sumw2();
  xjjroot::setthgrstyle(hyield_a, fitX::color_a, 20, 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
  hyield_a->SetBinContent(2, ysig1);
  hyield_a->SetBinError(2, ysig1err);
  TH1F* hyield_b = new TH1F("hyield_b", "", 5, 0, 5); hyield_b->Sumw2();
  xjjroot::setthgrstyle(hyield_b, fitX::color_b, 20, 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
  hyield_b->SetBinContent(4, ysig2);
  hyield_b->SetBinError(4, ysig2err);
  float ymax = 300.;
  TH2F* hempty = new TH2F("hempty", ";;Raw Yield", 5, 0, 5, 10, 0, ymax);
  xjjroot::sethempty(hempty, 0, 0.3);
  hempty->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hempty->GetXaxis()->SetBinLabel(4, "X(3872)");
  hempty->GetXaxis()->SetLabelSize(hempty->GetXaxis()->GetLabelSize()*1.5);
  xjjroot::setgstyle();
  TCanvas* cy = new TCanvas("cy", "", 600, 600);
  hempty->Draw();
  hyield_a->Draw("ple same");
  hyield_b->Draw("ple same");
  xjjroot::drawtex(0.42, ysig1/ymax + 0.2, Form("%.0f #pm %.0f", ysig1, ysig1err), 0.042, 22, 62, fitX::color_a);
  xjjroot::drawtex(0.72, ysig2/ymax + 0.2, Form("%.0f #pm %.0f", ysig2, ysig2err), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawtex(0.92, 0.84, "|y| < 1.5", 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.77, "p_{T} > 15 GeV/c", 0.042, 32, 62);
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  cy->SaveAs(Form("plots/cyield_%s.pdf", output.c_str()));

  // background
  TH1F* hbkg_a = new TH1F("hbkg_a", "", 5, 0, 5);
  xjjroot::setthgrstyle(hbkg_a, fitX::color_a, 20, 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
  hbkg_a->SetBinContent(2, ybkg1);
  hbkg_a->SetBinError(2, 0);
  TH1F* hbkg_b = new TH1F("hbkg_b", "", 5, 0, 5);
  xjjroot::setthgrstyle(hbkg_b, fitX::color_b, 20, 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
  hbkg_b->SetBinContent(4, ybkg2);
  hbkg_b->SetBinError(4, 0);
  float ymaxbkg = ybkg1 * 2.2;
  TH2F* hemptybkg = new TH2F("hemptybkg", ";;Bkg in Signal Region", 5, 0, 5, 10, 0, ymaxbkg);
  xjjroot::sethempty(hemptybkg, 0, 0.3);
  hemptybkg->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hemptybkg->GetXaxis()->SetBinLabel(4, "X(3872)");
  hemptybkg->GetXaxis()->SetLabelSize(hemptybkg->GetXaxis()->GetLabelSize()*1.5);
  xjjroot::setgstyle();
  TCanvas* cb = new TCanvas("cb", "", 600, 600);
  hemptybkg->Draw();
  hbkg_a->Draw("hist same");
  hbkg_b->Draw("hist same");
  xjjroot::drawtex(0.42, ybkg1/ymaxbkg + 0.2, Form("%.0f", ybkg1), 0.042, 22, 62, fitX::color_a);
  xjjroot::drawtex(0.72, ybkg2/ymaxbkg + 0.2, Form("%.0f", ybkg2), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawtex(0.92, 0.84, "|y| < 1.5", 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.77, "p_{T} > 15 GeV/c", 0.042, 32, 62);
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  cb->SaveAs(Form("plots/cbkg_%s.pdf", output.c_str()));

  TH1F* hyield = (TH1F*)hyield_a->Clone("hyield");
  hyield->Add(hyield_b);
  TH1F* hbkg = (TH1F*)hbkg_a->Clone("hbkg");
  hbkg->Add(hbkg_b);

  TFile* outf = new TFile(Form("rootfiles/root_fitXmc_%s.root", output.c_str()), "recreate");
  outf->cd();
  h->Write();
  h_ss->Write();
  hmc_a->Write();
  hmc_b->Write();
  hyield_a->Write();
  hyield_b->Write();
  hyield->Write();
  hbkg_a->Write();
  hbkg_b->Write();
  hbkg->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==7) { fitXmc(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); return 0; }
  return 1;
}
