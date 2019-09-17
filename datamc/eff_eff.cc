#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>

#include <string>

#include "fitX.h"
#include "project.h"
#include "MCefficiency.h"
#include "xjjcuti.h"

int n = 6;
void eff_eff(std::string input, std::string inputmcp_a, std::string inputmcp_b, 
             std::string cut, std::string cutgen, std::string output)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  std::cout<<cut<<std::endl;
  std::cout<<cutgen<<std::endl;

  std::string mcweight = "(pthatweight*Ncoll)";
  TFile* inf = TFile::Open(input.c_str());
  fitX::init(inf);
  TTree* ntf = (TTree*)inf->Get("ntf");
  std::vector<TString*> par_a(n), par_b(n), fun_a(n), fun_b(n);
  for(int i=0; i<n; i++)
    {
      ntf->SetBranchAddress(Form("par_a-%d", i+1), &(par_a[i]));
      ntf->SetBranchAddress(Form("par_b-%d", i+1), &(par_b[i]));
      ntf->SetBranchAddress(Form("fun_a-%d", i+1), &(fun_a[i]));
      ntf->SetBranchAddress(Form("fun_b-%d", i+1), &(fun_b[i]));
    }
  ntf->GetEntry(0);
  gROOT->cd();

  // TH1 must be defined after TTree declaration (some tricky issue) if no `gDirectory->cd("root:/");`
  //
  TTree* ntmixmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntmix"); if(!ntmixmcp_a) { return; }
  TTree* ntmixmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntmix"); if(!ntmixmcp_b) { return; }
  TTree* ntGenmcp_a = fitX::getnt(inputmcp_a, "Bfinder/ntGen"); if(!ntGenmcp_a) { return; }
  TTree* ntGenmcp_b = fitX::getnt(inputmcp_b, "Bfinder/ntGen"); if(!ntGenmcp_b) { return; }
  
  gDirectory->cd("root:/");

  std::string cutreco = Form("(%s) && Bpt>%f && Bpt<%f && TMath::Abs(By)>=%f && TMath::Abs(By)<%f && hiBin>=%f && hiBin<=%f", cut.c_str(), 
                             fitX::ptmincut, fitX::ptmaxcut, 
                             fitX::ymincut, fitX::ymaxcut, 
                             fitX::centmincut*2, fitX::centmaxcut*2);
  std::string cutmcreco = Form("%s && Bgen>=23333 && BgencollisionId==0", cutreco.c_str());
  std::string cutmcgen = Form("(%s) && Gpt>%f && Gpt<%f && TMath::Abs(Gy)>=%f && TMath::Abs(Gy)<%f && hiBin>=%f && hiBin<=%f && GisSignal==7 && GcollisionId==0", cutgen.c_str(), 
                              fitX::ptmincut, fitX::ptmaxcut, 
                              fitX::ymincut, fitX::ymaxcut, 
                              fitX::centmincut*2, fitX::centmaxcut*2);

  //
  std::vector<MCeff::MCefficiency*> mceff_a(n+1), mceff_b(n+1);
  for(int i=0; i<n+1; i++)
    {
      mceff_a[i] = new MCeff::MCefficiency(Form("_a-%d", i));
      mceff_b[i] = new MCeff::MCefficiency(Form("_b-%d", i));
      std::string mcweightBgenpt_a(mcweight), mcweightBgenpt_b(mcweight), mcweightGpt_a(mcweight), mcweightGpt_b(mcweight);

      if(i)
        {
          mceff_a[i]->heffmc()->SetTitle(fun_a[i-1]->Data());
          mceff_b[i]->heffmc()->SetTitle(fun_b[i-1]->Data());
          std::string ppar_a(*(par_a[i-1])), ppar_b(*(par_b[i-1]));
          if(xjjc::str_contains(ppar_a, "TMath::Sqrt")) { ppar_a = xjjc::str_replaceall(ppar_a, "TMath::Sqrt", "sqrt"); ppar_b = xjjc::str_replaceall(ppar_b, "TMath::Sqrt", "sqrt"); }
          if(xjjc::str_contains(ppar_a, "TMath::Exp")) { ppar_a = xjjc::str_replaceall(ppar_a, "TMath::Exp", "exp"); ppar_b = xjjc::str_replaceall(ppar_b, "TMath::Exp", "exp"); }
          mcweightBgenpt_a += ("*("+xjjc::str_replaceall(ppar_a, "x[0]", "Bgenpt")+")");
          mcweightBgenpt_b += ("*("+xjjc::str_replaceall(ppar_b, "x[0]", "Bgenpt")+")");
          mcweightGpt_a += ("*("+xjjc::str_replaceall(ppar_a, "x[0]", "Gpt")+")");
          mcweightGpt_b += ("*("+xjjc::str_replaceall(ppar_b, "x[0]", "Gpt")+")");
        }
      else
        {
          mceff_a[i]->heffmc()->SetTitle("PYTHIA+HYDJET");
          mceff_b[i]->heffmc()->SetTitle("PYTHIA+HYDJET");
        }

      std::cout << std::endl
                << " --> " << mcweightBgenpt_a << std::endl
                << " --> " << mcweightBgenpt_b << std::endl
                << " <-- " << mcweightGpt_a << std::endl
                << " <-- " << mcweightGpt_b << std::endl
                << std::endl;

      //
      std::cout<<" == mcp_a ==>"<<std::endl;
      ntmixmcp_a->Project(mceff_a[i]->heffmc()->GetName(), "Bpt", TCut(mcweightBgenpt_a.c_str())*TCut(cutmcreco.c_str()));
      fitX::printhist(mceff_a[i]->heffmc());
      ntmixmcp_a->Project(mceff_a[i]->heffmc_incl()->GetName(), Form("%d", fitX::ibin_a-1), TCut(mcweightBgenpt_a.c_str())*TCut(cutmcreco.c_str()));
      fitX::printhist(mceff_a[i]->heffmc_incl());
      std::cout<<" == mcp_b ==>"<<std::endl;
      ntmixmcp_b->Project(mceff_b[i]->heffmc()->GetName(), "Bpt", TCut(mcweightBgenpt_b.c_str())*TCut(cutmcreco.c_str()));
      fitX::printhist(mceff_b[i]->heffmc());
      ntmixmcp_b->Project(mceff_b[i]->heffmc_incl()->GetName(), Form("%d", fitX::ibin_b-1), TCut(mcweightBgenpt_b.c_str())*TCut(cutmcreco.c_str()));
      fitX::printhist(mceff_b[i]->heffmc_incl());

      //
      std::cout<<" == mcgenp_a ==>"<<std::endl;
      ntGenmcp_a->Project(mceff_a[i]->heffgen()->GetName(), "Gpt", TCut(mcweightGpt_a.c_str())*TCut(cutmcgen.c_str()));
      fitX::printhist(mceff_a[i]->heffgen());
      ntGenmcp_a->Project(mceff_a[i]->heffgen_incl()->GetName(), Form("%d", fitX::ibin_a-1), TCut(mcweightGpt_a.c_str())*TCut(cutmcgen.c_str()));
      fitX::printhist(mceff_a[i]->heffgen_incl());

      std::cout<<" == mcgenp_b ==>"<<std::endl;
      ntGenmcp_b->Project(mceff_b[i]->heffgen()->GetName(), "Gpt", TCut(mcweightGpt_b.c_str())*TCut(cutmcgen.c_str()));
      fitX::printhist(mceff_b[i]->heffgen());
      ntGenmcp_b->Project(mceff_b[i]->heffgen_incl()->GetName(), Form("%d", fitX::ibin_b-1), TCut(mcweightGpt_b.c_str())*TCut(cutmcgen.c_str()));
      fitX::printhist(mceff_b[i]->heffgen_incl());
    }

  std::string outputname = "rootfiles/"+output+"/eff_eff.root";
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  for(auto& mceff : mceff_a)
    {
      mceff->heffmc()->Write();
      mceff->heffmc_incl()->Write();
      mceff->heffgen()->Write();
      mceff->heffgen_incl()->Write();
    }
  for(auto& mceff : mceff_b)
    {
      mceff->heffmc()->Write();
      mceff->heffmc_incl()->Write();
      mceff->heffgen()->Write();
      mceff->heffgen_incl()->Write();
    }
  outf->cd();
  TTree* info = new TTree("info", "cut info");
  info->Branch("input", &input);
  info->Branch("inputmcp_a", &inputmcp_a);
  info->Branch("inputmcp_b", &inputmcp_b);
  info->Branch("cutmcreco", &cutmcreco);
  info->Branch("cutmcgen", &cutmcgen);
  info->Fill();
  info->Write();
  fitX::write();
  outf->Close();
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==7) { 
    eff_eff(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); return 0; }
  return 1;
}

