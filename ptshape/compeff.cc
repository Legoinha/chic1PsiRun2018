#include <TTree.h>
#include <TFile.h>
#include <TCut.h>
#include <TSystem.h>
#include "fitX.h"
#include "xjjcuti.h"
#include "project.h"

#include "MCefficiency.h"
#include "HEPData-ins1495026-v1-csv/ppATLAS.h"

void compeff(std::string inputmcp, std::string type, std::string cut, std::string cutgen, std::string outputname, std::string inputdir)
{
  int ibin = 0;
  if(xjjc::str_contains(type, "romptPsi")) ibin = fitX::ibin_a;
  if(xjjc::str_contains(type, "romptX")) ibin = fitX::ibin_b;
  std::cout<<cut<<std::endl;
  ppRef::ppATLAS getpp(inputdir);
  std::string mcweight = "(pthatweight*Ncoll)";
  std::string ptweight = getpp.formula[type];
  std::string ptweightreco = xjjc::str_replaceall(ptweight, "x", "Bgenpt"); ptweightreco = xjjc::str_replaceall(ptweightreco, "EBgenptp", "Exp");
  std::string ptweightgen = xjjc::str_replaceall(ptweight, "x", "Gpt"); ptweightgen = xjjc::str_replaceall(ptweightgen, "EGptp", "Exp");

  TTree* ntmixmcp = fitX::getnt(inputmcp, "ntmix"); if(!ntmixmcp) { return; }
  TTree* ntGenmcp = fitX::getnt(inputmcp, "ntGen"); if(!ntGenmcp) { return; }
  MCeff::MCefficiency mceff("", 1);
  MCeff::MCefficiency mceffweight("_weight", 1);
  std::string cutreco = Form("(%s) && Bpt>%f && Bpt<%f && TMath::Abs(By)<%f", cut.c_str(), fitX::ptmincut, fitX::ptmaxcut, fitX::ycut);
  std::string cutmcreco = Form("%s && Bgen>=23333 && BgencollisionId==0", cutreco.c_str());
  std::string cutmcgen = Form("(%s) && Gpt>%f && Gpt<%f && TMath::Abs(Gy)<%f && GisSignal==7 && GcollisionId==0", cutgen.c_str(), fitX::ptmincut, fitX::ptmaxcut, fitX::ycut);
  std::cout<<" == mcp ==>"<<std::endl;
  ntmixmcp->Project(mceff.heffmc->GetName(), "Bpt", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp->Project(mceff.heffmc_incl->GetName(), "Bpt", TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  // ntmixmcp->Project(mceff.heffmc_incl->GetName(), Form("%d", ibin-1), TCut(mcweight.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp->Project(mceffweight.heffmc->GetName(), "Bpt", TCut(mcweight.c_str())*TCut(ptweightreco.c_str())*TCut(cutmcreco.c_str()));
  ntmixmcp->Project(mceffweight.heffmc_incl->GetName(), "Bpt", TCut(mcweight.c_str())*TCut(ptweightreco.c_str())*TCut(cutmcreco.c_str()));
  // ntmixmcp->Project(mceffweight.heffmc_incl->GetName(), Form("%d", ibin-1), TCut(mcweight.c_str())*TCut(ptweightreco.c_str())*TCut(cutmcreco.c_str()));
  std::cout<<" == mcgenp ==>"<<std::endl;
  ntGenmcp->Project(mceff.heffgen->GetName(), "Gpt", TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));
  ntGenmcp->Project(mceff.heffgen_incl->GetName(), "Gpt", TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));
  // ntGenmcp->Project(mceff.heffgen_incl->GetName(), Form("%d", ibin-1), TCut(mcweight.c_str())*TCut(cutmcgen.c_str()));
  ntGenmcp->Project(mceffweight.heffgen->GetName(), "Gpt", TCut(mcweight.c_str())*TCut(ptweightgen.c_str())*TCut(cutmcgen.c_str()));
  ntGenmcp->Project(mceffweight.heffgen_incl->GetName(), "Gpt", TCut(mcweight.c_str())*TCut(ptweightgen.c_str())*TCut(cutmcgen.c_str()));
  // ntGenmcp->Project(mceffweight.heffgen_incl->GetName(), Form("%d", ibin-1), TCut(mcweight.c_str())*TCut(ptweightgen.c_str())*TCut(cutmcgen.c_str()));
  std::cout<<mceffweight.heffgen_incl->GetEntries()<<" "<<mceffweight.heffgen_incl->Integral()<<std::endl;

  gSystem->Exec(Form("mkdir -p %s", xjjc::str_replaceall(outputname, xjjc::str_divide(outputname, "/").back(), "").c_str()));
  TFile* outf = new TFile(Form("%s", outputname.c_str()), "recreate");
  outf->cd();
  mceff.heffmc->Write();
  mceff.heffmc_incl->Write();
  mceff.heffgen->Write();
  mceff.heffgen_incl->Write();
  mceffweight.heffmc->Write();
  mceffweight.heffmc_incl->Write();
  mceffweight.heffgen->Write();
  mceffweight.heffgen_incl->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==7) { compeff(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); return 0; }
  return 1;
}
