#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include <string>

#include "xjjrootuti.h"
#include "xjjcuti.h"
#include "fitX.h"

namespace acc
{
  float getacc(float p0, float p1, std::vector<float> &Gvar, std::vector<float> &weight, std::vector<int> &isAcc);
}
void acc_acc(std::string inputskim, std::string inputpars, std::string output, std::string type, std::string aorb)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if((aorb!="a" && aorb!="b") ||
     (type!="pt" && type!="absy"))
    { std::cout<<"error: invalid arguments."<<std::endl; return; }

  TFile* infskim = TFile::Open(inputskim.c_str());
  TTree* gnt = (TTree*)infskim->Get("gnt");
  int ng = gnt->GetEntries();
  TFile* infpars = TFile::Open(inputpars.c_str());
  fitX::init(infpars);
  TTree* par = (TTree*)infpars->Get("par");
  int np = par->GetEntries();

  gnt->SetBranchStatus("*", 0);
  float Gpt;    gnt->SetBranchStatus("Gpt", 1);    gnt->SetBranchAddress("Gpt", &Gpt);
  float Gy;     gnt->SetBranchStatus("Gy", 1);     gnt->SetBranchAddress("Gy", &Gy);
  float weight; gnt->SetBranchStatus("weight", 1); gnt->SetBranchAddress("weight", &weight);
  bool isAcc;   gnt->SetBranchStatus("isAcc", 1);  gnt->SetBranchAddress("isAcc", &isAcc);
  float p0;     par->SetBranchAddress(Form("%s0", aorb.c_str()), &p0);
  float p1;     par->SetBranchAddress(Form("%s1", aorb.c_str()), &p1);
  bool ispt = (type=="pt");

  // int nn = ng/10+1;
  std::vector<float> vGvar, vweight;
  std::vector<int> visAcc;
  for(int j=0; j<ng; j+=10) // +=
    {
      gnt->GetEntry(j);
      vGvar.push_back(ispt?Gpt:fabs(Gy));
      vweight.push_back(weight);
      visAcc.push_back((int)isAcc);
    }

  float nominal = acc::getacc(1, 0, vGvar, vweight, visAcc);
  float err = 0.07; //
  TH1F* hacc = new TH1F("hacc", ";#alpha;", 200, nominal-err, nominal+err);
  TH1F* hnominal = new TH1F("hnominal", ";;", 1, 0, 1);
  hnominal->Fill(0., nominal);
  for(int i=0; i<np; i++)
    {
      if(i%1000==0) xjjc::progressbar(i, np);
      par->GetEntry(i);
      hacc->Fill(acc::getacc(p0, p1, vGvar, vweight, visAcc));
    }
  xjjc::progressbar_summary(np);
  std::string outputname = ("rootfiles/"+output+"/acc_acc_"+aorb+".root");
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hnominal->Write();
  hacc->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==6) { acc_acc(argv[1], argv[2], argv[3], argv[4], argv[5]); return 0; }
  return 1;
}

float acc::getacc(float p0, float p1, std::vector<float> &Gvar, std::vector<float> &weight, std::vector<int> &isAcc)
{
  float num = 0, den = 0;
  int ng = Gvar.size();
  for(int j=0; j<ng; j++)
    {
      // float ratio = (p0+p1*(Gvar[j]));
      float ratio = exp(p0+p1*(Gvar[j]));
      float w = ratio * (weight[j]); // 
      den += w;
      num += (w*isAcc[j]);
    }
  return num/den;
}

