#include <iostream>
#include <iomanip>
#include <string>

#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>

namespace fitX
{
  class results
  {
  public:
    results(TFile* inf, std::string name = "");
    results() { ; }
    ~results() { ; }
    void print();
    float val(std::string opt) { return fval[opt]; }
    float err(std::string opt) { return ferr[opt]; }
    std::vector<std::string> vars() { return fvars; }
    std::string leg(std::string opt) { return flegs[opt]; }
    
  private:
    TFile* finf;
    std::string fname;
    TH1F *fhratio, *fhyield, *fhBenryield, *fhyieldpromptCorr;
    TEfficiency *fgrfprompt_a, *fgrfprompt_b, *fgreff_incl_a, *fgreff_incl_b;
    TH1D *fhscale_htnp_total_nominal_a, *fhscale_htnp_total_nominal_b;

    std::map<std::string, float> fval, ferr;
    void print_item(std::string var_, uint32_t precision, int stdwid, float scale = 1.);
    std::vector<std::string> fvars = {"yield", "Benryield", "eff", "tnp", "fprompt", "yieldpromptCorr"};
    std::map<std::string, std::string> flegs = 
      {
        std::pair<std::string, std::string>("yield", "Raw yield"),
        std::pair<std::string, std::string>("Benryield", "Raw yield (B-enr)"),
        std::pair<std::string, std::string>("eff", "Eff x Acc"),
        std::pair<std::string, std::string>("tnp", "TnP"),
        std::pair<std::string, std::string>("fprompt", "f_{prompt}"),
        std::pair<std::string, std::string>("yieldpromptCorr", "Corrected N_{prompt}"),
        std::pair<std::string, std::string>("ratio", "X(3872)/#psi(2S) ratio"),
      };
  };
}

fitX::results::results(TFile* inf, std::string name) : finf(inf), fname(name)
{
  fhratio = (TH1F*)inf->Get("hratio");
  fhyield = (TH1F*)inf->Get("hyield");
  fhBenryield = (TH1F*)inf->Get("hBenryield");
  fhyieldpromptCorr = (TH1F*)inf->Get("hyieldpromptCorr");
  fhscale_htnp_total_nominal_a = (TH1D*)inf->Get("hscale_htnp_total_nominal_a");
  fhscale_htnp_total_nominal_b = (TH1D*)inf->Get("hscale_htnp_total_nominal_b");
  fgrfprompt_a = (TEfficiency*)inf->Get("grfprompt_a");
  fgrfprompt_b = (TEfficiency*)inf->Get("grfprompt_b");
  fgreff_incl_a = (TEfficiency*)inf->Get("greff_incl_a");
  fgreff_incl_b = (TEfficiency*)inf->Get("greff_incl_b");

  fhratio->SetName(Form("%s%s", fhratio->GetName(), name.c_str()));
  fhyield->SetName(Form("%s%s", fhyield->GetName(), name.c_str()));
  fhBenryield->SetName(Form("%s%s", fhBenryield->GetName(), name.c_str()));
  fhyieldpromptCorr->SetName(Form("%s%s", fhyieldpromptCorr->GetName(), name.c_str()));
  fgrfprompt_a->SetName(Form("%s%s", fgrfprompt_a->GetName(), name.c_str()));
  fgrfprompt_b->SetName(Form("%s%s", fgrfprompt_b->GetName(), name.c_str()));
  fgreff_incl_a->SetName(Form("%s%s", fgreff_incl_a->GetName(), name.c_str()));
  fgreff_incl_b->SetName(Form("%s%s", fgreff_incl_b->GetName(), name.c_str()));

  fval["yield_a"] = fhyield->GetBinContent(fitX::ibin_a);
  ferr["yield_a"] = fhyield->GetBinError(fitX::ibin_a);
  fval["yield_b"] = fhyield->GetBinContent(fitX::ibin_b);
  ferr["yield_b"] = fhyield->GetBinError(fitX::ibin_b);
  
  fval["Benryield_a"] = fhBenryield->GetBinContent(fitX::ibin_a);
  ferr["Benryield_a"] = fhBenryield->GetBinError(fitX::ibin_a);
  fval["Benryield_b"] = fhBenryield->GetBinContent(fitX::ibin_b);
  ferr["Benryield_b"] = fhBenryield->GetBinError(fitX::ibin_b);

  fval["eff_a"] = fgreff_incl_a->GetEfficiency(fitX::ibin_a);
  ferr["eff_a u"] = fgreff_incl_a->GetEfficiencyErrorUp(fitX::ibin_a);
  ferr["eff_a d"] = fgreff_incl_a->GetEfficiencyErrorLow(fitX::ibin_a);
  ferr["eff_a"] = fgreff_incl_a->GetEfficiencyErrorLow(fitX::ibin_a);
  fval["eff_b"] = fgreff_incl_b->GetEfficiency(fitX::ibin_b);
  ferr["eff_b u"] = fgreff_incl_b->GetEfficiencyErrorUp(fitX::ibin_b);
  ferr["eff_b d"] = fgreff_incl_b->GetEfficiencyErrorLow(fitX::ibin_b);
  ferr["eff_b"] = fgreff_incl_b->GetEfficiencyErrorLow(fitX::ibin_b);

  fval["tnp_a"] = fhscale_htnp_total_nominal_a->GetBinContent(1);
  ferr["tnp_a"] = fhscale_htnp_total_nominal_a->GetBinError(1);
  fval["tnp_b"] = fhscale_htnp_total_nominal_b->GetBinContent(1);
  ferr["tnp_b"] = fhscale_htnp_total_nominal_b->GetBinError(1);
  
  fval["fprompt_a"] = fgrfprompt_a->GetEfficiency(fitX::ibin_a);
  ferr["fprompt_a u"] = fgrfprompt_a->GetEfficiencyErrorUp(fitX::ibin_a);
  ferr["fprompt_a d"] = fgrfprompt_a->GetEfficiencyErrorLow(fitX::ibin_a);
  ferr["fprompt_a"] = fgrfprompt_a->GetEfficiencyErrorLow(fitX::ibin_a);
  fval["fprompt_b"] = fgrfprompt_b->GetEfficiency(fitX::ibin_b);
  ferr["fprompt_b u"] = fgrfprompt_b->GetEfficiencyErrorUp(fitX::ibin_b);
  ferr["fprompt_b d"] = fgrfprompt_b->GetEfficiencyErrorLow(fitX::ibin_b);
  ferr["fprompt_b"] = fgrfprompt_b->GetEfficiencyErrorLow(fitX::ibin_b);

  fval["yieldpromptCorr_a"] = fhyieldpromptCorr->GetBinContent(fitX::ibin_a);
  ferr["yieldpromptCorr_a"] = fhyieldpromptCorr->GetBinError(fitX::ibin_a);
  fval["yieldpromptCorr_b"] = fhyieldpromptCorr->GetBinContent(fitX::ibin_b);
  ferr["yieldpromptCorr_b"] = fhyieldpromptCorr->GetBinError(fitX::ibin_b);

  fval["ratio"] = fhratio->GetBinContent(1);
  ferr["ratio"] = fhratio->GetBinError(1);
}

void fitX::results::print()
{
  int stdwid = 22;
  std::cout << std::left << "\e[33;1m --> [ \e[33;7m" << fname << "\e[0m\e[33;1m ]" << std::endl << std::string(stdwid*3+1, '-') << std::endl
            << std::setw(stdwid) << "| Var" << std::setw(stdwid) << "| psi(2S)" << std::setw(stdwid) << "| X(3872)" << "|" << std::endl 
            << std::string(stdwid*3+1, '-') << std::endl;
  print_item("yield", 0, stdwid);
  print_item("Benryield", 0, stdwid);
  print_item("eff", 1, stdwid, 1.e+3);
  print_item("tnp", 3, stdwid);
  print_item("fprompt", 2, stdwid);
  print_item("yieldpromptCorr", 0, stdwid);
  std::cout << std::setw(stdwid) << "| hratio" << std::setw(stdwid*2) << Form("| %.4f +- %.4f", fval["ratio"], ferr["ratio"]) << "|" << std::endl 
            << std::string(stdwid*3+1, '-') << std::endl
            << "\e[0m" << std::endl;
}

void fitX::results::print_item(std::string var_, uint32_t precision, int stdwid, float scale)
{
  std::string val_a_, val_b_;
  if(ferr.find(var_+" u") == ferr.end())
    {
      val_a_ = Form("| %.*f +- %.*f", precision, fval[var_+"_a"]*scale, precision, ferr[var_+"_a"]*scale);
      val_b_ = Form("| %.*f +- %.*f", precision, fval[var_+"_b"]*scale, precision, ferr[var_+"_b"]*scale);
    }
  else
    {
      val_a_ = Form("| %.*f + %.*f - %.*f", precision, fval[var_+"_a"]*scale, precision, ferr[var_+"_a u"]*scale, precision, ferr[var_+"_a d"]*scale);
      val_b_ = Form("| %.*f + %.*f - %.*f", precision, fval[var_+"_b"]*scale, precision, ferr[var_+"_b u"]*scale, precision, ferr[var_+"_b d"]*scale);
    }
  std::cout << std::setw(stdwid) << "| "+var_ << std::setw(stdwid) << val_a_ << std::setw(stdwid) << val_b_ << "|" << std::endl
            << std::string(stdwid*3+1, '-') << std::endl;
}
