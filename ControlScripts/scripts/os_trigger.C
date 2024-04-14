#include "RooGlobalFunc.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
//#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TH1F.h"
#include "TH2F.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>

using namespace RooFit;

void os_trigger () 
{
    TFile *TreeFile_h = new TFile("LOCAL_ANALYSIS_ztau3mutauh_skimmer_default.root","READ");
    TFile *TreeFile_mu = new TFile("LOCAL_ANALYSIS_ztau3mutaumu_skimmer_default.root","READ");
    TFile *TreeFile_e = new TFile("LOCAL_ANALYSIS_ztau3mutaue_skimmer_default.root","READ");
    
    TH2D * tauh_OS_vs_3mu  = (TH2D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_OS_vs_3mu_triggerMC4");
    TH2D * taumu_OS_vs_3mu  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_skimmer_default_OS_vs_3mu_triggerMC3");
    TH2D * taue_OS_vs_3mu  = (TH2D*)TreeFile_e->Get("ztau3mutaue_skimmer_default_OS_vs_3mu_triggerMC2");
    
    std::cout << "For tau_h: "<<std::endl;
    std::cout << "Percent of events where HLT fails: " << (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2))   /    (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2)+tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))
    << " of these, matched to trigger object " << (tauh_OS_vs_3mu->GetBinContent(1,2)) / (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2))
    << " overall os matched: " << (tauh_OS_vs_3mu->GetBinContent(1,2))   /    (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2)+tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))
    << " and percent of HLT pass: " << (tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))   /    (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2)+tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "For tau_mu: "<<std::endl;
    std::cout << "Percent of events where HLT fails: " << (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2))   /    (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2)+taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))
    << " of these, matched to trigger object " << (taumu_OS_vs_3mu->GetBinContent(1,2)) / (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2))
    << " overall os matched: " << (taumu_OS_vs_3mu->GetBinContent(1,2))   /    (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2)+taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))
    << " and percent of HLT pass: " << (taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))   /    (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2)+taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "For tau_e: "<<std::endl;
    std::cout << "Percent of events where HLT fails: " << (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2))   /    (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2)+taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))
    << " of these, matched to trigger object " << (taue_OS_vs_3mu->GetBinContent(1,2)) / (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2))
    << " overall os matched: " << (taue_OS_vs_3mu->GetBinContent(1,2))   /    (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2)+taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))
    << " and percent of HLT pass: " << (taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))   /    (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2)+taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))
    << std::endl;
    
    
    
    
}