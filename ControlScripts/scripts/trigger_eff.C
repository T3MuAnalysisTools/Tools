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

void trigger_eff () 
{
    TFile *TreeFile_h = new TFile("LOCAL_COMBINED_ztau3mutauh_prefc_default.root","READ");
    TFile *TreeFile_mu = new TFile("LOCAL_COMBINED_ztau3mutaumu_prefc_default.root","READ");
    TFile *TreeFile_e = new TFile("LOCAL_COMBINED_ztau3mutaue_prefc_default.root","READ");
    
    TH2D * tauh_OS_vs_3mu  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_OS_vs_3mu_triggerMC4");
    TH2D * taumu_OS_vs_3mu  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_OS_vs_3mu_triggerMC3");
    TH2D * taue_OS_vs_3mu  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_OS_vs_3mu_triggerMC2");
    
    TH1D * tauh_Whether_4object_reconstructed  = (TH1D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_4object_reconstructedMC4");
    TH2D * tauh_Whether_4object_fiducial_reco  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_4object_fiducial_and_reconstructedMC4");
    TH2D * tauh_Whether_triplet_fiducial_or_reco  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_triplet_fiducial_or_recoMC4");
    TH2D * tauh_Whether_reco_triplet_trigger_L1_and_HLT  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_reco_triplet_trigger_L1_and_HLTMC4");
    TH2D * tauh_Whether_Mu1_fiducial_and_reco  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_Mu1_fiducial_and_recoMC4");
    TH2D * tauh_Whether_Mu2_fiducial_and_reco  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_Mu2_fiducial_and_recoMC4");
    TH2D * tauh_Whether_Mu3_fiducial_and_reco  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_Mu3_fiducial_and_recoMC4");
    TH2D * tauh_Whether_Tau_h_fiducial_and_reco  = (TH2D*)TreeFile_h->Get("ztau3mutauh_prefc_default_Whether_Tau_h_fiducial_and_recoMC4");
    
    TH1D * taumu_Whether_4object_reconstructed  = (TH1D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_4object_reconstructedMC3");
    TH2D * taumu_Whether_4object_fiducial_reco  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_4object_fiducial_and_reconstructedMC3");
    TH2D * taumu_Whether_triplet_fiducial_or_reco  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_triplet_fiducial_or_recoMC3");
    TH2D * taumu_Whether_reco_triplet_trigger_L1_and_HLT  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_reco_triplet_trigger_L1_and_HLTMC3");
    TH2D * taumu_Whether_Mu1_fiducial_and_reco  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_Mu1_fiducial_and_recoMC3");
    TH2D * taumu_Whether_Mu2_fiducial_and_reco  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_Mu2_fiducial_and_recoMC3");
    TH2D * taumu_Whether_Mu3_fiducial_and_reco  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_Mu3_fiducial_and_recoMC3");
    TH2D * taumu_Whether_Tau_mu_fiducial_and_reco  = (TH2D*)TreeFile_mu->Get("ztau3mutaumu_prefc_default_Whether_Tau_mu_fiducial_and_recoMC3");
    
    TH1D * taue_Whether_4object_reconstructed  = (TH1D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_4object_reconstructedMC2");
    TH2D * taue_Whether_4object_fiducial_reco  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_4object_fiducial_and_reconstructedMC2");
    TH2D * taue_Whether_triplet_fiducial_or_reco  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_triplet_fiducial_or_recoMC2");
    TH2D * taue_Whether_reco_triplet_trigger_L1_and_HLT  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_reco_triplet_trigger_L1_and_HLTMC2");
    TH2D * taue_Whether_Mu1_fiducial_and_reco  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_Mu1_fiducial_and_recoMC2");
    TH2D * taue_Whether_Mu2_fiducial_and_reco  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_Mu2_fiducial_and_recoMC2");
    TH2D * taue_Whether_Mu3_fiducial_and_reco  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_Mu3_fiducial_and_recoMC2");
    TH2D * taue_Whether_Tau_e_fiducial_and_reco  = (TH2D*)TreeFile_e->Get("ztau3mutaue_prefc_default_Whether_Tau_e_fiducial_and_recoMC2");
    
    
    //For tau h
    std::cout << "For tau_h: "<<std::endl;
    std::cout << "Percent of events where HLT fails: " << (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2))   /    (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2)+tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))
    << " of these, matched to trigger object " << (tauh_OS_vs_3mu->GetBinContent(1,2)) / (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2))
    << " overall os matched: " << (tauh_OS_vs_3mu->GetBinContent(1,2))   /    (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2)+tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))
    << " and percent of HLT pass: " << (tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))   /    (tauh_OS_vs_3mu->GetBinContent(1,1)+tauh_OS_vs_3mu->GetBinContent(1,2)+tauh_OS_vs_3mu->GetBinContent(2,1)+tauh_OS_vs_3mu->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of mu1 in the fiducial volume: " << (tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))   /    (tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(1,1)+tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(1,2)+tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2) / ( tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1) + tauh_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of mu2 in the fiducial volume: " << (tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2))   /    (tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(1,1)+tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(1,2)+tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2) / ( tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1) + tauh_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2) )
    << std::endl;
    
    std::cout << "Percent of mu3 in the fiducial volume: " << (tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2))   /    (tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(1,1)+tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(1,2)+tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2) / ( tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1) + tauh_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2) )
    << std::endl;
    
    std::cout << "Percent of tauh in the fiducial volume: " << (tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,2))   /    (tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(1,1)+tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(1,2)+tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,1)+tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,2) / ( tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,1) + tauh_Whether_Tau_h_fiducial_and_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of events with 4 objects in fiducial phase space: " << (tauh_Whether_4object_fiducial_reco->GetBinContent(2,1)+tauh_Whether_4object_fiducial_reco->GetBinContent(2,2))   /    (tauh_Whether_4object_fiducial_reco->GetBinContent(1,1)+tauh_Whether_4object_fiducial_reco->GetBinContent(1,2)+tauh_Whether_4object_fiducial_reco->GetBinContent(2,1)+tauh_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << " of these percent (of 4 objects) reconstructed: " << tauh_Whether_4object_fiducial_reco->GetBinContent(2,2) / ( tauh_Whether_4object_fiducial_reco->GetBinContent(2,1) + tauh_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << " overall: " << tauh_Whether_4object_fiducial_reco->GetBinContent(2,2) /  (tauh_Whether_4object_fiducial_reco->GetBinContent(1,1)+tauh_Whether_4object_fiducial_reco->GetBinContent(1,2)+tauh_Whether_4object_fiducial_reco->GetBinContent(2,1)+tauh_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of events with 4 objects reconstructed (no fiducial cuts): " << tauh_Whether_4object_reconstructed->GetBinContent(2) / ( tauh_Whether_4object_reconstructed->GetBinContent(2) + tauh_Whether_4object_reconstructed->GetBinContent(1) )
    << std::endl;
    
    std::cout << "Percent of triplets in fiducial volume with cuts: " << (tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (tauh_Whether_triplet_fiducial_or_reco->GetBinContent(1,1)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(1,2)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << " of these, percent reconstructed: " << (tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << " with overall reco (percent of triplets in fiducial volume and reconstructed): " << (tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (tauh_Whether_triplet_fiducial_or_reco->GetBinContent(1,1)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(1,2)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+tauh_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "Percent of reco triplets Triggering L1T: " << (tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,1)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,2)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << " of these, percent triggering HLT: " << (tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << " Trigger Eff (wrt reco): " << (tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,1)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,2)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+tauh_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    std::cout << "Sanity check: "<< tauh_Whether_4object_reconstructed->GetBinContent(2) << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    
    
    
    
    //For tau mu
    
    std::cout << "For tau_mu: "<<std::endl;
    std::cout << "Percent of events where HLT fails: " << (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2))   /    (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2)+taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))
    << " of these, matched to trigger object " << (taumu_OS_vs_3mu->GetBinContent(1,2)) / (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2))
    << " overall os matched: " << (taumu_OS_vs_3mu->GetBinContent(1,2))   /    (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2)+taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))
    << " and percent of HLT pass: " << (taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))   /    (taumu_OS_vs_3mu->GetBinContent(1,1)+taumu_OS_vs_3mu->GetBinContent(1,2)+taumu_OS_vs_3mu->GetBinContent(2,1)+taumu_OS_vs_3mu->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "Percent of mu1 in the fiducial volume: " << (taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))   /    (taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(1,1)+taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(1,2)+taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2) / ( taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1) + taumu_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of mu2 in the fiducial volume: " << (taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2))   /    (taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(1,1)+taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(1,2)+taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2) / ( taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1) + taumu_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2) )
    << std::endl;
    
    std::cout << "Percent of mu3 in the fiducial volume: " << (taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2))   /    (taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(1,1)+taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(1,2)+taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2) / ( taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1) + taumu_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2) )
    << std::endl;
    
    std::cout << "Percent of taumu in the fiducial volume: " << (taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,2))   /    (taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(1,1)+taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(1,2)+taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,1)+taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,2) / ( taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,1) + taumu_Whether_Tau_mu_fiducial_and_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of events with 4 objects in fiducial phase space: " << (taumu_Whether_4object_fiducial_reco->GetBinContent(2,1)+taumu_Whether_4object_fiducial_reco->GetBinContent(2,2))   /    (taumu_Whether_4object_fiducial_reco->GetBinContent(1,1)+taumu_Whether_4object_fiducial_reco->GetBinContent(1,2)+taumu_Whether_4object_fiducial_reco->GetBinContent(2,1)+taumu_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << " of these percent (of 4 objects) reconstructed: " << taumu_Whether_4object_fiducial_reco->GetBinContent(2,2) / ( taumu_Whether_4object_fiducial_reco->GetBinContent(2,1) + taumu_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << " overall: " << taumu_Whether_4object_fiducial_reco->GetBinContent(2,2) /  (taumu_Whether_4object_fiducial_reco->GetBinContent(1,1)+taumu_Whether_4object_fiducial_reco->GetBinContent(1,2)+taumu_Whether_4object_fiducial_reco->GetBinContent(2,1)+taumu_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of events with 4 objects reconstructed (no fiducial cuts): " << taumu_Whether_4object_reconstructed->GetBinContent(2) / ( taumu_Whether_4object_reconstructed->GetBinContent(2) + taumu_Whether_4object_reconstructed->GetBinContent(1) )
    << std::endl;
    
    std::cout << "Percent of triplets in fiducial volume with cuts: " << (taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (taumu_Whether_triplet_fiducial_or_reco->GetBinContent(1,1)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(1,2)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << " of these, percent reconstructed: " << (taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << " with overall reco (percent of triplets in fiducial volume and reconstructed): " << (taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (taumu_Whether_triplet_fiducial_or_reco->GetBinContent(1,1)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(1,2)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taumu_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "Percent of reco triplets Triggering L1T: " << (taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,1)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,2)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << " of these, percent triggering HLT: " << (taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << " Trigger Eff (wrt reco): " << (taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,1)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,2)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taumu_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    
    
    
    
    
    //For tau e
    
    std::cout << "For tau_e: "<<std::endl;
    std::cout << "Percent of events where HLT fails: " << (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2))   /    (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2)+taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))
    << " of these, matched to trigger object " << (taue_OS_vs_3mu->GetBinContent(1,2)) / (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2))
    << " overall os matched: " << (taue_OS_vs_3mu->GetBinContent(1,2))   /    (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2)+taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))
    << " and percent of HLT pass: " << (taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))   /    (taue_OS_vs_3mu->GetBinContent(1,1)+taue_OS_vs_3mu->GetBinContent(1,2)+taue_OS_vs_3mu->GetBinContent(2,1)+taue_OS_vs_3mu->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "Percent of mu1 in the fiducial volume: " << (taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))   /    (taue_Whether_Mu1_fiducial_and_reco->GetBinContent(1,1)+taue_Whether_Mu1_fiducial_and_reco->GetBinContent(1,2)+taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2) / ( taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,1) + taue_Whether_Mu1_fiducial_and_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of mu2 in the fiducial volume: " << (taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2))   /    (taue_Whether_Mu2_fiducial_and_reco->GetBinContent(1,1)+taue_Whether_Mu2_fiducial_and_reco->GetBinContent(1,2)+taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2) / ( taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,1) + taue_Whether_Mu2_fiducial_and_reco->GetBinContent(2,2) )
    << std::endl;
    
    std::cout << "Percent of mu3 in the fiducial volume: " << (taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2))   /    (taue_Whether_Mu3_fiducial_and_reco->GetBinContent(1,1)+taue_Whether_Mu3_fiducial_and_reco->GetBinContent(1,2)+taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2) / ( taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,1) + taue_Whether_Mu3_fiducial_and_reco->GetBinContent(2,2) )
    << std::endl;
    
    std::cout << "Percent of taue in the fiducial volume: " << (taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,2))   /    (taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(1,1)+taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(1,2)+taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,1)+taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,2))
    << " of these percent reconstructed: " << taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,2) / ( taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,1) + taue_Whether_Tau_e_fiducial_and_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of events with 4 objects in fiducial phase space: " << (taue_Whether_4object_fiducial_reco->GetBinContent(2,1)+taue_Whether_4object_fiducial_reco->GetBinContent(2,2))   /    (taue_Whether_4object_fiducial_reco->GetBinContent(1,1)+taue_Whether_4object_fiducial_reco->GetBinContent(1,2)+taue_Whether_4object_fiducial_reco->GetBinContent(2,1)+taue_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << " of these percent (of 4 objects) reconstructed: " << taue_Whether_4object_fiducial_reco->GetBinContent(2,2) / ( taue_Whether_4object_fiducial_reco->GetBinContent(2,1) + taue_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << " overall: " << taue_Whether_4object_fiducial_reco->GetBinContent(2,2) /  (taue_Whether_4object_fiducial_reco->GetBinContent(1,1)+taue_Whether_4object_fiducial_reco->GetBinContent(1,2)+taue_Whether_4object_fiducial_reco->GetBinContent(2,1)+taue_Whether_4object_fiducial_reco->GetBinContent(2,2))
    << std::endl;
    
    std::cout << "Percent of events with 4 objects reconstructed (no fiducial cuts): " << taue_Whether_4object_reconstructed->GetBinContent(2) / ( taue_Whether_4object_reconstructed->GetBinContent(2) + taue_Whether_4object_reconstructed->GetBinContent(1) )
    << std::endl;
    
    std::cout << "Percent of triplets in fiducial volume with cuts: " << (taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (taue_Whether_triplet_fiducial_or_reco->GetBinContent(1,1)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(1,2)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << " of these, percent reconstructed: " << (taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << " with overall reco (percent of triplets in fiducial volume and reconstructed): " << (taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))   /    (taue_Whether_triplet_fiducial_or_reco->GetBinContent(1,1)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(1,2)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,1)+taue_Whether_triplet_fiducial_or_reco->GetBinContent(2,2))
    << std::endl;
    
    
    std::cout << "Percent of reco triplets Triggering L1T: " << (taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,1)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,2)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << " of these, percent triggering HLT: " << (taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << " Trigger Eff (wrt reco): " << (taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))   /    (taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,1)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(1,2)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,1)+taue_Whether_reco_triplet_trigger_L1_and_HLT->GetBinContent(2,2))
    << std::endl;
    
    
    
    
}