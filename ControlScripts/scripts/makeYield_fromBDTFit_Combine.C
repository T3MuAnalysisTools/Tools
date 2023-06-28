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

void makeYield_fromBDTFit_Combine () 
{
    system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv");
    
    //Category names: 0=tauh, 1=taumu, 2=taue
    TString cat_name[3];
    cat_name[0] = "ztau3mutauh_default_";
    cat_name[1] = "ztau3mutaumu_default_";
    cat_name[2] = "ztau3mutaue_default_";
    
    TString cat_label[3];
    cat_label[0] = "{h}";
    cat_label[1] = "{#mu}";
    cat_label[2] = "{e}";
    
    TString print_label[3];
    print_label[0] = "h";
    print_label[1] = "mu";
    print_label[2] = "e";
    
    TString card_modifier_name[3];
    card_modifier_name[0] = "ZTT_tauh_test";
    card_modifier_name[1] = "ZTT_taumu_test";
    card_modifier_name[2] = "ZTT_taue_test";
    
    TString combined_card_name[3];
    combined_card_name[0] = "ZTT_tauh_Combined_Mod";
    combined_card_name[1] = "ZTT_taumu_Combined_Mod";
    combined_card_name[2] = "ZTT_taue_Combined_Mod";
    
    TString hname;
    
    float signal_region_min(1.4);
    float signal_region_max(2.1);
    
    float signal_peak_region_min(1.73);
    float signal_peak_region_max(1.82);
    
    //Filename and histograms
    TFile * file_tau[3];
    
    TH1D  * tau_T3Mu[3];
    TH1D  * tau_T3Mu_Dat[3];
    TH1D  * tau_BDT_Output_Data[3];
    TH1D  * tau_BDT_Output_MC[3];
    TH2D  * tau_T3Mu_vs_BDT_Output_Data[3];
    TH1D  * tau_T3Mu_vs_BDT_Output_Data_Projection[3];
    
    TH2D  * tau_cut1_vs_cut2_vs_sig[3];
    TH2D  * tau_cut1_vs_cut2_vs_limit[3];
    
    
    TFile *TreeFile = new TFile("Combine_Tree_ztau3mutau.root","READ");
    TTree *tree[3];
    
    Float_t tripletMass;
    Float_t bdt_cv;
    Float_t weight;
    Float_t isMC;
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      tree[i] = (TTree *) TreeFile->Get(cat_base[i]);
      
      tau_T3Mu[i] = new TH1D("tau_T3Mu","tau_T3Mu_"+hname,40,signal_region_min,signal_region_max);
      tau_T3Mu_Dat[i] = new TH1D("tau_T3Mu_Dat","tau_T3Mu_Dat_"+hname,40,signal_region_min,signal_region_max);
      tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,100,-0.9,0.9);
      tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,100,-0.9,0.9);
      
      tree[i]->SetBranchAddress("tripletMass",&tripletMass);
      tree[i]->SetBranchAddress("bdt_cv",&bdt_cv);
      tree[i]->SetBranchAddress("weight",&weight);
      tree[i]->SetBranchAddress("isMC",&isMC);
      
      Long64_t nentries = tree[i]->GetEntries();
      for (Long64_t j=0;j<nentries;j++) {
        tree[i]->GetEntry(j);
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max){
                if(isMC>0){
                  tau_T3Mu[i]->Fill(tripletMass,weight);
                  tau_BDT_Output_MC[i]->Fill(bdt_cv,weight);
                }
                if(isMC==0 && (tripletMass<=signal_region_min || tripletMass>=signal_peak_region_max) ){//blinded
                  tau_T3Mu_Dat[i]->Fill(tripletMass);
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      //tau_T3Mu_vs_BDT_Output_Data[i] = (TH2D*)file_tau[i]->Get(cat_name[i]+"BDT_2Dscan_TripletMassData");
      
      
    }
    
    
    
    /*
    //For fitting BDT Output: Using crystal ball
    //General
    RooRealVar * BDTOutput_x[3];
    RooRealVar * cbmean[3];
    RooRealVar * cbsigma[3];
    RooRealVar * n[3];
    RooRealVar * alpha[3];
    
    RooCBShape * cball[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, #tau_"+cat_label[i],-0.9,0.9);
      BDTOutput_x[i]->setRange("R1",-0.9,0.9);
      cbmean[i] = new RooRealVar("cbmean"+hname, "cbmean" , -0.5, -0.9,0.0) ;
      cbsigma[i] = new RooRealVar("cbsigma"+hname, "cbsigma" , 0.2, 0.000001, 1.0) ;
      n[i] = new RooRealVar("n"+hname, "n", 15.0, 0.5, 20);
      alpha[i] = new RooRealVar("alpha"+hname,"alpha value CB",5.0,0.01,10);
      
      cball[i] = new RooCBShape("cball"+hname, "crystal ball", *BDTOutput_x[i], *cbmean[i], *cbsigma[i], *alpha[i], *n[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*cball[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("R1"), Save());
      
    }
    */
    
    
    
    
    
    
    //For fitting BDT Output in Data: Using Bifurgauss
    //General
    RooRealVar * BDTOutput_x[3];
    RooRealVar * bgausmean[3];
    RooRealVar * bgaussigma_a[3];
    RooRealVar * bgaussigma_b[3];
    
    
    RooBifurGauss * bgaus_dist[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, #tau_"+cat_label[i],-0.9,0.9);
      BDTOutput_x[i]->setRange("R1",-0.9,0.9);
      bgausmean[i] = new RooRealVar("bgausmean"+hname, "bgausmean" , -0.5, -0.9,0.0) ;
      bgaussigma_a[i] = new RooRealVar("bgaussigma_a"+hname, "bgaussigma_a" , 0.2, 0.000001, 1.0) ;
      bgaussigma_b[i] = new RooRealVar("bgaussigma_b"+hname, "bgaussigma_b" , 0.2, 0.000001, 1.0);
      
      bgaus_dist[i] = new RooBifurGauss("bgaus_dist"+hname, "bgaus dist", *BDTOutput_x[i], *bgausmean[i], *bgaussigma_a[i], *bgaussigma_b[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*bgaus_dist[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("R1"), Save());
      
    }
    
    
    
    /*
    //BDT Fit Plots
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1800, 600);
    canvas2->Divide(3, 1);
    
    RooPlot * xFrame_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas2->cd( i+1 );
      xFrame_bdt[i] = BDTOutput_x[i]->frame();
      bdt_data[i]->plotOn(xFrame_bdt[i]);
      BDT_distribution[i]->plotOn(xFrame_bdt[i],LineColor(4),LineWidth(2), Normalization(bdt_data[i]->sumEntries("1", "R1"), RooAbsReal::NumEvent),ProjectionRange("R1"));
      xFrame_bdt[i]->SetTitle("BDT Output, #tau_"+cat_label[i]);
      xFrame_bdt[i]->SetXTitle("BDT Score");
      xFrame_bdt[i]->SetYTitle("Events");
      xFrame_bdt[i]->Draw();
    }
    */
    
    
    
    
    
    
    
    
    //For fitting BDT Output in MC: Using Bifurgauss
    //General
    RooRealVar * bgausmeanMC[3];
    RooRealVar * bgaussigmaMC_a[3];
    RooRealVar * bgaussigmaMC_b[3];
    
    
    RooBifurGauss * bgaus_distMC[3];
    RooDataHist * bdt_MC[3];
    RooRealVar * BDTNormMC[3];
    RooAddPdf * BDT_distributionMC[3];
    RooFitResult * fitresult_bdtMC[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      bgausmeanMC[i] = new RooRealVar("bgausmeanMC"+hname, "bgausmeanMC" , 0.5, 0.0,0.9) ;
      bgaussigmaMC_a[i] = new RooRealVar("bgaussigmaMC_a"+hname, "bgaussigmaMC_a" , 0.2, 0.000001, 1.0) ;
      bgaussigmaMC_b[i] = new RooRealVar("bgaussigmaMC_b"+hname, "bgaussigmaMC_b" , 0.2, 0.000001, 1.0);
      
      bgaus_distMC[i] = new RooBifurGauss("bgaus_distMC"+hname, "bgaus dist MC", *BDTOutput_x[i], *bgausmeanMC[i], *bgaussigmaMC_a[i], *bgaussigmaMC_b[i]);
      bdt_MC[i] = new RooDataHist("bdt_MC"+hname, "bdt_MC", *BDTOutput_x[i], Import(*tau_BDT_Output_MC[i]));
      BDTNormMC[i] = new RooRealVar("BDTNormMC"+hname, "BDTNormMC", 5.0,0.0,50);
      BDT_distributionMC[i] = new RooAddPdf("BDT_distributionMC"+hname, "BDT_distributionMC", RooArgList(*bgaus_distMC[i]), RooArgList(*BDTNormMC[i]));
      fitresult_bdtMC[i] = BDT_distributionMC[i]->fitTo(*bdt_MC[i], Range("R1"), Save());
      
    }
    
    
    
    /*
    //MC BDT Fit Plots
    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 1800, 600);
    canvas3->Divide(3, 1);
    
    RooPlot * xFrame_bdtMC[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas3->cd( i+1 );
      xFrame_bdtMC[i] = BDTOutput_x[i]->frame();
      bdt_MC[i]->plotOn(xFrame_bdtMC[i]);
      BDT_distributionMC[i]->plotOn(xFrame_bdtMC[i],LineColor(4),LineWidth(2), Normalization(bdt_MC[i]->sumEntries("1", "R1"), RooAbsReal::NumEvent),ProjectionRange("R1"));
      xFrame_bdtMC[i]->SetTitle("BDT Output, #tau_"+cat_label[i]);
      xFrame_bdtMC[i]->SetXTitle("BDT Score");
      xFrame_bdtMC[i]->SetYTitle("Events");
      xFrame_bdtMC[i]->Draw();
    }
    */
    
    
    
    
    
    
    

    //cout << "tauh MC count: " << tauh_Vis_Mass->Integral() << " tauh Data count: " << tauh_Vis_Mass_Dat->Integral() << endl;
    //cout << "taumu MC count: " << taumu_Vis_Mass->Integral() << " taumu Data count: " << taumu_Vis_Mass_Dat->Integral() << endl;
    //cout << "taue MC count: " << taue_Vis_Mass->Integral() << " taue Data count: " << taue_Vis_Mass_Dat->Integral() << endl;
    
    
    
    
    
    //Triplet Mass Fits
    RooRealVar * InvMass[3];
    
    RooPolynomial * poly[3];
    RooDataHist * data[3];
    RooRealVar * LineNorm[3];
    RooAddPdf * pdf[3];
    RooFitResult * fitresult[3];
    
    RooRealVar * mean[3];
    RooRealVar * sigma[3];
    RooGaussian * Gauss[3];
    RooDataHist * mc[3];
    RooRealVar * GaussNorm[3];
    RooAddPdf * mc_pdf[3];
    RooFitResult * mc_fitresult[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      InvMass[i] = new RooRealVar("InvMass"+hname,"InvMass, #tau_"+cat_label[i],signal_region_min,signal_region_max);
      InvMass[i]->setRange("R1",signal_region_min,signal_peak_region_min); //background   
      InvMass[i]->setRange("R2",signal_peak_region_max,signal_region_max); //background
      InvMass[i]->setRange("R3",1.70,1.85); //signal range for fitting
      InvMass[i]->setRange("R4",signal_peak_region_min,signal_peak_region_max); //signal range for yield
      
      
      //Flat fit for data
      poly[i] = new RooPolynomial("poly"+hname, "poly dist", *InvMass[i]);
      data[i] = new RooDataHist("data"+hname, "data", *InvMass[i], Import(*tau_T3Mu_Dat[i]));
      LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,15);
      pdf[i] = new RooAddPdf("pdf"+hname, "pdf", RooArgList(*poly[i]), RooArgList(*LineNorm[i]));
      fitresult[i] = pdf[i]->fitTo(*data[i], Range("R1,R2"), Save());
      
      
      //Gaussian fit for MC
      mean[i] = new RooRealVar("mean"+hname, "mean" , 1.776,0.,5.) ;
      sigma[i] = new RooRealVar("sigma"+hname, "sigma" , 0.5,0.001,10) ;
      
      Gauss[i] = new RooGaussian("Gauss"+hname, "Gauss dist", *InvMass[i], *mean[i], *sigma[i]);
      mc[i] = new RooDataHist("mc"+hname, "mc", *InvMass[i], Import(*tau_T3Mu[i]));
      GaussNorm[i] = new RooRealVar("GaussNorm"+hname, "GaussNorm",  0.5,0.001,1.0);
      mc_pdf[i] = new RooAddPdf("mc_pdf"+hname, "mc_pdf", RooArgList(*Gauss[i]), RooArgList(*GaussNorm[i]));
      fitresult_bdt[i] = mc_pdf[i]->fitTo(*mc[i], Range("R3"), Save());
      
    }
    
    
    
    
    
    /*
    // This gives the integral from the fits.
    double pdf1_integral_restricted = pdf1.createIntegral(InvMass1,NormSet(InvMass1),Range("R4"))->getVal();
    double mc_pdf1_integral_restricted = mc_pdf1.createIntegral(InvMass1,NormSet(InvMass1),Range("R4"))->getVal();
    double pdf2_integral_restricted = pdf2.createIntegral(InvMass2,NormSet(InvMass2),Range("R4"))->getVal();
    double mc_pdf2_integral_restricted = mc_pdf2.createIntegral(InvMass2,NormSet(InvMass2),Range("R4"))->getVal();
    double pdf3_integral_restricted = pdf3.createIntegral(InvMass3,NormSet(InvMass3),Range("R4"))->getVal();
    double mc_pdf3_integral_restricted = mc_pdf3.createIntegral(InvMass3,NormSet(InvMass3),Range("R4"))->getVal();
    
    // Normalizations need to be added manually. The pdfs are normalized to 1 and scaled to the data plotted. nData1 and nSignal1 are for normalization (same region as fit). "R4" is for yields.
    const double nData1 = data1.sumEntries("1", "R1,R2");
    const double nSignal1 = mc1.sumEntries("1", "R3");
    const double nSignal1_restricted = mc1.sumEntries("1", "R4");// used a separate range for getting the yield and a different range for fitting
    const double nData2 = data2.sumEntries("1", "R1,R2");
    const double nSignal2 = mc2.sumEntries("1", "R3");
    const double nSignal2_restricted = mc2.sumEntries("1", "R4");
    const double nData3 = data3.sumEntries("1", "R1,R2");
    const double nSignal3 = mc3.sumEntries("1", "R3");
    const double nSignal3_restricted = mc3.sumEntries("1", "R4");
    
    cout << "nData1: " << nData1 << " nSignal1: " << nSignal1_restricted << endl;
    cout << "nData2: " << nData2 << " nSignal2: " << nSignal2_restricted << endl;
    cout << "nData3: " << nData3 << " nSignal3: " << nSignal3_restricted << endl;
    
    cout << " LineNorm1: " << LineNorm1.getValV() << " guessed data in signal region: " << (pdf1_integral_restricted/(1-pdf1_integral_restricted))*nData1 << endl;
    cout << " LineNorm1: " << LineNorm2.getValV() << " guessed data in signal region: " << (pdf2_integral_restricted/(1-pdf2_integral_restricted))*nData2 << endl;
    cout << " LineNorm1: " << LineNorm3.getValV() << " guessed data in signal region: " << (pdf3_integral_restricted/(1-pdf3_integral_restricted))*nData3 << endl;
    
    cout << "mc_pdf1_integral: " << GaussNorm1.getValV()*mc_pdf1_integral_restricted << " pdf1_integral: " << LineNorm1.getValV()*pdf1_integral_restricted << endl;
    cout << "mc_pdf2_integral: " << GaussNorm2.getValV()*mc_pdf2_integral_restricted << " pdf2_integral: " << LineNorm2.getValV()*pdf2_integral_restricted << endl;
    cout << "mc_pdf3_integral: " << GaussNorm3.getValV()*mc_pdf3_integral_restricted << " pdf3_integral: " << LineNorm3.getValV()*pdf3_integral_restricted << endl;
    
    double scaling1 = mc1.sumEntries("1")/tauh_T3Mu->GetEntries();
    double scaling2 = mc2.sumEntries("1")/taumu_T3Mu->GetEntries();
    double scaling3 = mc3.sumEntries("1")/taue_T3Mu->GetEntries();
    
    cout << "Unscaled mc1: " << nSignal1_restricted/scaling1 << " scaling 1: " << scaling1 << endl;// Getting the unweighted content of the signal histogram
    cout << "Unscaled mc2: " << nSignal2_restricted/scaling2 << " scaling 2: " << scaling2 << endl;
    cout << "Unscaled mc3: " << nSignal3_restricted/scaling3 << " scaling 3: " << scaling3 << endl;
    */
    
    
    double pdf_integral_restricted[3];
    double mc_pdf_integral_restricted[3];
    double nData[3];
    double nSignal[3];
    double nSignal_restricted[3];
    double scaling[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      // This gives the integral from the fits.
      pdf_integral_restricted[i] = pdf[i]->createIntegral(*InvMass[i],NormSet(*InvMass[i]),Range("R4"))->getVal();
      mc_pdf_integral_restricted[i] = mc_pdf[i]->createIntegral(*InvMass[i],NormSet(*InvMass[i]),Range("R4"))->getVal();
      
      // Normalizations need to be added manually. The pdfs are normalized to 1 and scaled to the data plotted. nData1 and nSignal1 are for normalization (same region as fit). "R4" is for yields.
      nData[i] = data[i]->sumEntries("1", "R1,R2");
      nSignal[i] = mc[i]->sumEntries("1", "R3");
      nSignal_restricted[i] = mc[i]->sumEntries("1", "R4");// used a separate range for getting the yield and a different range for fitting
      
      
      
      cout << "  " << endl;
      
      cout << "nData "+print_label[i]+" : " << nData[i] << " nSignal "+print_label[i]+" : " << nSignal_restricted[i] << endl;
      cout << " LineNorm "+print_label[i]+" : " << LineNorm[i]->getValV() << " guessed data in signal region: " << (pdf_integral_restricted[i]/(1-pdf_integral_restricted[i]))*nData[i] << endl;
      cout << "mc_pdf "+print_label[i]+" _integral: " << GaussNorm[i]->getValV()*mc_pdf_integral_restricted[i] << " pdf "+print_label[i]+" _integral: " << LineNorm[i]->getValV()*pdf_integral_restricted[i] << endl;
      
      scaling[i] = mc[i]->sumEntries("1")/(tau_T3Mu[i]->GetEntries());
      
      cout << "Unscaled mc "+print_label[i]+" : " << nSignal_restricted[i]/scaling[i] << " scaling  "+print_label[i]+" : " << scaling[i] << endl;// Getting the unweighted content of the signal histogram
      
      cout << "  " << endl;
      
    }
    
    
    
    
    /*
    //Triplet Mass Fit Plots
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 1800, 600);
    canvas1->Divide(3, 1);
    
    RooPlot * xFrame[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas1->cd( i+1 );
      xFrame[i] = InvMass[i]->frame();
      data[i]->plotOn(xFrame[i]);
      pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), #tau_"+cat_label[i]);
      xFrame[i]->SetXTitle("3#mu inv. mass (GeV)");
      xFrame[i]->SetYTitle("Events");
      xFrame[i]->Draw();
    }
    */
    
    
    
    
    
    
    
    /*
    
    
    double TightBDTCutFraction1 = BDT_distribution1.createIntegral(BDTOutput_x1,NormSet(BDTOutput_x1),Range("R_Small"))->getVal() / BDT_distribution1.createIntegral(BDTOutput_x1,NormSet(BDTOutput_x1),Range("R_Large"))->getVal();
    double TightBDTCutFraction2 = BDT_distribution2.createIntegral(BDTOutput_x2,NormSet(BDTOutput_x2),Range("R_Small"))->getVal() / BDT_distribution2.createIntegral(BDTOutput_x2,NormSet(BDTOutput_x2),Range("R_Large"))->getVal();
    double TightBDTCutFraction3 = BDT_distribution3.createIntegral(BDTOutput_x3,NormSet(BDTOutput_x3),Range("R_Small"))->getVal() / BDT_distribution3.createIntegral(BDTOutput_x3,NormSet(BDTOutput_x3),Range("R_Large"))->getVal();
    
    cout << "Yield mc1: " << BDT_distributionMC1.createIntegral(BDTOutput_x1,NormSet(BDTOutput_x1),Range("R_Small"))->getVal() * BDTNormMC1.getValV() << " Yield Background 1: " << (pdf1_integral_restricted/(1-pdf1_integral_restricted))*nData1*TightBDTCutFraction1 << endl;
    cout << "Yield mc2: " << BDT_distributionMC2.createIntegral(BDTOutput_x2,NormSet(BDTOutput_x2),Range("R_Small"))->getVal() * BDTNormMC2.getValV() << " Yield Background 2: " << (pdf2_integral_restricted/(1-pdf2_integral_restricted))*nData2*TightBDTCutFraction2 << endl;
    cout << "Yield mc3: " << BDT_distributionMC3.createIntegral(BDTOutput_x3,NormSet(BDTOutput_x3),Range("R_Small"))->getVal() * BDTNormMC3.getValV() << " Yield Background 3: " << (pdf3_integral_restricted/(1-pdf3_integral_restricted))*nData3*TightBDTCutFraction3 << endl;
    */
    
    
    
    
    
    
    
    //From BDT Output:
    BDTOutput_x[0]->setRange("R_Small",0.309056,100); BDTOutput_x[0]->setRange("R_Large",0.05,100);
    BDTOutput_x[1]->setRange("R_Small",0.329983,100); BDTOutput_x[1]->setRange("R_Large",0.05,100);
    BDTOutput_x[2]->setRange("R_Small",0.333186,100); BDTOutput_x[2]->setRange("R_Large",0.05,100);
    
    double TightBDTCutFraction[3];
    
    RooAddPdf BDT_distributionTest = *BDT_distribution[0];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      TightBDTCutFraction[i] = BDT_distribution[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Small"))->getVal() / BDT_distribution[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Large"))->getVal();
      cout << "Yield mc "+print_label[i]+" : " << BDT_distributionMC[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Small"))->getVal() * (bdt_MC[i]->sumEntries("1", "R1")) << " Yield Background "+print_label[i]+" : " << (pdf_integral_restricted[i]/(1-pdf_integral_restricted[i]))*nData[i]*TightBDTCutFraction[i] << endl;
      
    }
    
    
    
    
    
    
        //Calculate limits
        double a[3],b[3];
        double N_s_1[3], N_s_2[3];
        double N_b_1[3], N_b_2[3];
        double S1[3], S2[3], S[3];
        std::vector<double> S1_list[3], S2_list[3], S_list[3], a_list[3], b_list[3], N_s_1_yield_list[3], N_s_2_yield_list[3], N_b_1_yield_list[3], N_b_2_yield_list[3];
        
        TString command_a[3];
        TString command_b[3];
        TString command_a_and_b[3];
        TString command_run[3];
        
        //double X_min = std::min(tau_BDT_Output_Data[0]->GetXaxis()->GetXmin(), tau_BDT_Output_MC[0]->GetXaxis()->GetXmin());
        //double X_max = std::max(tau_BDT_Output_Data[0]->GetXaxis()->GetXmax(), tau_BDT_Output_MC[0]->GetXaxis()->GetXmax());
        
        double X_min = 0.2;
        double X_max = 0.5;
        
        //Loop on both cuts in [X_min;X_max]
        Int_t dim = 0;
        //Increase N to increase (a,b) scan granularity!
        Int_t N = 5; double step = (X_max - X_min)/N;
        
        for(int i=0; i<3; i++){
          hname=to_string(i+1);
          
          //tau_cut1_vs_cut2_vs_sig[i] = new TH2D("tau_cut1_vs_cut2_vs_sig","tau_cut1_vs_cut2_vs_sig_"+hname,N,X_min,X_max,N,X_min,X_max,"a","b");
          tau_cut1_vs_cut2_vs_limit[i] = new TH2D("tau_cut1_vs_cut2_vs_limit","tau_cut1_vs_cut2_vs_limit_"+hname,N,X_min,X_max,N,X_min,X_max);
              
        }
        
        for(int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                        
                        for(int k=0; k<3; k++){
                                
                                a[k] = X_min + i * step;
                                b[k] = X_min + j * step;
                                if(a[k]<b[k]) continue;
                                
                                cout<<"i: "<<i<<" j: "<<j<<endl;
                                
                                
                                BDTOutput_x[k]->setRange("R_a",a[k],100); BDTOutput_x[k]->setRange("R_b",b[k],a[k]);
                                
                                
                                //computing areas in range [a;X_max]
                                N_s_1[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() * (BDTNormMC[k]->getValV());
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal());
                                
                                //computing areas in range [b;a]
                                N_s_2[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() * (BDTNormMC[k]->getValV());
                                N_b_2[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal());
                                
                                /*
                                N_s_1[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() * (BDTNormMC[k]->getValV()) * 4500.0/59.0;
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * 4500.0/59.0;
                                N_s_2[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() * (BDTNormMC[k]->getValV()) * 4500.0/59.0;
                                N_b_2[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * 4500.0/59.0;
                                */
                                
                                
                                if(N_b_1[k]>0.0&&N_b_2[k]>0.0){
                                        //S1[k] = N_s_1[k] / sqrt(N_s_1[k] + N_b_1[k]);
                                        //S2[k] = N_s_2[k] / sqrt(N_s_2[k] + N_b_2[k]);
                                        
                                        //S1, S2, S is now used to store limits, not significances
                                        S1[k] = sqrt( (2 * (N_s_1[k] + N_b_1[k]) * log(1 + (N_s_1[k]/N_b_1[k])) ) - 2 * N_s_1[k] );
                                        S2[k] = sqrt( (2 * (N_s_2[k] + N_b_2[k]) * log(1 + (N_s_2[k]/N_b_2[k])) ) - 2 * N_s_2[k] );
                                        
                                        //Combined significance
                                        //S[k] = sqrt(S1[k]*S1[k] + S2[k]*S2[k]);
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        command_a[k] = "python Combine/jupyternb/"+ card_modifier_name[k] +".py --luminosity " + std::to_string(59.0) + " --s " + std::to_string(N_s_1[k]) + " --b " + std::to_string(N_b_1[k]) + " --cuttype a";
                                        command_b[k] = "python Combine/jupyternb/"+ card_modifier_name[k] +".py --luminosity " + std::to_string(59.0) + " --s " + std::to_string(N_s_2[k]) + " --b " + std::to_string(N_b_2[k]) + " --cuttype b";
                                        system(command_a[k]);
                                        system(command_b[k]);
                                        
                                        command_a_and_b[k] = "combineCards.py "+combined_card_name[k]+"_a.txt "+combined_card_name[k]+"_b.txt > "+combined_card_name[k]+".txt";
                                        system(command_a_and_b[k]);
                                        command_run[k] = "combine -M AsymptoticLimits "+combined_card_name[k]+".txt --cl 0.9 -t -1  > out"+ to_string(k+1) +".txt";
                                        system(command_run[k]);
                                        
                                        
                                        
                                        std::ifstream f1("out"+ to_string(k+1) +".txt");
                                        std::string line;
                                        while (std::getline(f1, line)) {
                                          if (line.find("50.0%") != std::string::npos) {
                                            std::vector<std::string> linsp;
                                            std::istringstream iss(line);
                                            for (std::string s; iss >> s;) linsp.push_back(s);
                                            //cout << "linsp mu: " << linsp.back() << endl;
                                            S[k] = std::stod(linsp.back());
                                          }
                                        }
                                        
                                        //cout<<"N_s_1[k]: "<<N_s_1[k]<<" N_b_1[k]: "<<N_b_1[k]<<" N_s_2[k]: "<<N_s_2[k]<<" N_b_2[k]: "<<N_b_2[k]<<endl;
                                        //cout<<"S: "<<S[k]<<endl;
                                        
                                        
                                        if(!isnan(S[k])){
                                        
                                                //cout<<"N_s_1[k]: "<<N_s_1[k]<<" N_b_1[k]: "<<N_b_1[k]<<" N_s_2[k]: "<<N_s_2[k]<<" N_b_2[k]: "<<N_b_2[k]<<endl;
                                                //cout<<"S: "<<S[k]<<endl;
                                                
                                                //S1_list[k].push_back(S1[k]);
                                                //S2_list[k].push_back(S2[k]);
                                                a_list[k].push_back(a[k]);
                                                b_list[k].push_back(b[k]);
                                                S_list[k].push_back(S[k]);
                                                
                                                
                                                N_s_1_yield_list[k].push_back(N_s_1[k]);
                                                N_s_2_yield_list[k].push_back(N_s_2[k]);
                                                N_b_1_yield_list[k].push_back(N_b_1[k]);
                                                N_b_2_yield_list[k].push_back(N_b_2[k]);
                                                
                                                dim++;
                                                tau_cut1_vs_cut2_vs_limit[k]->Fill(a[k]+0.00000001,b[k]+0.00000001,S[k] );
                                        
                                        }
                                }
                        
                        }
                        
                }
        }
        
    //Taking absolute maximum of the combined significance
    
    double S_max[3];
    int S_maxIndex[3];
    float a_max[3];
    float b_max[3];
    double * a_array[3];
    double * sig_a_array[3];
    
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      S_max[i] = *min_element(S_list[i].begin(), S_list[i].end());
      S_maxIndex[i] = std::min_element(S_list[i].begin(),S_list[i].end()) - S_list[i].begin();
      cout<<"S_min[i]: "<<S_max[i]<<" S_minIndex[i]: "<<S_maxIndex[i]<<endl;
      
      a_max[i] = a_list[i].at(S_maxIndex[i]);
      b_max[i] = b_list[i].at(S_maxIndex[i]);
      
      cout<<"b cut: "<<b_max[i]<<" a cut: "<<a_max[i]<<endl;
      cout<<"a signal yield: "<< N_s_1_yield_list[i].at(S_maxIndex[i]) <<" a bkg yield: "<< N_b_1_yield_list[i].at(S_maxIndex[i]) <<endl;
      cout<<"b,a signal yield: "<< N_s_2_yield_list[i].at(S_maxIndex[i]) <<" b,a bkg yield: "<< N_b_2_yield_list[i].at(S_maxIndex[i]) <<endl;
      
      //a_array[i] = &a_list[i][0];
      //sig_a_array[i] = &S1_list[i][0];
      
    }
    
    
    /*
    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 1800, 600);
    canvas4->Divide(3, 1);
    
    TGraph * Sig_Plot[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas4->cd( i+1 );
      Sig_Plot[i] = new TGraph(a_list[i].size(),a_array[i],sig_a_array[i]);
      Sig_Plot[i]->SetTitle("Significance vs. BDT Cut Value, #tau_"+cat_label[i]+";BDT Cut Value;Significance");
      Sig_Plot[i]->Draw();
      
    }
    */
    
    
    
    /*
    //std::vector<double> lumi = {97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0};
    std::vector<double> lumi = {97.7, 129.0};
    std::vector<double> ZTT_taumu_lim;
    std::vector<double> ZTT_tauh_lim;
    std::vector<double> ZTT_taue_lim;
    TString command1 = "";
    TString command2 = "";
    TString command3 = "";
    for(int i=0; i<lumi.size(); i++){
      command1 = "python Combine/jupyternb/ZTT_taumu_test.py --luminosity " + std::to_string(lumi[i]);
      command2 = "python Combine/jupyternb/ZTT_tauh_test.py --luminosity " + std::to_string(lumi[i]);
      command3 = "python Combine/jupyternb/ZTT_taue_test.py --luminosity " + std::to_string(lumi[i]);
      system(command1);
      system(command2);
      system(command3);
      system("combine -M AsymptoticLimits ZTT_taumu_Combined_Mod.txt --cl 0.9 -t -1  > out1.txt");
      system("combine -M AsymptoticLimits ZTT_tauh_Combined_Mod.txt --cl 0.9 -t -1 > out2.txt");
      system("combine -M AsymptoticLimits ZTT_taue_Combined_Mod.txt --cl 0.9 -t -1 > out3.txt");
      
      
      std::ifstream f1("out1.txt");
      std::string line_mu;
      while (std::getline(f1, line_mu)) {
        if (line_mu.find("50.0%") != std::string::npos) {
          std::vector<std::string> linsp;
          std::istringstream iss(line_mu);
          for (std::string s; iss >> s;) linsp.push_back(s);
          cout << "linsp mu: " << linsp.back() << endl;
          ZTT_taumu_lim.push_back(std::stod(linsp.back()));
        }
      }
      
      std::ifstream f2("out2.txt");
      std::string line_h;
      while (std::getline(f2, line_h)) {
        if (line_h.find("50.0%") != std::string::npos) {
          std::vector<std::string> linsp;
          std::istringstream iss(line_h);
          for (std::string s; iss >> s;) linsp.push_back(s);
          cout << "linsp h: " << linsp.back() << endl;
          ZTT_tauh_lim.push_back(std::stod(linsp.back()));
        }
      }
      
      std::ifstream f3("out3.txt");
      std::string line_e;
      while (std::getline(f3, line_e)) {
        if (line_e.find("50.0%") != std::string::npos) {
          std::vector<std::string> linsp;
          std::istringstream iss(line_e);
          for (std::string s; iss >> s;) linsp.push_back(s);
          cout << "linsp e: " << linsp.back() << endl;
          ZTT_taue_lim.push_back(std::stod(linsp.back()));
        }
      }
      
    }
    */
    
    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 1800, 600);
    canvas5->Divide(3, 1);
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas5->cd( i+1 );
      tau_cut1_vs_cut2_vs_limit[i]->SetTitle("Limit vs. BDT Cut Values, #tau_"+cat_label[i]+";a;b");
      tau_cut1_vs_cut2_vs_limit[i]->SetStats(0);
      tau_cut1_vs_cut2_vs_limit[i]->SetMinimum(S_max[i]);
      tau_cut1_vs_cut2_vs_limit[i]->Draw("colz");
      
    }
    
    
}