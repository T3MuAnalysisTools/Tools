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
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace RooFit;

void makeYield () 
{
    //histo MC
    TFile * file_tauh = TFile::Open("LOCAL_COMBINED_ztau3mutauh_default_LumiScaled.root");
    TH1D * tauh_Vis_Mass   = (TH1D*)file_tauh->Get("ztau3mutauh_default_Cut_11_Nminus0_VisMass_MC4");
    TH1D * tauh_Vis_Mass_Dat   = (TH1D*)file_tauh->Get("ztau3mutauh_default_Cut_11_Nminus0_VisMass_Data");
    TH1D * tauh_T3Mu   = (TH1D*)file_tauh->Get("ztau3mutauh_default_TripletMassMC4");
    TH1D * tauh_T3Mu_Dat   = (TH1D*)file_tauh->Get("ztau3mutauh_default_TripletMassData");
    
    TFile * file_taumu = TFile::Open("LOCAL_COMBINED_ztau3mutaumu_default_LumiScaled.root");
    TH1D * taumu_Vis_Mass   = (TH1D*)file_taumu->Get("ztau3mutaumu_default_Cut_9_Nminus0_VisMass_MC3");
    TH1D * taumu_Vis_Mass_Dat   = (TH1D*)file_taumu->Get("ztau3mutaumu_default_Cut_9_Nminus0_VisMass_Data");
    TH1D * taumu_T3Mu   = (TH1D*)file_taumu->Get("ztau3mutaumu_default_TripletMassMC3");
    TH1D * taumu_T3Mu_Dat   = (TH1D*)file_taumu->Get("ztau3mutaumu_default_TripletMassData");
    
    TFile * file_taue = TFile::Open("LOCAL_COMBINED_ztau3mutaue_default_LumiScaled.root");
    TH1D * taue_Vis_Mass   = (TH1D*)file_taue->Get("ztau3mutaue_default_Cut_9_Nminus0_VisMass_MC2");
    TH1D * taue_Vis_Mass_Dat   = (TH1D*)file_taue->Get("ztau3mutaue_default_Cut_9_Nminus0_VisMass_Data");
    TH1D * taue_T3Mu   = (TH1D*)file_taue->Get("ztau3mutaue_default_TripletMassMC2");
    TH1D * taue_T3Mu_Dat   = (TH1D*)file_taue->Get("ztau3mutaue_default_TripletMassData");

    //cout << "tauh MC count: " << tauh_Vis_Mass->Integral() << " tauh Data count: " << tauh_Vis_Mass_Dat->Integral() << endl;
    //cout << "taumu MC count: " << taumu_Vis_Mass->Integral() << " taumu Data count: " << taumu_Vis_Mass_Dat->Integral() << endl;
    //cout << "taue MC count: " << taue_Vis_Mass->Integral() << " taue Data count: " << taue_Vis_Mass_Dat->Integral() << endl;
    
    double signal_low = 1.7233333;
    double signal_high = 1.8333333;
    //For tau h
    RooRealVar InvMass1("InvMass1","3#mu inv. mass (GeV), #tau_{h}",1.3,2.2);
    InvMass1.setRange("R1",1.3,signal_low); //background    
    InvMass1.setRange("R2",signal_high,2.2); //background
    InvMass1.setRange("R3",1.73,1.82); //signal
    InvMass1.setRange("R4",signal_low,signal_high); //signal range used to obtain signal and background yields.
    RooPolynomial poly1("poly1","poly(InvMass1)",InvMass1);
    RooDataHist data1("data1", "data1", InvMass1, Import(*tauh_T3Mu_Dat));
    RooDataHist mc1("mc1", "mc1", InvMass1, Import(*tauh_T3Mu));
    RooRealVar LineNorm1("LineNorm1", "LineNorm1", 0.5,0.01,150);
    RooAddPdf pdf1("pdf1", "pdf1", RooArgList(poly1), RooArgList(LineNorm1));
    RooFitResult *fitresult1;
    fitresult1 = pdf1.fitTo(data1, Range("R1,R2"), Save());
    
    RooRealVar mean1("mean1", "mean1", 1.776,0.,5.);
    RooRealVar sigma1("sigma1", "sigma1", 0.5,0.001,10);
    RooRealVar GaussNorm1("GaussNorm1", "GaussNorm1", 0.5,0.001,1.0);
    RooGaussian Gauss1("Gauss1", "Gauss1", InvMass1, mean1, sigma1);
    RooAddPdf mc_pdf1("mc_pdf1", "mc_pdf1", RooArgList(Gauss1), RooArgList(GaussNorm1));
    RooFitResult *mc_fitresult1;
    mc_fitresult1 = mc_pdf1.fitTo(mc1, Range("R3"), Save());
    
    
    //For tau mu
    RooRealVar InvMass2("InvMass2","3#mu inv. mass (GeV), #tau_{#mu}",1.3,2.2);
    InvMass2.setRange("R1",1.3,signal_low); //background    
    InvMass2.setRange("R2",signal_high,2.2); //background
    InvMass2.setRange("R3",1.73,1.82); //signal
    InvMass2.setRange("R4",signal_low,signal_high); //signal range for yield
    RooPolynomial poly2("poly2","poly(InvMass2)",InvMass2);
    RooDataHist data2("data2", "data2", InvMass2, Import(*taumu_T3Mu_Dat));
    RooDataHist mc2("mc2", "mc2", InvMass2, Import(*taumu_T3Mu));
    RooRealVar LineNorm2("LineNorm2", "LineNorm2", 50.0,0.01,150);
    RooAddPdf pdf2("pdf2", "pdf2", RooArgList(poly2), RooArgList(LineNorm2));
    RooFitResult *fitresult2;
    fitresult2 = pdf2.fitTo(data2, Range("R1,R2"), Save());
    
    RooRealVar mean2("mean2", "mean2", 1.776,0.,5.);
    RooRealVar sigma2("sigma2", "sigma2", 0.5,0.001,10);
    RooRealVar GaussNorm2("GaussNorm2", "GaussNorm2", 0.5,0.001,1.0);
    RooGaussian Gauss2("Gauss2", "Gauss2", InvMass2, mean2, sigma2);
    RooAddPdf mc_pdf2("mc_pdf2", "mc_pdf2", RooArgList(Gauss2), RooArgList(GaussNorm2));
    RooFitResult *mc_fitresult2;
    mc_fitresult2 = mc_pdf2.fitTo(mc2, Range("R3"), Save());
    
    
    //For tau e
    RooRealVar InvMass3("InvMass3","3#mu inv. mass (GeV), #tau_{e}",1.3,2.2);
    InvMass3.setRange("R1",1.3,signal_low); //background    
    InvMass3.setRange("R2",signal_high,2.2); //background
    InvMass3.setRange("R3",1.71,1.82); //signal
    InvMass3.setRange("R4",signal_low,signal_high); //signal range for yield
    RooPolynomial poly3("poly3","poly(InvMass3)",InvMass3);
    RooDataHist data3("data3", "data3", InvMass3, Import(*taue_T3Mu_Dat));
    RooDataHist mc3("mc3", "mc3", InvMass3, Import(*taue_T3Mu));
    RooRealVar LineNorm3("LineNorm3", "LineNorm3", 50.0,0.01,150);
    RooAddPdf pdf3("pdf3", "pdf3", RooArgList(poly3), RooArgList(LineNorm3));
    RooFitResult *fitresult3;
    fitresult3 = pdf3.fitTo(data3, Range("R1,R2"), Save());
    
    RooRealVar mean3("mean3", "mean3", 1.776,0.,5.);// Not working?
    RooRealVar sigma3("sigma3", "sigma3", 0.5,0.001,10);
    RooRealVar GaussNorm3("GaussNorm3", "GaussNorm3", 0.5,0.001,1.0);
    RooGaussian Gauss3("Gauss3", "Gauss3", InvMass3, mean3, sigma3);
    RooAddPdf mc_pdf3("mc_pdf3", "mc_pdf3", RooArgList(Gauss3), RooArgList(GaussNorm3));
    RooFitResult *mc_fitresult3;
    mc_fitresult3 = mc_pdf3.fitTo(mc3, Range("R3"), Save());
    
    
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 1800, 600);
    canvas1->Divide(3, 1);
    
    
    
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
    
    cout << "mc_pdf1_integral: " << GaussNorm1.getValV()*mc_pdf1_integral_restricted << " pdf1_integral: " << LineNorm1.getValV()*pdf1_integral_restricted << endl;
    cout << "mc_pdf2_integral: " << GaussNorm2.getValV()*mc_pdf2_integral_restricted << " pdf2_integral: " << LineNorm2.getValV()*pdf2_integral_restricted << endl;
    cout << "mc_pdf3_integral: " << GaussNorm3.getValV()*mc_pdf3_integral_restricted << " pdf3_integral: " << LineNorm3.getValV()*pdf3_integral_restricted << endl;
    
    double scaling1 = mc1.sumEntries("1")/tauh_T3Mu->GetEntries();
    double scaling2 = mc2.sumEntries("1")/taumu_T3Mu->GetEntries();
    double scaling3 = mc3.sumEntries("1")/taue_T3Mu->GetEntries();
    
    cout << "Unscaled mc1: " << nSignal1_restricted/scaling1 << " scaling 1: " << scaling1 << endl;// Getting the unweighted content of the signal histogram
    cout << "Unscaled mc2: " << nSignal2_restricted/scaling2 << " scaling 2: " << scaling2 << endl;
    cout << "Unscaled mc3: " << nSignal3_restricted/scaling3 << " scaling 3: " << scaling3 << endl;
    
    canvas1->cd( 1 );
    RooPlot *xFrame1 = InvMass1.frame();
    data1.plotOn(xFrame1);
    //mc1.plotOn(xFrame1);
    pdf1.plotOn(xFrame1,LineColor(4),LineWidth(2), Normalization(nData1, RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
    mc_pdf1.plotOn(xFrame1,LineColor(1),LineWidth(2), Normalization(nSignal1, RooAbsReal::NumEvent),ProjectionRange("R3"));
    xFrame1->SetTitle("3#mu inv. mass (GeV), #tau_{h}");
    xFrame1->SetXTitle("3#mu inv. mass (GeV)");
    xFrame1->SetYTitle("Events");
    xFrame1->Draw();
    
    canvas1->cd( 2 );
    RooPlot *xFrame2 = InvMass2.frame();
    data2.plotOn(xFrame2);
    //mc2.plotOn(xFrame2);
    pdf2.plotOn(xFrame2,LineColor(4),LineWidth(2), Normalization(nData2, RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
    mc_pdf2.plotOn(xFrame2,LineColor(1),LineWidth(2), Normalization(nSignal2, RooAbsReal::NumEvent),ProjectionRange("R3"));
    xFrame2->SetTitle("3#mu inv. mass (GeV), #tau_{#mu}");
    xFrame2->SetXTitle("3#mu inv. mass (GeV)");
    xFrame2->SetYTitle("Events");
    xFrame2->Draw();
    
    canvas1->cd( 3 );
    RooPlot *xFrame3 = InvMass3.frame();
    data3.plotOn(xFrame3);
    //mc3.plotOn(xFrame3);
pdf3.plotOn(xFrame3);
    pdf3.plotOn(xFrame3,LineColor(4),LineWidth(2), Normalization(nData3, RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
    mc_pdf3.plotOn(xFrame3,LineColor(1),LineWidth(2), Normalization(nSignal3, RooAbsReal::NumEvent),ProjectionRange("R3"));
    xFrame3->SetTitle("3#mu inv. mass (GeV), #tau_{e}");
    xFrame3->SetXTitle("3#mu inv. mass (GeV)");
    xFrame3->SetYTitle("Events");
    xFrame3->Draw();
    
    
}