#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1F.h"
#include <cmath>
#include "TRandom.h"


using namespace RooFit;

void DecayFit()
{

  TFile *f = new TFile("/afs/cern.ch/work/n/nimenend/Analysis/workdir_Decay_Length_2018_Feb_27_2020/LOCAL_COMBINED_dstophipi_default.root", "READ");
  TH1F *template1, *template2, *h_data;
  template1 = (TH1F*)f->Get("dstophipi_default_DecayLength_promptMC1");
  template2 = (TH1F*)f->Get("dstophipi_default_DecayLength_non_promptMC1");
  h_data = (TH1F*)f->Get("dstophipi_default_DecayLengthData");

  
  RooRealVar DecayLength("DecayLength","Proper Decay Length",-20,20.0);

  RooDataHist data("data","data",DecayLength,h_data);
  RooDataHist mc1("mc1","MC 1",DecayLength,template1);
  RooDataHist mc2("mc2","MC 2",DecayLength,template2);

  RooRealVar fraction("fraction","fraction",0.1,  0,  1); // These are variables for output

  // Make PDF from MC histograms
  RooHistPdf modelmc1("modelmc1","modelmc1",DecayLength, mc1);
  RooHistPdf modelmc2("modelmc2","modelmc2",DecayLength, mc2);


  RooAddPdf model("model","model",RooArgList(modelmc1,modelmc2),fraction);  // define the model as a sum of two templates with their fraction which is a free parameter


  // Plot the pre-fit histograms
  RooPlot* dframe = DecayLength.frame(Title("Data"));
  data.plotOn(dframe);

  RooPlot* mc1frame = DecayLength.frame(Title("Prompt"));
  mc1.plotOn(mc1frame);

  RooPlot* mc2frame = DecayLength.frame(Title("Non_Prompt"));
  mc2.plotOn(mc2frame);


  TCanvas* c = new TCanvas("c","Templates Example",800,600);
  c->Divide(1,3);
  gROOT->SetStyle("Plain"); // Removes gray background from plot
  c->cd(1) ; gPad->SetLeftMargin(0.15) ;   dframe->GetYaxis()->SetTitleOffset(1.4) ;   dframe->Draw();
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; mc1frame->GetYaxis()->SetTitleOffset(1.4) ; mc1frame->Draw();
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; mc2frame->GetYaxis()->SetTitleOffset(1.4) ; mc2frame->Draw();
  

  model.fitTo(data); //  Do the fit of  model 



  // plot fit resutls
  RooPlot* fitFrame=DecayLength.frame(Bins(50));
  model.paramOn(fitFrame);
  data.plotOn(fitFrame, RooFit::LineColor(kRed));
  model.plotOn(fitFrame, LineStyle(kDashed));
  model.plotOn(fitFrame, Components("modelmc1"), LineColor(kGreen));
  model.plotOn(fitFrame, Components("modelmc2"), LineColor(kBlue));
  fitFrame->chiSquare() ;



  TCanvas* c6 = new TCanvas("c6","Fit Model",800,600);
  gROOT->SetStyle("Plain"); // Removes gray background from plots
  gPad->SetLeftMargin(0.15) ;   

  c6->SetFrameLineWidth(3);
  c6->SetTickx();
  c6->SetTicky();

  fitFrame->GetYaxis()->SetTitleOffset(1.4) ;   
  fitFrame->Draw();  

 
  // ---  add legend 
  TLegend *legmc = new TLegend(0.16,0.70,0.43,0.86);
  legmc->AddEntry(fitFrame->getObject(1),"Data","LPE");
  legmc->AddEntry(fitFrame->getObject(2),"Fit","LPE");
  legmc->AddEntry(fitFrame->getObject(3),"Prompt","LPE");
  legmc->AddEntry(fitFrame->getObject(4),"Non-Prompt","LPE");
  legmc->Draw();  
  
} 
