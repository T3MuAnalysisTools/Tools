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

void TemplateFit()
{
  //  make fake data--------------
  TH1F * fraction1 = new TH1F("fraction1","fraction1",40,-20.,20);
  TH1F * fraction2 = new TH1F("fraction2","fraction2",40,-20.,20);

  for(int ij=0;ij<10000;ij++) {// Filling the histograms with no of Gaussian Distribution 
    fraction1->Fill(gRandom->Gaus(2.,3.));  // mean  =2, sigma = 3.
    fraction2->Fill(gRandom->Gaus(-1.,.8)); // mean = -1 sigma = 0.8
  }


  TH1F *h_data = new TH1F("h_data","  ",40,-20,20); // Create 'data' which contains template1 and 0.5 of template2
  h_data->Add(fraction1);
  h_data->Add(fraction2,0.5);  // so 'data' consists of 15000 events with 2:1 fraction of MC1 and MC2  so the fraction1   is 2/3 = 0.666666
  //-------------------------------




  //   make templates -----------
  TH1F * template1 = new TH1F("template1","template1",40,-20.,20);
  TH1F * template2 = new TH1F("template2","template2",40,-20.,20);

  for (int i=0;i<50000;i++) // 
    {
      Double_t life1=fraction1->GetRandom();
      template1->Fill(life1); //  fill template 1 according to fraction1 
      Double_t life2=fraction2->GetRandom();
      template2->Fill(life2);//  fill template 2 according to fraction1 
    }
  //-------------------------------




  
  RooRealVar x("x","x",-20,20.0);

  RooDataHist data("data","data",x,h_data);
  RooDataHist mc1("mc1","MC 1",x,template1);
  RooDataHist mc2("mc2","MC 2",x,template2);

  RooRealVar fraction("fraction","fraction",0.1,  0,  1); // These are variables for output

  // Make PDF from MC histograms
  RooHistPdf modelmc1("modelmc1","modelmc1",x, mc1);
  RooHistPdf modelmc2("modelmc2","modelmc2",x, mc2);


  RooAddPdf model("model","model",RooArgList(modelmc1,modelmc2),fraction);  // define the model as a sum of two templates with their fraction which is a free parameter


  // Plot the pre-fit histograms
  RooPlot* dframe = x.frame(Title("Data"));
  data.plotOn(dframe);

  RooPlot* mc1frame = x.frame(Title("Template  1"));
  mc1.plotOn(mc1frame);

  RooPlot* mc2frame = x.frame(Title("Template  2"));
  mc2.plotOn(mc2frame);


  TCanvas* c = new TCanvas("roofit_example","RooFit FractionFit Example",800,600);
  c->Divide(1,3);
  gROOT->SetStyle("Plain"); // Removes gray background from plot
  c->cd(1) ; gPad->SetLeftMargin(0.15) ;   dframe->GetYaxis()->SetTitleOffset(1.4) ;   dframe->Draw();
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; mc1frame->GetYaxis()->SetTitleOffset(1.4) ; mc1frame->Draw();
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; mc2frame->GetYaxis()->SetTitleOffset(1.4) ; mc2frame->Draw();
  

  model.fitTo(data); //  Do the fit of  model 



  // plot fit resutls
  RooPlot* fitFrame=x.frame(Bins(50), Title("Fit Model"));
  model.paramOn(fitFrame);
  data.plotOn(fitFrame, RooFit::LineColor(kRed));
  model.plotOn(fitFrame, LineStyle(kDashed));
  model.plotOn(fitFrame, Components("modelmc1"), LineColor(kGreen));
  model.plotOn(fitFrame, Components("modelmc2"), LineColor(kBlue));

  TCanvas* c6 = new TCanvas("c6","Fit Model",800,600);
  gROOT->SetStyle("Plain"); // Removes gray background from plots
  gPad->SetLeftMargin(0.15) ;   
  fitFrame->GetYaxis()->SetTitleOffset(1.4) ;   
  fitFrame->Draw();  

 
  // ---  add legend 
  TLegend *legmc = new TLegend(0.12,0.70,0.43,0.86);
  legmc->AddEntry(fitFrame->getObject(1),"Data","LPE");
  legmc->AddEntry(fitFrame->getObject(2),"Fit","LPE");
  legmc->AddEntry(fitFrame->getObject(3),"Template 1","LPE");
  legmc->AddEntry(fitFrame->getObject(4),"Template 2","LPE");
  legmc->Draw();  
  
} 
