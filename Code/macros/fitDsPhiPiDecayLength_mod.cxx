// include add header files
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
#include "tdrstyle.C"
#include "CMSStyle_cmslumi.C"

using namespace RooFit;

void fitDsPhiPiDecayLength_mod(TString inputfile, TString inputfile_mc){

		  double DsScale = 0.0423*1.12;
		  double lumiScale = (1000*1.449*354.50*59.6/900506);
		  double eraAScale = 4042.23/7368.01;
		  double eraBScale = 1610.00/3090.22;
		  double eraCScale =  2014.95/3851.56;
		  double eraDScale = 5560.00/10185.06;
		  

		  TFile* rootfile = new TFile(inputfile, "READ");
		  TFile* rootfile_mc = new TFile(inputfile_mc, "READ");
		  TTree* tree = (TTree*)rootfile->Get("tree");
		  TTree* tree_mc = (TTree*)rootfile_mc->Get("tree");
		  // get data
		  TH1F *h_decay_length_data = new TH1F("h_decay_length_data", "Data", 50, 0.0, 0.1);
		  TH1F *h_decay_length_prompt_mc = new TH1F("h_decay_length_prompt_mc", "Prompt", 50, 0.0, 0.1);
		  TH1F *h_decay_length_non_prompt_mc = new TH1F("h_decay_length_non_prompt_mc", "Non-Prompt", 50, 0.0, 0.1);
		  TH1F *h_decay_length_fake_mc = new TH1F("h_decay_length_fake_mc", "Fakes", 50, 0.0, 0.1);

		  int nentries = tree->GetEntriesFast();
		  int nentries_mc = tree_mc->GetEntriesFast();

		  int era = -1;
		  Float_t DecayLength = -1;
		  Float_t TripleMass = 0.;
		  Float_t eventWeight = 0.;
		  Float_t PVSV_dxy = 0.0;
		  Float_t PVSV_dz = 0.0;
		  int var_isPrompt = -1;
		  int ds_motherPdgId = -1;
		  int dstar_mom = -1;

		  float ds_eta;
		  float ds_pt;
		  float ds_phi;

		  tree->SetBranchAddress("era", &era);
		  tree->SetBranchAddress("TripleMass", &TripleMass);
		  tree->SetBranchAddress("DecayLength", &DecayLength);
		  tree->SetBranchAddress("eventWeight", &eventWeight);
		  tree->SetBranchAddress("var_isPrompt", &var_isPrompt);
		  tree->SetBranchAddress("PVSV_dxy", &PVSV_dxy);
		  tree->SetBranchAddress("PVSV_dz", &PVSV_dz);
		  tree->SetBranchAddress("ds_eta", &ds_eta);
		  
		  tree_mc->SetBranchAddress("era", &era);
		  tree_mc->SetBranchAddress("TripleMass", &TripleMass);
		  tree_mc->SetBranchAddress("DecayLength", &DecayLength);
		  tree_mc->SetBranchAddress("eventWeight", &eventWeight);
		  tree_mc->SetBranchAddress("var_isPrompt", &var_isPrompt);
		  tree_mc->SetBranchAddress("dstar_mom", &dstar_mom);
		  tree_mc->SetBranchAddress("ds_motherPdgId", &ds_motherPdgId);
		  tree_mc->SetBranchAddress("PVSV_dxy", &PVSV_dxy);
		  tree_mc->SetBranchAddress("PVSV_dz", &PVSV_dz);
		  tree_mc->SetBranchAddress("ds_eta", &ds_eta);

		  for (int i=0; i<nentries; ++i){
					 if (i%1000==0) cout<<"Processing "<<i<<"/"<<nentries<<endl;
					 tree->GetEntry(i);
					 double scale = 0.0;
					 if (DecayLength>0.099) DecayLength = 0.099;
					 double distance3d = sqrt(pow(PVSV_dxy,2)+pow(PVSV_dz,2));
					 if (distance3d<0.01) continue;
					 if (ds_eta>0.8) continue;
					 if (era>-1){
								if (era==0) scale = eraAScale;
								else if (era==1) scale = eraBScale;
								else if (era==2) scale = eraCScale;
								else if (era==3) scale = eraDScale;
								if (TripleMass<1.80 && TripleMass>=1.70) h_decay_length_data->Fill(DecayLength, -scale);
								else if (TripleMass<2.01 && TripleMass>=1.93) { h_decay_length_data->Fill(DecayLength, 1.0); }
					 }

		  }

		  for (int i=0; i<nentries_mc; ++i){
					 if (i%1000==0) cout<<"Processing "<<i<<"/"<<nentries<<endl;
					 tree_mc->GetEntry(i);
					 double scale = 0.0;
					 if (DecayLength>0.099) DecayLength = 0.099;
					 double distance3d = sqrt(pow(PVSV_dxy,2)+pow(PVSV_dz,2));
					 if (distance3d<0.01) continue;
					 if (ds_eta>0.8) continue;
					 //if (DecayLength<0.01) continue;
					 if (TripleMass<2.01 && TripleMass>=1.93){
								if ( ( ds_motherPdgId<9 && ds_motherPdgId>0 ) ||
													 (( ds_motherPdgId==433 || ds_motherPdgId==10431 ) && ( (abs(dstar_mom)>400 && abs(dstar_mom)<500) || abs(dstar_mom)==2212 )) ||
													 ( ds_motherPdgId>2100 && ds_motherPdgId<2205 ) ||
													 (( ds_motherPdgId==433 || ds_motherPdgId==10431 ) && ( abs(dstar_mom)>2100 && abs(dstar_mom)<2205 ))
									) {
										  h_decay_length_prompt_mc->Fill(DecayLength, DsScale*lumiScale*eventWeight);
								}
								else if (
													 ( ds_motherPdgId>500 && ds_motherPdgId<600 ) ||
													 (( ds_motherPdgId==433 || ds_motherPdgId==10431 ) && ( abs(dstar_mom)>500 && abs(dstar_mom)<600 )) ||
													 (ds_motherPdgId>5100) ||
													 (( ds_motherPdgId==433 || ds_motherPdgId==10431 ) && ( abs(dstar_mom)>5100))
										  ) {
										  h_decay_length_non_prompt_mc->Fill(DecayLength, DsScale*lumiScale*eventWeight);
								}
								else if (ds_motherPdgId==0) h_decay_length_fake_mc->Fill(DecayLength, DsScale*lumiScale*eventWeight);
								else {
										  cout<<ds_motherPdgId<<endl;
										  cout<<dstar_mom<<endl;
								}
					 }
		  }

		  cout<<"Fitting histograms"<<endl;
		  RooRealVar x("x","x",0,0.1);

		  RooDataHist data("data","data",x,h_decay_length_data);
		  RooDataHist prompt("prompt","prompt",x,h_decay_length_prompt_mc);
		  RooDataHist non_prompt("non_prompt","non_prompt",x,h_decay_length_non_prompt_mc);
		  RooDataHist fake("fake","fake",x,h_decay_length_fake_mc);

		  cout<<"Models"<<endl;
		  RooRealVar fraction("fraction","fraction", 0.21, 0, 1);
		  RooRealVar fakeRate("fakeRate","fakeRate", 0.02, 0, 1);
		  // Make PDF from MC histograms
		  RooHistPdf model_prompt("model_prompt","model_prompt",x, prompt);
		  RooHistPdf model_non_prompt("model_non_prompt","model_non_prompt",x, non_prompt);
		  RooHistPdf model_fake("fake","fake",x, fake);


		  RooAddPdf model("model","model",RooArgList(model_non_prompt,model_fake, model_prompt), RooArgList(fraction,fakeRate));  // define the model as a sum of two templates with their fraction which is a free parameter
		  cout<<"Done pdf; frame"<<endl;
		  // Plot the pre-fit histograms
		  RooPlot* dframe = x.frame(Title("Data"));
		  data.plotOn(dframe);

		  RooPlot* promptframe = x.frame(Title("Prompt D_{s}"));
		  prompt.plotOn(promptframe);

		  RooPlot* non_promptframe = x.frame(Title("Non-Prompt D_{s}"));
		  non_prompt.plotOn(non_promptframe);

		  RooPlot* fakeframe = x.frame(Title("Fakes"));
		  fake.plotOn(fakeframe);

		  TCanvas* c = new TCanvas("c","",660,660);
		  c->Divide(1,3);
		  gROOT->SetStyle("Plain"); // Removes gray background from plot
		  c->cd(1) ; gPad->SetLeftMargin(0.15) ;   dframe->GetYaxis()->SetTitleOffset(1.4) ;   dframe->Draw();
		  c->cd(2) ; gPad->SetLeftMargin(0.15) ; promptframe->GetYaxis()->SetTitleOffset(1.4) ; promptframe->Draw();
		  c->cd(3) ; gPad->SetLeftMargin(0.15) ; non_promptframe->GetYaxis()->SetTitleOffset(1.4) ; non_promptframe->Draw();


		  auto fitResult = model.fitTo(data, SumW2Error(true), Save(true)); //  Do the fit of  model 

		  setTDRStyle();
		  CMS_lumi(c, 4, 0);
		  //CMSStyle();
		  // plot fit resutls
		  RooPlot* fitFrame=x.frame(Bins(50));

		  //model.paramOn(fitFrame);
		  data.plotOn(fitFrame, RooFit::LineColor(kRed));
		  model.plotOn(fitFrame, LineStyle(kDashed));
		  model.plotOn(fitFrame, Components("model_prompt"), LineColor(kGreen));
		  model.plotOn(fitFrame, Components("model_non_prompt"), LineColor(kBlue));
		  model.plotOn(fitFrame, Components("fake"), LineColor(kOrange+10));

		  TCanvas* c6 = new TCanvas("c6","Fit Model",660,660);
		  gROOT->SetStyle("Plain"); // Removes gray background from plots

		  fitFrame->GetYaxis()->SetTitleOffset(1.7);   
		  fitFrame->GetXaxis()->SetTitle("Decay Length (cm)");
		  fitFrame->Draw();  

		  setTDRStyle();
		  CMS_lumi(c6, 4, 0);
		  // ---  add legend 
		  TLegend *legmc = new TLegend(0.61,0.65,0.88,0.81);
		  legmc->AddEntry(fitFrame->getObject(0),"Data","LPE");
		  legmc->AddEntry(fitFrame->getObject(1),"Fit","LE");
		  legmc->AddEntry(fitFrame->getObject(2),"Prompt D_{s}","LE");
		  legmc->AddEntry(fitFrame->getObject(3),"D_{s} from B decays","LE");
		  legmc->AddEntry(fitFrame->getObject(4),"Fakes","LE");
		  legmc->Draw("p");
		  gStyle->SetStatStyle(0);
		  gStyle->SetTitleStyle(0);
		  gROOT->ForceStyle();
		  TGaxis().SetMaxDigits(3);
		  TGaxis().SetExponentOffset(-0.075, -0.01, "y");

		  TLatex latex;
		  latex.SetTextFont(42);
		  latex.SetTextAlign(13);
		  latex.SetTextSize(0.03);
		  string fraction_str = "Fraction = "+to_string(fraction.getVal()).substr(0,5) + 
					 " +/- " +to_string(fraction.getPropagatedError(*fitResult)).substr(0,5);
		  latex.DrawLatex(0.05,1000, fraction_str.c_str());

		  c->SaveAs("decay_length_components_2018_corrected.pdf"); 
		  c6->SaveAs("decay_length_2018_corrected.pdf");
		  rootfile->Close();
}
