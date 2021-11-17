#ifndef BDT_OPTIMAL_CUT_H_
#define BDT_OPTIMAL_CUT_H_

#include "TH1F.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TCanvas.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <array>

using namespace std;

struct BDTcut3d { 
  float a; 
  float b; 
  float c; 
};

double TH1_integral(TH1F *h, float xmin, float xmax){
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);

    return integral;
}

double log_significance(double S, double B){
    double significance = 0;
    significance = sqrt(2*( (S+B) * log( 1+S/B ) - S));
    //cout<<"log sign is "<<significance<<" while S/sqrt(S + B) gives "<< S/sqrt(S + B) <<endl;
    return significance;
}


BDTcut3d Get_BDT_cut_3D(TString categ) 
{
    TString file_name = workdir+TMVA_inputpath+"outputTree.root";
    TH1F *h_test_signal;
    TH1F *h_test_bkg;
    TH1F *h_test_signal2;
    TH1F *h_test_bkg2;

    //open input files
    TChain *t = new TChain("outputTree");
    t->Add(file_name); 
    std::cout<<"Opened input file: "<<file_name<<std::endl;

    TString cat = ""; 
    if (categ=="A") cat = "0";
    if (categ=="B") cat = "1";
    if (categ=="C") cat = "2";
    TString signal = "weight*(isMC>0  && category =="+cat+")";
    TString bkg    = "weight*(isMC==0 && category =="+cat+" && isSB==1)";

    Int_t N = 400;
    //bdt score distribution
    TString binning = "(400,-1,1)"; 
    t->Draw("bdt>>h_test_bkg"+binning, bkg);
    h_test_bkg = (TH1F *)gDirectory->Get("h_test_bkg");
    h_test_bkg2 = (TH1F *)gDirectory->Get("h_test_bkg");
    t->Draw("bdt>>h_test_signal"+binning, signal);
    h_test_signal = (TH1F *)gDirectory->Get("h_test_signal");
    h_test_signal2 = (TH1F *)gDirectory->Get("h_test_signal");

    h_test_signal->SetDirectory(0);
    h_test_bkg->SetDirectory(0);
    h_test_signal2->SetDirectory(0);
    h_test_bkg2->SetDirectory(0);

    
    TH3F *h3 = new TH3F("h3","test",N,-1,1,N,-1,1,N,-1,1);
    float a,b,c;
    double N_s_1, N_b_1;
    double N_s_2, N_b_2;
    double N_s_3, N_b_3;
    double S1, S2, S3, S;

    //Bkg is scaled depending on categ to fit signal mass window
    double bkg_scale = 1;
    if(categ.Contains("A")) bkg_scale = 4. * 12./(380.);
    if(categ.Contains("B")) bkg_scale = 4. * 19./(380.);
    if(categ.Contains("C")) bkg_scale = 4. * 25./(380.);
    cout<<"bkg_scale = "<<bkg_scale<<endl;

    double X_min = std::min(h_test_signal->GetXaxis()->GetXmin(), h_test_signal->GetXaxis()->GetXmin());
    double X_max = std::max(h_test_signal->GetXaxis()->GetXmax(), h_test_signal->GetXaxis()->GetXmax());
    h_test_signal->GetXaxis()->SetRangeUser(X_min,X_max);
    //Loop on both cuts in [X_min;X_max]
    Int_t dim = 0;
    //Increase N to increase (a,b) scan granularity!
    double step = (X_max - X_min)/N;
    for(int i=0; i<N; i++){
	a = X_min + i * step;
	for(int j=0; j<N; j++){
	    b = X_min + j * step;
            if(a<b) continue;
	    for(int k=0; k<N; k++){
	        c = X_min + k * step;
                if(b<c) continue;
                //computing areas in range [a;X_max]
	        N_s_1 = TH1_integral(h_test_signal,a,X_max);
	        N_b_1 = TH1_integral(h_test_bkg,a,X_max);
                //skip iteration if integral in the tails is < 0.01% of total (sensitive to fluctuations!)
                if(N_s_1 < TH1_integral(h_test_signal,X_min,X_max)*0.0001) continue;
                if(N_b_1 < TH1_integral(h_test_bkg,X_min,X_max)*0.0001) continue;
                //computing areas in range [b;a]
	        N_s_2 = TH1_integral(h_test_signal,b,a);
                N_b_2 = TH1_integral(h_test_bkg,b,a);
                //skip iteration if integral in cat 2 is < 0.01% of total (sensitive to fluctuations!)
                if(N_s_2 < TH1_integral(h_test_signal,X_min,X_max)*0.0001) continue;
                if(N_b_2 < TH1_integral(h_test_bkg,X_min,X_max)*0.0001) continue;
                //computing areas in range [c;b]
	        N_s_3 = TH1_integral(h_test_signal,c,b);
                N_b_3 = TH1_integral(h_test_bkg,c,b);
                //skip iteration if integral in cat 3 is < 0.01% of total (sensitive to fluctuations!)
                if(N_s_3 < TH1_integral(h_test_signal,X_min,X_max)*0.0001) continue;
                if(N_b_3 < TH1_integral(h_test_bkg,X_min,X_max)*0.0001) continue;
                //scale n(b) depending on categ
                N_b_1 = N_b_1*bkg_scale;
                N_b_2 = N_b_2*bkg_scale;
                N_b_3 = N_b_3*bkg_scale;
 	        if ( (N_b_1)>0 && (N_b_2)>0 && (N_b_3)>0) {
                    double S = 0;
                    S1 = log_significance(N_s_1, N_b_1);
                    S2 = log_significance(N_s_2, N_b_2);
                    S3 = log_significance(N_s_3, N_b_3);
                    //S1 = N_s_1 / sqrt(N_s_1 + N_b_1);
	            //S2 = N_s_2 / sqrt(N_s_2 + N_b_2);
	            //Combined significance
	            S = sqrt(S1*S1 + S2*S2 + S3*S3);
	            h3->SetBinContent(i,j,k,S);
   	            dim++;
	        }
	    }
	}
    }
    //Taking absolute maximum of the combined significance

    Int_t nbinx=0, nbiny=0, nbinz=0;
    cout<<"h3->GetMaximumBin(nbinx, nbiny, nbinz) "<<endl;
    h3->GetMaximumBin(nbinx, nbiny, nbinz);
    cout<<"nbinx="<<nbinx<<", nbiny="<<nbiny<<", nbinz="<<nbinz<<endl;

    float bcx = ((TAxis*)h3->GetXaxis())->GetBinCenter(nbinx);
    float bcy = ((TAxis*)h3->GetYaxis())->GetBinCenter(nbiny);
    float bcz = ((TAxis*)h3->GetZaxis())->GetBinCenter(nbinz);
    cout<<"bcx="<<bcx<<" bcy="<<bcy<<" bcz="<<bcz<<endl;

    float S_max = h3->GetBinContent(h3->GetMaximumBin());
    cout<<"S_max="<<S_max<<endl;

    //Computing cut efficiency on signal
    Double_t N_S_12 = TH1_integral(h_test_signal,bcz,X_max);
    Double_t N_S_tot = TH1_integral(h_test_signal,X_min,X_max);
    cout<<"Signal events kept by BDT "<<N_S_12<<" over "<<N_S_tot<<" ratio: "<<N_S_12/N_S_tot<<endl;
    //Computing cut efficiency on backgroup
    Double_t N_B_12 = TH1_integral(h_test_bkg,bcz,X_max);
    Double_t N_B_tot = TH1_integral(h_test_bkg,X_min,X_max);
    cout<<"Background events kept by BDT "<<N_B_12<<" over "<<N_B_tot<<" ratio: "<<N_B_12/N_B_tot<<endl;
    
    return {bcx,bcy,bcz};
}



TGraph* plot_ROC(TString categ, TH1F *h_signal, TH1F *h_bkg){

    float bdt_min = std::min(h_signal->GetXaxis()->GetXmin(), h_bkg->GetXaxis()->GetXmin());
    float bdt_max = std::max(h_signal->GetXaxis()->GetXmax(), h_bkg->GetXaxis()->GetXmax());
    int n = 200; //number of points
    float step = abs(bdt_max-bdt_min)/(n*1.0);

    Double_t x[n], y[n];
    double n_true_positive = 0;
    double n_false_positive = 0;
    double n_true_negative_rej = 0;
    double n_false_negative_rej = 0;
    float bdt_wp = 0;

    for(int i = 0; i<n; i++){
        bdt_wp = bdt_min + step*i;

        n_true_positive  = TH1_integral(h_signal, bdt_wp, bdt_max);
        n_false_positive = TH1_integral(h_bkg, bdt_wp, bdt_max);
        n_true_negative_rej  = TH1_integral(h_bkg,    bdt_min, bdt_wp);
        n_false_negative_rej = TH1_integral(h_signal, bdt_min, bdt_wp);

        x[i] = n_true_positive/(n_true_positive+n_false_negative_rej);
        y[i] = n_true_negative_rej/(n_false_positive+n_true_negative_rej);
    }

    TGraph* gr = new TGraph(n,x,y);
    gr->SetTitle("ROC;Signal Efficiency;Background rejection");
    return gr;
}

array<double, 6> get_wp_TMVA(TChain *t, BDTcut3d BDTcutvalues){

   TString bdt_branch = "BDT";
   TString target_branch = "classID";
   TTreeReader reader (t);
   
   TTreeReaderValue<float> reader_bdt (reader, bdt_branch);
   TTreeReaderValue<int> reader_classid (reader, target_branch);

   double n_c_true_positive = 0;
   double n_c_false_positive = 0;
   double n_c_true_negative_rej = 0;
   double n_c_false_negative_rej = 0;
   double n_b_true_positive = 0;
   double n_b_false_positive = 0;
   double n_b_true_negative_rej = 0;
   double n_b_false_negative_rej = 0;
   double n_a_true_positive = 0;
   double n_a_false_positive = 0;
   double n_a_true_negative_rej = 0;
   double n_a_false_negative_rej = 0;

   while (reader.Next()) {
       if ((*reader_bdt) >= BDTcutvalues.c and (*reader_classid) == 0) {
           n_c_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.c and (*reader_classid) == 1){
           n_c_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.c and (*reader_classid) == 1) {
           n_c_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.c and (*reader_classid) == 0) {
           n_c_false_negative_rej += 1.;
       }

       if ((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 0) {
           n_b_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 1){
           n_b_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 1) {
           n_b_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 0) {
           n_b_false_negative_rej += 1.;
       }

       if ((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 0) {
           n_a_true_positive += 1.;
       } else if((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 1){
           n_a_false_positive += 1.;
       }
       if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 1) {
           n_a_true_negative_rej += 1.;
       } else if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 0) {
           n_a_false_negative_rej += 1.;
       }
    }

    array<double, 6> c;
    c[0] = n_c_true_positive/    (n_c_true_positive +n_c_false_negative_rej); //x[0]
    c[1] = n_c_true_negative_rej/(n_c_false_positive+n_c_true_negative_rej); //y[0]
    c[2] = n_b_true_positive/    (n_b_true_positive +n_b_false_negative_rej); //x[1]
    c[3] = n_b_true_negative_rej/(n_b_false_positive+n_b_true_negative_rej); //y[1]
    c[4] = n_a_true_positive/    (n_a_true_positive +n_a_false_negative_rej); //x[2]
    c[5] = n_a_true_negative_rej/(n_a_false_positive+n_a_true_negative_rej); //y[2]
    return c;
}

void plot_wp_ROC(TString categ, TH1F *h_test_signal, TH1F *h_test_bkg){
 
     BDTcut3d BDTcutvalues = Get_BDT_cut_3D(categ);

     TGraph* ROC =  plot_ROC(categ, h_test_signal, h_test_bkg);
     TCanvas *c1 = new TCanvas("c1","ROC "+categ,150,10,990,660);
 //    ROC->GetYaxis()->SetRangeUser(0.5,1.05);
     ROC->Draw();
 
     TChain *tTrain = new TChain("tree");
     tTrain->Add("TMVA_"+TMVA_inputpath+categ+".root/"+TMVA_inputpath+categ+"/TrainTree");

     tTrain->LoadTree(-1);
     int nentries = tTrain->GetEntries();
     cout<<"nentries "<<nentries<<endl;

     array<double, 6> coordinates;
     coordinates = get_wp_TMVA(tTrain, BDTcutvalues);

     Double_t x[3], y[3];
     x[0] = coordinates[0];
     y[0] = coordinates[1];
     x[1] = coordinates[2];
     y[1] = coordinates[3];
     x[2] = coordinates[4];
     y[2] = coordinates[5];

     TGraph* gr = new TGraph(3,x,y);
     TLatex *l = new TLatex(0.5, 0.5, "label");
     l->SetTextSize(0.025);
     l->SetTextFont(42);
     l->SetTextAlign(21);
     l->SetTextColor(kBlue);
     l->DrawLatex(x[0],y[0]+0.01,Form("bdt.c %4.3f",BDTcutvalues.c));
     l->DrawLatex(x[0],y[0]-0.02,Form("(%4.2f, %4.2f)",x[0], y[0]));
     l->DrawLatex(x[1],y[1]+0.01,Form("bdt.b %4.3f",BDTcutvalues.b));
     l->DrawLatex(x[1],y[1]-0.02,Form("(%4.2f, %4.2f)",x[1], y[1]));
     l->DrawLatex(x[2],y[2]+0.01,Form("bdt.a %4.3f",BDTcutvalues.a));
     l->DrawLatex(x[2],y[2]-0.02,Form("(%4.2f, %4.2f)",x[2], y[2]));
     gr->Draw("same *p");

     c1->Update();
     c1->SaveAs(TMVA_inputpath+categ+"/"+method+"_ROC_"+categ+".png");
     delete tTrain;
}
#endif // BDT_OPTIMAL_CUT_H ///:~
