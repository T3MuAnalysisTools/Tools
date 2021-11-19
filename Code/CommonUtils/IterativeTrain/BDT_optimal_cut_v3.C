#include "TH1F.h"
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "T3M_common.h"
#include "BDT_optimal_cut_v3.h"

using namespace std;

void BDT_optimal_cut_v3() 
{

    int ncat = sizeof(cat_label) / sizeof(cat_label[0]);

    for(int k=0; k<ncat; k++)
    {
        cout<<"category "<<cat_label[k]<<endl;
        TH1F *h_test_signal;
        TH1F *h_test_bkg;
        TH1F *h_test_signal2;
        TH1F *h_test_bkg2;

        TString file_name = TMVA_inputpath+"outputTree.root";

        //open input files
        TChain *t = new TChain("outputTree");
        t->Add(file_name); 
        std::cout<<"Opened input file: "<<file_name<<std::endl;

        TString c = std::to_string(k); 
        TString signal = "weight*(isMC>0  && category =="+c+")";
        TString bkg    = "weight*(isMC==0 && category =="+c+") && isSB == 1";
        
        //bdt score distribution
        TString binning = "(240, -0.6, 0.6)"; 
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
       
        //     h_test_signal->GetXaxis()->SetRangeUser(0,1.0);

        //Make up on plots
        h_test_bkg->SetLineColor(kBlack);
        h_test_signal->SetLineColor(kRed);

        double X_min = std::min(h_test_signal->GetXaxis()->GetXmin(), h_test_signal->GetXaxis()->GetXmin());
        double X_max = std::max(h_test_signal->GetXaxis()->GetXmax(), h_test_signal->GetXaxis()->GetXmax());
       
        X_min = -0.4; X_max = 0.4;
        //Compute cut and make colz plot
        BDTcut3d cut_value = Get_BDT_cut_3D(cat_label[k]);
        //Draw BDT working point on ROC curve
        plot_wp_ROC(cat_label[k], h_test_signal, h_test_bkg);

        TCanvas *c1 = new TCanvas("c1","c1",150,10,800,800);
        h_test_bkg->Draw("HISTE");
        h_test_signal->Draw("same HISTE");
        c1->Update();
        TLine l;
        l.DrawLine(cut_value.a,0,cut_value.a,0.1);
        l.DrawLine(cut_value.b,0,cut_value.b,0.1);
        l.DrawLine(cut_value.c,0,cut_value.c,0.1);

        TLegend*leg = new TLegend(0.1,0.75,0.45,0.9);
        leg->AddEntry(h_test_signal,cat_label[k]+"_signal","f");
        leg->AddEntry(h_test_bkg,cat_label[k]+"_bkg","f");
        leg->Draw();
        c1->Update();
        c1->SaveAs(TMVA_inputpath+cat_label[k]+"/"+method+"_"+cat_label[k]+"normalizedSignal_v3.png");

        //Drawing BDT score from scratch without signal normalization
        TCanvas *c2 = new TCanvas("c2","c2",150,10,800,800);
        h_test_signal2->GetXaxis()->SetRangeUser(X_min,X_max);
        h_test_bkg2->GetXaxis()->SetRangeUser(X_min,X_max);
        h_test_bkg2->SetLineColor(kBlack);
        h_test_signal2->SetLineColor(kRed);
        h_test_signal2->Scale(1/h_test_signal2->Integral());
        h_test_bkg2->Scale(1/h_test_bkg2->Integral());
        h_test_signal2->Rebin(4);
        h_test_bkg2->Rebin(4);

        double Y_max;
        Y_max = 0.12;

        h_test_bkg2->Draw("HISTE");
        h_test_bkg2->GetYaxis()->SetRangeUser(0, Y_max);
        h_test_bkg2->GetXaxis()->SetRangeUser(X_min,X_max);
        h_test_signal2->Draw("same HISTE");
        c2->Update();
        l.DrawLine(cut_value.a,0,cut_value.a,Y_max);
        l.DrawLine(cut_value.b,0,cut_value.b,Y_max);
        l.DrawLine(cut_value.c,0,cut_value.c,Y_max);
        TLatex ta(cut_value.a,0,"a");
        ta.Draw();
        TLatex tb(cut_value.b,0,"b");
        tb.Draw();
        TLatex tc(cut_value.c,0,"c");
        tc.Draw();

        TLegend*leg2 = new TLegend(0.1,0.76,0.45,0.9);
        leg2->AddEntry(h_test_signal2,cat_label[k]+"_signal","f");
        leg2->AddEntry(h_test_bkg2,cat_label[k]+"_bkg","f");
        leg2->Draw();
        c2->Update();
        c2->SaveAs(TMVA_inputpath+cat_label[k]+"/"+method+"_"+cat_label[k]+"_v3.png");

     }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}
