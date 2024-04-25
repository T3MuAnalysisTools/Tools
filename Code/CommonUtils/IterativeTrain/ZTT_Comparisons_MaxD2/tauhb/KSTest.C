#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"

void KSTest () 
{
        
        //cat_type = "_ZTT_tau_NoCV_3mu";
        //cat_type = "_ZTT_tau_CV_3mu"
        //cat_type = "_ZTT_mu3mu";
        //cat_type = "_ZTT_e3mu"
        
        TString cat_type = "_ZTT_tau_CV_3mu";
            
        const int no_cats(18+1);
        
        TString cat_label[4];
        cat_label[0] = "{h,A}";
        cat_label[1] = "{h,B}";
        cat_label[2] = "{#mu}";
        cat_label[3] = "{e}";
        
        TString cat_l=cat_label[1];
        
        TFile *TreeFile[no_cats];
        TString hname;
        
        TTree *TreeTest[no_cats];
        
        TH1D  * Train_Sig[no_cats];
        TH1D  * Train_Bkg[no_cats];
        TH1D  * Test_Sig[no_cats];
        TH1D  * Test_Bkg[no_cats];
        
        double ks_score_sig[no_cats];
        double ks_score_bkg[no_cats];
        
        size_t sigBins, bkgBins;
        Double_t sigLeft, sigRight, bkgLeft, bkgRight;
        
        
        
        for(int i=0; i<no_cats; i++){
            
            hname=to_string(i);
            TreeFile[i] = new TFile("category_"+hname+cat_type+".root","READ");
            
            TreeTest[i] = (TTree *) TreeFile[i]->Get("output_"+hname+cat_type+"/TestTree");
            
            Train_Sig[i] = (TH1D *)TreeFile[i]->Get("output_"+hname+cat_type+"/"+ "Method_BDT/BDT/MVA_BDT_Train_S");
            Train_Bkg[i] = (TH1D *)TreeFile[i]->Get("output_"+hname+cat_type+"/"+ "Method_BDT/BDT/MVA_BDT_Train_B");
            
            sigBins = Train_Sig[i]->GetNbinsX();
            sigLeft = Train_Sig[i]->GetBinLowEdge(1);
            sigRight = Train_Sig[i]->GetBinLowEdge(sigBins + 1);
            bkgBins = Train_Bkg[i]->GetNbinsX();
            bkgLeft = Train_Bkg[i]->GetBinLowEdge(1);
            bkgRight = Train_Bkg[i]->GetBinLowEdge(bkgBins + 1);
            
            Test_Sig[i] = new TH1D("Test_Sig"+hname, "Test_Sig", sigBins, sigLeft, sigRight);
            Test_Bkg[i] = new TH1D("Test_Bkg"+hname, "Test_Bkg", bkgBins, bkgLeft, bkgRight);
            TreeTest[i]->Project("Test_Sig"+hname, "BDT", "classID==0");
            TreeTest[i]->Project("Test_Bkg"+hname, "BDT", "classID==1");
            
            ks_score_sig[i]=Train_Sig[i]->KolmogorovTest(Test_Sig[i],"X");
            ks_score_bkg[i]=Train_Bkg[i]->KolmogorovTest(Test_Bkg[i],"X");
            
            cout << "Sig.: " << ks_score_sig[i] << '\n';
            cout << "Bkg.: " << ks_score_bkg[i] << endl;
            
            if(ks_score_sig[i]>0.20&&ks_score_sig[i]<0.80&&ks_score_bkg[i]>0.20&&ks_score_bkg[i]<0.80){
                    cout << "Use: " << i << " to plot ROCs." << endl;
            }
        
        }
        
        
        double cat_list[no_cats];
        
        for (Int_t i = 0; i < no_cats; ++i) {
                cat_list[i] = no_cats+5-i-1;// No of vars
        }
        
        TGraph *graph1 = new TGraph(no_cats, cat_list, ks_score_sig);
        graph1->SetMarkerStyle(20); // Set marker style to a filled circle
        graph1->SetMarkerSize(1.5); // Set marker size
        graph1->SetMarkerColor(kBlue);
        graph1->SetTitle("KS Test Scores, #tau_"+cat_l+";No of variables;KS Test Score");
        
        TGraph *graph2 = new TGraph(no_cats, cat_list, ks_score_bkg);
        graph2->SetMarkerStyle(21); // Set marker style to a filled square
        graph2->SetMarkerSize(1.5); // Set marker size
        graph2->SetMarkerColor(kRed);
        
        graph1->GetYaxis()->SetRangeUser(0, 1.0); // Set the range for graph1
        graph2->GetYaxis()->SetRangeUser(0, 1.0); // Set the range for graph2
        
        TCanvas *canvas = new TCanvas("canvas", "Multiple Graphs", 800, 600);
        graph1->Draw("AP");
        graph2->Draw("P SAME");
        
        TLine *ref_line = new TLine(graph1->GetXaxis()->GetXmin(), 0.5, graph1->GetXaxis()->GetXmax(), 0.5);
        ref_line->SetLineColor(kBlack); // Set line color
        ref_line->SetLineStyle(2); // Set line style (dashed)
        ref_line->Draw(); // Draw the line on the canvas
        
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Specify legend position (x1, y1, x2, y2)
        // Add entries to the legend
        legend->AddEntry(graph1, "S", "lp"); // "lp" option indicates line and marker for the entry
        legend->AddEntry(graph2, "B", "lp");
        //legend->AddEntry(ref_line, "Ref", "lp");
        // Draw the legend
        legend->Draw();
        
        //canvas->Draw();
        canvas->SaveAs("tauhb.png");
        delete canvas;
    




}


