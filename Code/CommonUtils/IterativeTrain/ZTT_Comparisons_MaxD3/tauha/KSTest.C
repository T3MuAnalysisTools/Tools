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
        
        TString cat_type = "_ZTT_tau_NoCV_3mu";
            
        const int no_cats(10+1);
        
        TString cat_label[4];
        cat_label[0] = "{h,A}";
        cat_label[1] = "{h,B}";
        cat_label[2] = "{#mu}";
        cat_label[3] = "{e}";
        
        TString cat_l=cat_label[0];
        
        TFile *TreeFile[no_cats];
        TString hname;
        
        TTree *TreeTest[no_cats];
        TTree *TreeTrain[no_cats];
        
        TH1D  * Train_Sig[no_cats];
        TH1D  * Train_Bkg[no_cats];
        TH1D  * Test_Sig[no_cats];
        TH1D  * Test_Bkg[no_cats];
        
        TH1D  * Train_Bkg_Restricted[no_cats];//Look only at the tail of the 
        TH1D  * Test_Bkg_Restricted[no_cats];
        
        double ks_score_sig[no_cats];
        double ks_score_bkg[no_cats];
        double ks_score_bkg_restr[no_cats];
        
        size_t sigBins, bkgBins, bkgBins_restricted, bkg_restricted_startBin, bkg_restricted_endBin;
        Double_t sigLeft, sigRight, bkgLeft, bkgRight, bkgLeft_restricted, bkgRight_restricted;
        
        vector<double> bdt_vals_sig_train[no_cats];
        vector<double> bdt_vals_sig_test[no_cats];
        vector<double> bdt_vals_bkg_train[no_cats];
        vector<double> bdt_vals_bkg_test[no_cats];
        vector<double> bdt_vals_bkg_train_tail[no_cats];
        vector<double> bdt_vals_bkg_test_tail[no_cats];
        
        
        double ks_score_sig_unbinned[no_cats];
        double ks_score_bkg_unbinned[no_cats];
        double ks_score_bkg_restr_unbinned[no_cats];
        
        
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
            
            bkg_restricted_startBin=Test_Bkg[i]->FindBin(-0.0);
            bkgLeft_restricted=Test_Bkg[i]->GetBinLowEdge(bkg_restricted_startBin);
            bkg_restricted_endBin=bkgBins;
            
            for (int j = 1; j <= Test_Bkg[i]->GetNbinsX(); ++j) {//end bin: largest occupied bin
                    if(Train_Bkg[i]->GetBinContent(j)>0.00001 ||  Test_Bkg[i]->GetBinContent(j)>0.00001 ){
                            bkg_restricted_endBin=j;
                    }
            }
            
            
            bkgRight_restricted=Test_Bkg[i]->GetBinLowEdge(bkg_restricted_endBin+1);
            bkgBins_restricted=bkg_restricted_endBin-bkg_restricted_startBin+1;
            
            //cout << "Bkg bins: " << bkgBins << " Bkg bins restricted: " << bkgBins_restricted << endl;
            //cout << "Bkg left: " << bkgLeft << " Bkg left restricted: " << bkgLeft_restricted << " Bkg right: " << bkgRight << " Bkg right restricted: " << bkgRight_restricted << endl;
            
            Train_Bkg_Restricted[i] = new TH1D("Train_Bkg_Restricted"+hname, "Train_Bkg_Restricted", bkgBins_restricted, bkgLeft_restricted, bkgRight_restricted);
            Test_Bkg_Restricted[i] = new TH1D("Test_Bkg_Restricted"+hname, "Test_Bkg_Restricted", bkgBins_restricted, bkgLeft_restricted, bkgRight_restricted);
            
            double train_integral(0.0);
            double train_integral_tail(0.0);
            for (int j = 1; j <= Train_Bkg[i]->GetNbinsX(); ++j) {
                    
                    train_integral+=Train_Bkg[i]->GetBinContent(j);
                    
                    //cout << "Training bin content for bin: " << j << " is: " << Train_Bkg[i]->GetBinContent(j) << " with error percent: " << 100*Train_Bkg[i]->GetBinError(j)/(Train_Bkg[i]->GetBinContent(j)+0.0000001) << endl;
                    
                    if(j>=bkg_restricted_startBin&&j<=bkg_restricted_endBin){
                            
                            train_integral_tail+=Train_Bkg[i]->GetBinContent(j);
                            
                            Train_Bkg_Restricted[i]->SetBinContent(j-bkg_restricted_startBin+1, Train_Bkg[i]->GetBinContent(j));
                            Train_Bkg_Restricted[i]->SetBinError(j-bkg_restricted_startBin+1, Train_Bkg[i]->GetBinError(j));
                            Test_Bkg_Restricted[i]->SetBinContent(j-bkg_restricted_startBin+1, Test_Bkg[i]->GetBinContent(j));
                            Test_Bkg_Restricted[i]->SetBinError(j-bkg_restricted_startBin+1, Test_Bkg[i]->GetBinError(j));
                    }
            }
            
            ks_score_sig[i]=Train_Sig[i]->KolmogorovTest(Test_Sig[i],"X=10000");
            ks_score_bkg[i]=Train_Bkg[i]->KolmogorovTest(Test_Bkg[i],"X=10000");
            ks_score_bkg_restr[i]=Train_Bkg_Restricted[i]->KolmogorovTest(Test_Bkg_Restricted[i],"X=10000");
            
            
            //More accurate K-S Test (unbinned)
            
            TreeTrain[i] = (TTree *) TreeFile[i]->Get("output_"+hname+cat_type+"/TrainTree");
            
            Float_t bdt;  // Variable to hold the BDT value for each entry
            Int_t classID; // Variable to hold the classID value for each entry
        
            // Set branch addresses to read the BDT and classID branches
            TreeTest[i]->SetBranchAddress("BDT", &bdt);
            TreeTest[i]->SetBranchAddress("classID", &classID);
            
            TreeTrain[i]->SetBranchAddress("BDT", &bdt);
            TreeTrain[i]->SetBranchAddress("classID", &classID);
        
            // Loop over all entries in the TTree
            Long64_t nentries2 = TreeTrain[i]->GetEntries();
            for (Long64_t j = 0; j < nentries2; ++j) {
                TreeTrain[i]->GetEntry(j);  // Get the j-th entry
                
                // classID==0 for signal, classID==1 for background
                if (classID == 0) {
                    //cout<< "Train BDT Signal: " << bdt << endl;
                    bdt_vals_sig_train[i].push_back(bdt);
                }
                if (classID == 1) {
                    //cout<< "Train BDT Background: " << bdt << endl;
                    bdt_vals_bkg_train[i].push_back(bdt);
                    if(bdt>0.0){
                        bdt_vals_bkg_train_tail[i].push_back(bdt);
                    }
                }
            }
            
            // Loop over all entries in the TTree
            Long64_t nentries1 = TreeTest[i]->GetEntries();
            for (Long64_t j = 0; j < nentries1; ++j) {
                TreeTest[i]->GetEntry(j);  // Get the j-th entry
                
                // classID==0 for signal, classID==1 for background
                if (classID == 0) {
                    //cout<< "Test BDT Signal: " << bdt << endl;
                    bdt_vals_sig_test[i].push_back(bdt);
                    
                }
                if (classID == 1) {
                    //cout<< "Test BDT Background: " << bdt << endl;
                    bdt_vals_bkg_test[i].push_back(bdt);
                    if(bdt>0.0){
                        bdt_vals_bkg_test_tail[i].push_back(bdt);
                    }
                }
            }
            
            
            //Sort them and copy them to an array
            sort(bdt_vals_sig_train[i].begin(), bdt_vals_sig_train[i].end());
            sort(bdt_vals_bkg_train[i].begin(), bdt_vals_bkg_train[i].end());
            sort(bdt_vals_bkg_train_tail[i].begin(), bdt_vals_bkg_train_tail[i].end());
            sort(bdt_vals_sig_test[i].begin(), bdt_vals_sig_test[i].end());
            sort(bdt_vals_bkg_test[i].begin(), bdt_vals_bkg_test[i].end());
            sort(bdt_vals_bkg_test_tail[i].begin(), bdt_vals_bkg_test_tail[i].end());
            
            double bdt_vals_sig_train_array[bdt_vals_sig_train[i].size()];
            double bdt_vals_bkg_train_array[bdt_vals_bkg_train[i].size()];
            double bdt_vals_bkg_train_tail_array[bdt_vals_bkg_train_tail[i].size()];
            double bdt_vals_sig_test_array[bdt_vals_sig_test[i].size()];
            double bdt_vals_bkg_test_array[bdt_vals_bkg_test[i].size()];
            double bdt_vals_bkg_test_tail_array[bdt_vals_bkg_test_tail[i].size()];
            
            copy(bdt_vals_sig_train[i].begin(), bdt_vals_sig_train[i].end(), bdt_vals_sig_train_array);
            copy(bdt_vals_bkg_train[i].begin(), bdt_vals_bkg_train[i].end(), bdt_vals_bkg_train_array);
            copy(bdt_vals_bkg_train_tail[i].begin(), bdt_vals_bkg_train_tail[i].end(), bdt_vals_bkg_train_tail_array);
            copy(bdt_vals_sig_test[i].begin(), bdt_vals_sig_test[i].end(), bdt_vals_sig_test_array);
            copy(bdt_vals_bkg_test[i].begin(), bdt_vals_bkg_test[i].end(), bdt_vals_bkg_test_array);
            copy(bdt_vals_bkg_test_tail[i].begin(), bdt_vals_bkg_test_tail[i].end(), bdt_vals_bkg_test_tail_array);
            
            
            ks_score_sig_unbinned[i]=TMath::KolmogorovTest(bdt_vals_sig_train[i].size(),bdt_vals_sig_train_array,bdt_vals_sig_test[i].size(),bdt_vals_sig_test_array,"");
            ks_score_bkg_unbinned[i]=TMath::KolmogorovTest(bdt_vals_bkg_train[i].size(),bdt_vals_bkg_train_array,bdt_vals_bkg_test[i].size(),bdt_vals_bkg_test_array,"");
            ks_score_bkg_restr_unbinned[i]=TMath::KolmogorovTest(bdt_vals_bkg_train_tail[i].size(),bdt_vals_bkg_train_tail_array,bdt_vals_bkg_test_tail[i].size(),bdt_vals_bkg_test_tail_array,"");
            
            cout << "i : " << i << endl;
            //These are binned
            
            //cout << "Train integral: " << train_integral << " train integral ratio: : " << train_integral_tail/train_integral << endl;
            //cout << "Sig.: " << ks_score_sig[i] << '\n';
            //cout << "Bkg.: " << ks_score_bkg[i] << endl;
            //cout << "Bkg Restr.: " << ks_score_bkg_restr[i] << endl;
            
            //cout << "Bkg training events: " << Train_Bkg[i]->GetEntries() << " Bkg training tail events: " << Train_Bkg_Restricted[i]->GetEntries() << endl;
            
            
            cout<< "K-S Probabibility Signal Unbinned: " << ks_score_sig_unbinned[i] << endl;
            cout<< "K-S Probabibility Background Unbinned: " << ks_score_bkg_unbinned[i] << endl;
            cout<< "K-S Probabibility Background Tail Unbinned: " << ks_score_bkg_restr_unbinned[i] << endl;
            
            //if(ks_score_sig[i]>0.20&&ks_score_sig[i]<0.80 && ks_score_bkg[i]>0.20&&ks_score_bkg[i]<0.80 && ks_score_bkg_restr[i]>0.40&&ks_score_bkg_restr[i]<0.60){
            if(ks_score_sig_unbinned[i]>0.2 && ks_score_bkg_unbinned[i]>0.2 && ks_score_bkg_restr_unbinned[i]>0.40){
                    cout << "Use: " << i << " to plot ROC curves." << endl;
            }
        
        }
        
        
        double cat_list[no_cats];
        
        for (Int_t i = 0; i < no_cats; ++i) {
                cat_list[i] = no_cats+5-i-1;// No of vars
        }
        
        TGraph *graph1 = new TGraph(no_cats, cat_list, ks_score_sig_unbinned);
        graph1->SetMarkerStyle(20); // Set marker style to a filled circle
        graph1->SetMarkerSize(1.5); // Set marker size
        graph1->SetMarkerColor(kBlue);
        graph1->SetTitle("KS Test Scores, #tau_"+cat_l+";Number of variables;Kolmogorov test score");
        
        TGraph *graph2 = new TGraph(no_cats, cat_list, ks_score_bkg_unbinned);
        graph2->SetMarkerStyle(21); // Set marker style to a filled square
        graph2->SetMarkerSize(1.5); // Set marker size
        graph2->SetMarkerColor(kRed);
        
        TGraph *graph3 = new TGraph(no_cats, cat_list, ks_score_bkg_restr_unbinned);
        graph3->SetMarkerStyle(22); // Set marker style to a filled square
        graph3->SetMarkerSize(1.5); // Set marker size
        graph3->SetMarkerColor(kGreen);
        
        graph1->GetYaxis()->SetRangeUser(0, 1.0); // Set the range for graph1
        graph2->GetYaxis()->SetRangeUser(0, 1.0); // Set the range for graph2
        graph3->GetYaxis()->SetRangeUser(0, 1.0); // Set the range for graph2
        
        TCanvas *canvas = new TCanvas("canvas", "Multiple Graphs", 800, 600);
        graph1->Draw("AP");
        graph2->Draw("P SAME");
        graph3->Draw("P SAME");
        
        TLine *ref_line = new TLine(graph1->GetXaxis()->GetXmin(), 0.5, graph1->GetXaxis()->GetXmax(), 0.5);
        ref_line->SetLineColor(kBlack); // Set line color
        ref_line->SetLineStyle(2); // Set line style (dashed)
        //ref_line->Draw(); // Draw the line on the canvas
        
        TLegend *legend = new TLegend(0.81, 0.4, 0.9, 0.6); // Specify legend position (x1, y1, x2, y2)
        // Add entries to the legend
        legend->AddEntry(graph1, "S", "lp"); // "lp" option indicates line and marker for the entry
        legend->AddEntry(graph2, "B", "lp");
        legend->AddEntry(graph3, "B_{tail}", "lp");
        //legend->AddEntry(ref_line, "Ref", "lp");
        // Draw the legend
        legend->Draw();
        
        //canvas->Draw();
        canvas->SaveAs("tauha.png");
        delete canvas;
        
    




}


