#include <TCanvas.h>
#include <TH2F.h>
#include <TROOT.h>
#include <iostream>
#include <TSystem.h>
#include <TLatex.h>

bool CreateDirectoryIfNotExists(const char* path) {
    if (gSystem->AccessPathName(path) == 0) {
        std::cout << "Directory already exists: " << path << std::endl;
        return true;
    } else {
        if (gSystem->mkdir(path) == 0) {
            std::cout << "Directory created successfully: " << path << std::endl;
            return true;
        } else {
            std::cerr << "Error: Failed to create directory: " << path << std::endl;
            return false;
        }
    }
}

void ResizeCanvasAndHistogram(const char* canvasName, const char* histogramName, Int_t width, Int_t height, Int_t cat) {
    
    TString cat_label1[4];
    cat_label1[0] = "{h,A}";
    cat_label1[1] = "{h,B}";
    cat_label1[2] = "{#mu}";
    cat_label1[3] = "{e}";
    
    // Attempt to find the canvas and histogram objects
    TObject* obj = gROOT->FindObject(canvasName);
    TCanvas* canvas = dynamic_cast<TCanvas*>(obj);
    if (!canvas) {
        std::cerr << "Error: Canvas not found." << std::endl;
        return;
    }

    TH2F* histogram = dynamic_cast<TH2F*>(canvas->GetPrimitive(histogramName));
    if (!histogram) {
        std::cerr << "Error: Histogram not found." << std::endl;
        return;
    }

    
    // Replace the original histogram with the new one
    canvas->Clear();
    
    // Resize the canvas
    canvas->SetCanvasSize(width, height);
    
    histogram->Draw("colz");
    canvas->Update();
    
    TPaletteAxis* paletteAxis = (TPaletteAxis*)histogram->GetListOfFunctions()->FindObject( "palette" );
    paletteAxis->SetLabelSize( 0.03 );
    paletteAxis->SetX1NDC( paletteAxis->GetX1NDC() + 0.02 );
    
    histogram->Draw("textsame");  // add text
    
    TLatex* t = new TLatex( 0.4, 0.89, "Linear correlation coefficients in %, Category #tau_"+cat_label1[cat] );
    t->SetNDC();
    t->SetTextSize( 0.026 );
    t->AppendPad();  
    
    canvas->Update();

    std::cout << "Canvas and histogram resized successfully!" << std::endl;
}

void MakeCorrelation () 
{
        
        TString cat_label[4];
        cat_label[0] = "correlation/tauha";
        cat_label[1] = "correlation/tauhb";
        cat_label[2] = "correlation/taumu";
        cat_label[3] = "correlation/taue";
        
        TString hname;
        
        TString datasetname;
        
        TString cat_type[4];
        cat_type[0] = "_ZTT_tau_NoCV_3mu";
        cat_type[1] = "_ZTT_tau_CV_3mu";
        cat_type[2] = "_ZTT_mu3mu";
        cat_type[3] = "_ZTT_e3mu";
        
        double left_margin[4];
        left_margin[0] = 0.23;
        left_margin[1] = 0.23;
        left_margin[2] = 0.27;
        left_margin[3] = 0.23;
        
        double right_margin[4];
        right_margin[0] = 0.11;
        right_margin[1] = 0.11;
        right_margin[2] = 0.22;
        right_margin[3] = 0.16;
        
        double bottom_margin[4];
        bottom_margin[0] = 0.15;
        bottom_margin[1] = 0.15;
        bottom_margin[2] = 0.17;
        bottom_margin[3] = 0.15;
        
        Bool_t isRegression=kFALSE; // If classification
        
        const int no_cats(4);
        
        for(int i=0; i<no_cats; i++){
        
            hname=to_string(0);
        
            TMVA::correlations("output_"+hname+cat_type[i], "category_"+hname+cat_type[i]+".root", isRegression, kFALSE, kTRUE);
        
            TCanvas & c1 = *dynamic_cast<TCanvas *>(gROOT->FindObject("CorrelationMatrixS"));
            c1.SetLeftMargin( left_margin[i] );
            c1.SetRightMargin( right_margin[i] );
            c1.SetBottomMargin( bottom_margin[i] );
            ResizeCanvasAndHistogram("CorrelationMatrixS", "CorrelationMatrixS", 1500, 1000, i);
        
            if (CreateDirectoryIfNotExists(cat_label[i])) {
                c1.SaveAs(cat_label[i]+"/CorrelationMatrixS.png");
            }
        
            TCanvas & c2 = *dynamic_cast<TCanvas *>(gROOT->FindObject("CorrelationMatrixB"));
            c2.SetLeftMargin( left_margin[i] );
            c2.SetRightMargin( right_margin[i] );
            c2.SetBottomMargin( bottom_margin[i] );
            ResizeCanvasAndHistogram("CorrelationMatrixB", "CorrelationMatrixB", 1500, 1000, i);
        
            if (CreateDirectoryIfNotExists(cat_label[i])) {
                c2.SaveAs(cat_label[i]+"/CorrelationMatrixB.png");
            }
        
        
        }



}


