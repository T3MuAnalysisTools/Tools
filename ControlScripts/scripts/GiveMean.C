void GiveMean()
{
  TFile *file0=new TFile("LOCAL_COMBINED_signalvertexselector_default_LumiScaled.root");
  
  
  // MC1: B0, MC2: B+-  (i =0 and i=1)
  
  TH1F  * hists[2][10]; // first is MC type, second is number of histos
  TString hname;
  
  
  for(int i=0; i<2;i++){
    hname=to_string(i+1);
    
    hists[i][0] = (TH1F*)file0->Get("signalvertexselector_default_PairVtxPV_dRto_SVPV_2CorMC"+hname);
    
    hists[i][1] = (TH1F*)file0->Get("signalvertexselector_default_PairVtxPV_dRto_PairLV_2CorMC"+hname);
    
    hists[i][2] = (TH1F*)file0->Get("signalvertexselector_default_PairVtxPV_SV_Significance_2CorMC"+hname);
    
    hists[i][3] = (TH1F*)file0->Get("signalvertexselector_default_PairAngle_2CorMC"+hname);
    
    hists[i][4] = (TH1F*)file0->Get("signalvertexselector_default_PairVtxDist_2CorMC"+hname);
    
    hists[i][5] = (TH1F*)file0->Get("signalvertexselector_default_PairProduct_2CorMC"+hname);
    
    //hists[i][6] = (TH1F*)file0->Get("signalvertexselector_default_Rank_Correct_2iso_Chi2MC"+hname);
    
    //hists[i][7] = (TH1F*)file0->Get("signalvertexselector_default_Rank_Correct_3iso_Chi2MC"+hname);
  }
  
  std::cout << "Mean of PairVtxPV_dRto_SVPV_2Cor for MC1/B0: " << hists[0][0]->GetXaxis()->GetBinCenter(hists[0][0]->GetMaximumBin()) << " with RMS: " << hists[0][0]->GetStdDev() << " Mean of PairVtxPV_dRto_SVPV_2Cor for MC2/B+-: " << hists[1][0]->GetXaxis()->GetBinCenter(hists[1][0]->GetMaximumBin()) << " with RMS: " << hists[1][0]->GetStdDev() << std::endl;
  
  std::cout << "Mean of PairVtxPV_dRto_PairLV_2Cor for MC1/B0: " << hists[0][1]->GetXaxis()->GetBinCenter(hists[0][1]->GetMaximumBin()) << " with RMS: " << hists[0][1]->GetStdDev() << " Mean of PairVtxPV_dRto_PairLV_2Cor for MC2/B+-: " << hists[1][1]->GetXaxis()->GetBinCenter(hists[1][1]->GetMaximumBin()) << " with RMS: " << hists[1][1]->GetStdDev() << std::endl;
  
  std::cout << "Mean of PairVtxPV_SV_Significance_2Cor for MC1/B0: " << hists[0][2]->GetXaxis()->GetBinCenter(hists[0][2]->GetMaximumBin()) << " with RMS: " << hists[0][2]->GetStdDev() << " Mean of PairVtxPV_SV_Significance_2Cor for MC2/B+-: " << hists[1][2]->GetXaxis()->GetBinCenter(hists[1][2]->GetMaximumBin()) << " with RMS: " << hists[1][2]->GetStdDev() << std::endl;
  
  std::cout << "Mean of PairAngle_2Cor for MC1/B0: " << hists[0][3]->GetXaxis()->GetBinCenter(hists[0][3]->GetMaximumBin()) << " with RMS: " << hists[0][3]->GetStdDev() << " Mean of PairAngle_2Cor for MC2/B+-: " << hists[1][3]->GetXaxis()->GetBinCenter(hists[1][3]->GetMaximumBin()) << " with RMS: " << hists[1][3]->GetStdDev() << std::endl;
  
  std::cout << "Mean of PairVtxDist_2CorMC for MC1/B0: " << hists[0][4]->GetXaxis()->GetBinCenter(hists[0][4]->GetMaximumBin()) << " with RMS: " << hists[0][4]->GetStdDev() << " Mean of PairVtxDist_2CorMC for MC2/B+-: " << hists[1][4]->GetXaxis()->GetBinCenter(hists[1][4]->GetMaximumBin()) << " with RMS: " << hists[1][4]->GetStdDev() << std::endl;
  
  std::cout << "Mean of PairProduct_2CorMC for MC1/B0: " << hists[0][5]->GetXaxis()->GetBinCenter(hists[0][5]->GetMaximumBin()) << " with RMS: " << hists[0][5]->GetStdDev() << " Mean of PairProduct_2CorMC for MC2/B+-: " << hists[1][5]->GetXaxis()->GetBinCenter(hists[1][5]->GetMaximumBin()) << " with RMS: " << hists[1][5]->GetStdDev() << std::endl;
  
}