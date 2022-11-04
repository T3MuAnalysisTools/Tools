void GiveNum()
{
  TFile *file0=new TFile("LOCAL_COMBINED_signalvertexselector_default_LumiScaled.root");
  
  
  // MC1: B0, MC2: B+-  (i =0 and i=1)
  
  TH1F  * hists[2][10]; // first is MC type, second is number of histos
  float WhetherdRMatch[2][2]; // first is MC type, second is bin no
  float WhetherSVMassMatch1[2][2]; // first is MC type, second is bin no
  float WhetherSVThreeMassMatch1[2][2]; // first is MC type, second is bin no
  TString hname;
  
  float Percentages[2][3];// No of matched tracks in iso. 0: less than 2, 1: greater than or equal to two, 2: greater than or equal to three
  float ChargedParticleEvents[2];//percent of events where there are atleast one charged particle
  
  
  for(int i=0; i<2;i++){
    hname=to_string(i+1);
    
    hists[i][0] = (TH1F*)file0->Get("signalvertexselector_default_WhetherdRMatchMC"+hname);
    
    hists[i][1] = (TH1F*)file0->Get("signalvertexselector_default_IsoTrackMatchedToSV_MassMatch1MC"+hname);
    
    hists[i][2] = (TH1F*)file0->Get("signalvertexselector_default_IsoTrackMatchedToSV_ThreeMassMatch1MC"+hname);
    
    hists[i][3] = (TH1F*)file0->Get("signalvertexselector_default_NumberOfFS_ChargedParticles_RecoMatchMC"+hname);
    
    hists[i][4] = (TH1F*)file0->Get("signalvertexselector_default_NumberOfFS_ChargedParticlesMC"+hname);
    
    hists[i][5] = (TH1F*)file0->Get("signalvertexselector_default_Rank_Correct_1iso3mu_Chi2MC"+hname);
    
    hists[i][6] = (TH1F*)file0->Get("signalvertexselector_default_Rank_Correct_2iso_Chi2MC"+hname);
    
    hists[i][7] = (TH1F*)file0->Get("signalvertexselector_default_Rank_Correct_3iso_Chi2MC"+hname);
    
    hists[i][8] = (TH1F*)file0->Get("signalvertexselector_default_Whether_Lowest_Chi2_is_Correct_2iso_varProduct1MC"+hname);
  }
  
  for(int i=0; i<2;i++){
    for(int j=0; j<2;j++){
      WhetherdRMatch[i][j]=hists[i][0]->GetBinContent(j+1);
      WhetherSVMassMatch1[i][j]=hists[i][1]->GetBinContent(j+1);
      WhetherSVThreeMassMatch1[i][j]=hists[i][2]->GetBinContent(j+1);
      //std::cout << "Bin content is : " <<hists[i][2]->GetBinContent(j+1)<< std::endl;
    }
    
    Percentages[i][0]=(hists[i][3]->GetBinContent(1)+hists[i][3]->GetBinContent(2))/(hists[i][3]->GetBinContent(1)+hists[i][3]->GetBinContent(2)+hists[i][3]->GetBinContent(3)+hists[i][3]->GetBinContent(4)+hists[i][3]->GetBinContent(5)+hists[i][3]->GetBinContent(6));
    
    Percentages[i][1]=(hists[i][3]->GetBinContent(3)+hists[i][3]->GetBinContent(4)+hists[i][3]->GetBinContent(5)+hists[i][3]->GetBinContent(6))/(hists[i][3]->GetBinContent(1)+hists[i][3]->GetBinContent(2)+hists[i][3]->GetBinContent(3)+hists[i][3]->GetBinContent(4)+hists[i][3]->GetBinContent(5)+hists[i][3]->GetBinContent(6));
    
    Percentages[i][2]=(hists[i][3]->GetBinContent(4)+hists[i][3]->GetBinContent(5)+hists[i][3]->GetBinContent(6))/(hists[i][3]->GetBinContent(1)+hists[i][3]->GetBinContent(2)+hists[i][3]->GetBinContent(3)+hists[i][3]->GetBinContent(4)+hists[i][3]->GetBinContent(5)+hists[i][3]->GetBinContent(6));
    
    ChargedParticleEvents[i]=1.0-hists[i][4]->GetBinContent(1)/(hists[i][4]->GetBinContent(1)+hists[i][4]->GetBinContent(2)+hists[i][4]->GetBinContent(3)+hists[i][4]->GetBinContent(4)+hists[i][4]->GetBinContent(5)+hists[i][4]->GetBinContent(6));
    
  }
  
  float var1=0.0;
  float var2=0.0;
  float var3=0.0;
  
  for(int i=1; i<8;i++){
    var1+=hists[0][5]->GetBinContent(i);
  }
  
  for(int i=1; i<22;i++){
    var2+=hists[0][6]->GetBinContent(i);
  }
  
  for(int i=1; i<36;i++){
    var3+=hists[0][7]->GetBinContent(i);
  }
  
  
  
  
  std::cout << "Efficiency of reconstructing tracks MC1: " << (WhetherdRMatch[0][1]/(WhetherdRMatch[0][1]+WhetherdRMatch[0][0])) << " Efficiency of reconstructing tracks MC2: " << (WhetherdRMatch[1][1]/(WhetherdRMatch[1][1]+WhetherdRMatch[1][0])) << std::endl;
  
  
  std::cout << "Efficiency of finding 2 prong in SV MC1: " << (WhetherSVMassMatch1[0][1]/(WhetherSVMassMatch1[0][1]+WhetherSVMassMatch1[0][0])) << " Efficiency of finding 2 prong in SV MC2: " << (WhetherSVMassMatch1[1][1]/(WhetherSVMassMatch1[1][1]+WhetherSVMassMatch1[1][0])) << std::endl;
  
  std::cout << "Efficiency of finding 3 prong in SV MC1: " << (WhetherSVThreeMassMatch1[0][1]/(WhetherSVThreeMassMatch1[0][1]+WhetherSVThreeMassMatch1[0][0])) << " Efficiency of finding 3 prong in SV MC2: " << (WhetherSVThreeMassMatch1[1][1]/(WhetherSVThreeMassMatch1[1][1]+WhetherSVThreeMassMatch1[1][0])) << std::endl;
  
  std::cout << "Percentage of MC1 events where matched tracks are less than 2: "<<Percentages[0][0]<<" Percentage of MC1 events where matched tracks are greater than or equal to two: "<<Percentages[0][1]<<" Percentage of MC1 events where matched tracks are greater than or equal to three: "<<Percentages[0][2]<<std::endl;
  
  std::cout << "Percentage of MC2 events where matched tracks are less than 2: "<<Percentages[1][0]<<" Percentage of MC2 events where matched tracks are greater than or equal to two: "<<Percentages[1][1]<<" Percentage of MC2 events where matched tracks are greater than or equal to three: "<<Percentages[1][2]<<std::endl;
  
  std::cout << "Percent of MC1 events where there are atleast one charged particle: "<<ChargedParticleEvents[0]<< " Percent of MC2 events where there are atleast one charged particle: "<<ChargedParticleEvents[1]<<std::endl;
  
  std::cout << " Efficiency of lowest Product1 being correct pair MC1: " << (hists[0][8]->GetBinContent(2)/(hists[0][8]->GetBinContent(1)+hists[0][8]->GetBinContent(2))) << " Efficiency of lowest Product1 being correct pair MC2: " << (hists[1][8]->GetBinContent(2)/(hists[1][8]->GetBinContent(1)+hists[1][8]->GetBinContent(2))) << std::endl;
  
  std::cout << "Percent of MC1 events where the lowest chi2 is correct (1-prong): "<<hists[0][5]->GetBinContent(1)/var1<<std::endl;
  std::cout << "Percent of MC1 events where the lowest chi2 is correct (2-prong): "<<hists[0][6]->GetBinContent(1)/var2<<std::endl;
  std::cout << "Percent of MC1 events where the lowest chi2 is correct (3-prong): "<<hists[0][7]->GetBinContent(1)/var3<<std::endl;
}