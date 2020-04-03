#include "TH1F.h"
#include "TH1D.h"
#include "TRandom.h"
#include <vector>
#include <iostream>

void comp2017vs2018_peak()
{

  std::vector<TH1D*> hists2018;
  std::vector<TH1D*> hists2017;

  TFile *f2018 = new TFile("../../LOCAL_COMBINED_dstophipi_default.root", "READ"); 
  TFile *f2017 = new TFile("../../../workdir_Decay_Length_2017_Mar_30_2020/LOCAL_COMBINED_dstophipi_default.root", "READ");

  TH1D *Muon1_Eta_2018, *Muon1_Eta_2017, *Muon1_Phi_2018, *Muon1_Phi_2017, *Muon1_Pt__2018, *Muon1_Pt__2017, 
       *Muon2_Eta_2018, *Muon2_Eta_2017, *Muon2_Phi_2018, *Muon2_Phi_2017, *Muon2_Pt__2018, *Muon2_Pt__2017,
       *Track_Eta_2018, *Track_Eta_2017, *Track_Phi_2018, *Track_Phi_2017, *Track_Pt__2018, *Track_Pt__2017,
       *DecayLength_2018, *DecayLength_2017, *Ds_L_2018, *Ds_L_2017, 
       *Iso1_2018, *Iso1_2017, *Iso1Mu1_2018, *Iso1Mu1_2017, *Iso8Mu1_2018, *Iso8Mu1_2017,
       *MaxD0SigPV_2018, *MaxD0SigPV_2017, *MaxD0SigSV_2018, *MaxD0SigSV_2017,
       *MinD0SigPV_2018, *MinD0SigPV_2017, *MinD0SigSV_2018, *MinD0SigSV_2017,
       *MaxDca_2018, *MaxDca_2017, *MinDca_2018, *MinDca_2017,
       *MaxdeltaMuZ_2018, *MaxdeltaMuZ_2017, *MaxVertexPairQuality_2018, *MaxVertexPairQuality_2017,
       *MindcaTrackSV_2018, *MindcaTrackSV_2017, *MinMuon_chi2LocalPosition_2018, *MinMuon_chi2LocalPosition_2017,
       *NSV_2018, *NSV_2017, *NtracksClose_2018, *NtracksClose_2017,
       *SVPVDsDirAngle_2018, *SVPVDsDirAngle_2017, *VertexKFChi2_2018, *VertexKFChi2_2017;

  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  c1->SetRightMargin(0.09);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  Muon1_Eta_2018 = (TH1D*)f2018->Get("dstophipi_default_Muon1_Eta_peakData"); hists2018.push_back(Muon1_Eta_2018);
  Muon1_Eta_2017 = (TH1D*)f2017->Get("dstophipi_default_Muon1_Eta_peakData"); hists2017.push_back(Muon1_Eta_2017);
  Muon1_Phi_2018 = (TH1D*)f2018->Get("dstophipi_default_Muon1_Phi_peakData"); hists2018.push_back(Muon1_Phi_2018);
  Muon1_Phi_2017 = (TH1D*)f2017->Get("dstophipi_default_Muon1_Phi_peakData"); hists2017.push_back(Muon1_Phi_2017);
  Muon1_Pt__2018 = (TH1D*)f2018->Get("dstophipi_default_Muon1_Pt_peakData");  hists2018.push_back(Muon1_Pt__2018);
  Muon1_Pt__2017 = (TH1D*)f2017->Get("dstophipi_default_Muon1_Pt_peakData");  hists2017.push_back(Muon1_Pt__2017);
  Muon2_Eta_2018 = (TH1D*)f2018->Get("dstophipi_default_Muon2_Eta_peakData"); hists2018.push_back(Muon2_Eta_2018);
  Muon2_Eta_2017 = (TH1D*)f2017->Get("dstophipi_default_Muon2_Eta_peakData"); hists2017.push_back(Muon2_Eta_2017);
  Muon2_Phi_2018 = (TH1D*)f2018->Get("dstophipi_default_Muon2_Phi_peakData"); hists2018.push_back(Muon2_Phi_2018);
  Muon2_Phi_2017 = (TH1D*)f2017->Get("dstophipi_default_Muon2_Phi_peakData"); hists2017.push_back(Muon2_Phi_2017);
  Muon2_Pt__2018 = (TH1D*)f2018->Get("dstophipi_default_Muon2_Pt_peakData");  hists2018.push_back(Muon2_Pt__2018);
  Muon2_Pt__2017 = (TH1D*)f2017->Get("dstophipi_default_Muon2_Pt_peakData");  hists2017.push_back(Muon2_Pt__2017);
  Track_Eta_2018 = (TH1D*)f2018->Get("dstophipi_default_Track_Eta_peakData"); hists2018.push_back(Track_Eta_2018);
  Track_Eta_2017 = (TH1D*)f2017->Get("dstophipi_default_Track_Eta_peakData"); hists2017.push_back(Track_Eta_2017);
  Track_Phi_2018 = (TH1D*)f2018->Get("dstophipi_default_Track_Phi_peakData"); hists2018.push_back(Track_Phi_2018);
  Track_Phi_2017 = (TH1D*)f2017->Get("dstophipi_default_Track_Phi_peakData"); hists2017.push_back(Track_Phi_2017);
  Track_Pt__2018 = (TH1D*)f2018->Get("dstophipi_default_Track_Pt_peakData");  hists2018.push_back(Track_Pt__2018);
  Track_Pt__2017 = (TH1D*)f2017->Get("dstophipi_default_Track_Pt_peakData");  hists2017.push_back(Track_Pt__2017);
  DecayLength_2018 = (TH1D*)f2018->Get("dstophipi_default_DecayLength_peakData"); hists2018.push_back(DecayLength_2018);
  DecayLength_2017 = (TH1D*)f2017->Get("dstophipi_default_DecayLength_peakData"); hists2017.push_back(DecayLength_2017);
  Ds_L_2018 = (TH1D*)f2018->Get("dstophipi_default_DsL_peakData"); hists2018.push_back(Ds_L_2018);
  Ds_L_2017 = (TH1D*)f2017->Get("dstophipi_default_DsL_peakData"); hists2017.push_back(Ds_L_2017);
  Iso1_2018 = (TH1D*)f2018->Get("dstophipi_default_Iso1_peakData"); hists2018.push_back(Iso1_2018);
  Iso1_2017 = (TH1D*)f2017->Get("dstophipi_default_Iso1_peakData"); hists2017.push_back(Iso1_2017);
  Iso1Mu1_2018 = (TH1D*)f2018->Get("dstophipi_default_Iso1Mu1_peakData"); hists2018.push_back(Iso1Mu1_2018);
  Iso1Mu1_2017 = (TH1D*)f2017->Get("dstophipi_default_Iso1Mu1_peakData"); hists2017.push_back(Iso1Mu1_2017);
  Iso8Mu1_2018 = (TH1D*)f2018->Get("dstophipi_default_Iso8Mu1_peakData"); hists2018.push_back(Iso8Mu1_2018);
  Iso8Mu1_2017 = (TH1D*)f2017->Get("dstophipi_default_Iso8Mu1_peakData"); hists2017.push_back(Iso8Mu1_2017);
  MaxD0SigPV_2018 = (TH1D*)f2018->Get("dstophipi_default_MaxD0SigPV_peakData"); hists2018.push_back(MaxD0SigPV_2018);
  MaxD0SigPV_2017 = (TH1D*)f2017->Get("dstophipi_default_MaxD0SigPV_peakData"); hists2017.push_back(MaxD0SigPV_2017);
  MaxD0SigSV_2018 = (TH1D*)f2018->Get("dstophipi_default_MaxD0SigSV_peakData"); hists2018.push_back(MaxD0SigSV_2018);
  MaxD0SigSV_2017 = (TH1D*)f2017->Get("dstophipi_default_MaxD0SigSV_peakData"); hists2017.push_back(MaxD0SigSV_2017);
  MinD0SigPV_2018 = (TH1D*)f2018->Get("dstophipi_default_MinD0SigPV_peakData"); hists2018.push_back(MinD0SigPV_2018);
  MinD0SigPV_2017 = (TH1D*)f2017->Get("dstophipi_default_MinD0SigPV_peakData"); hists2017.push_back(MinD0SigPV_2017);
  MinD0SigSV_2018 = (TH1D*)f2018->Get("dstophipi_default_MinD0SigSV_peakData"); hists2018.push_back(MinD0SigSV_2018);
  MinD0SigSV_2017 = (TH1D*)f2017->Get("dstophipi_default_MinD0SigSV_peakData"); hists2017.push_back(MinD0SigSV_2017);
  MaxDca_2018 = (TH1D*)f2018->Get("dstophipi_default_MaxDca_peakData"); hists2018.push_back(MaxDca_2018);
  MaxDca_2017 = (TH1D*)f2017->Get("dstophipi_default_MaxDca_peakData"); hists2017.push_back(MaxDca_2017);
  MinDca_2018 = (TH1D*)f2018->Get("dstophipi_default_MinDca_peakData"); hists2018.push_back(MinDca_2018);
  MinDca_2017 = (TH1D*)f2017->Get("dstophipi_default_MinDca_peakData"); hists2017.push_back(MinDca_2017);
  MaxdeltaMuZ_2018 = (TH1D*)f2018->Get("dstophipi_default_MaxdeltaMuZ_peakData"); hists2018.push_back(MaxdeltaMuZ_2018);
  MaxdeltaMuZ_2017 = (TH1D*)f2017->Get("dstophipi_default_MaxdeltaMuZ_peakData"); hists2017.push_back(MaxdeltaMuZ_2017);
  MaxVertexPairQuality_2018 = (TH1D*)f2018->Get("dstophipi_default_MaxVertexPairQuality_peakData"); hists2018.push_back(MaxVertexPairQuality_2018);
  MaxVertexPairQuality_2017 = (TH1D*)f2017->Get("dstophipi_default_MaxVertexPairQuality_peakData"); hists2017.push_back(MaxVertexPairQuality_2017);
  MindcaTrackSV_2018 = (TH1D*)f2018->Get("dstophipi_default_MindcaTrackSV_peakData"); hists2018.push_back(MindcaTrackSV_2018);
  MindcaTrackSV_2017 = (TH1D*)f2017->Get("dstophipi_default_MindcaTrackSV_peakData"); hists2017.push_back(MindcaTrackSV_2017);
  MinMuon_chi2LocalPosition_2018 = (TH1D*)f2018->Get("dstophipi_default_MinMuon_chi2LocalPosition_peakData"); hists2018.push_back(MinMuon_chi2LocalPosition_2018);
  MinMuon_chi2LocalPosition_2017 = (TH1D*)f2017->Get("dstophipi_default_MinMuon_chi2LocalPosition_peakData"); hists2017.push_back(MinMuon_chi2LocalPosition_2017);
  NSV_2018 = (TH1D*)f2018->Get("dstophipi_default_NSV_peakData"); hists2018.push_back(NSV_2018);
  NSV_2017 = (TH1D*)f2017->Get("dstophipi_default_NSV_peakData"); hists2017.push_back(NSV_2017);
  NtracksClose_2018 = (TH1D*)f2018->Get("dstophipi_default_NtracksClose_peakData"); hists2018.push_back(NtracksClose_2018);
  NtracksClose_2017 = (TH1D*)f2017->Get("dstophipi_default_NtracksClose_peakData"); hists2017.push_back(NtracksClose_2017);
  SVPVDsDirAngle_2018 = (TH1D*)f2018->Get("dstophipi_default_SVPVDsDirAngle_peakData"); hists2018.push_back(SVPVDsDirAngle_2018);
  SVPVDsDirAngle_2017 = (TH1D*)f2017->Get("dstophipi_default_SVPVDsDirAngle_peakData"); hists2017.push_back(SVPVDsDirAngle_2017);
  //VertexKFChi2_2018 = (TH1D*)f2018->Get("dstophipi_default_VertexKFChi2Data"); hists2018.push_back(VertexKFChi2_2018);
  //VertexKFChi2_2017 = (TH1D*)f2017->Get("dstophipi_default_VertexKFChi2Data"); hists2017.push_back(VertexKFChi2_2017);


  TLegend * legend = new TLegend(0.7,0.9,1,1);
  legend->AddEntry(hists2017[0],"2017 Data","lep");
  legend->AddEntry(hists2018[0],"2018 Data","lep");
  legend->Draw();

  c1->Print("2017vs2018comp.pdf[");

  for(int i=0;i<hists2018.size();i++){
    hists2018[i]->SetMarkerColor(kBlack);
    hists2017[i]->SetMarkerColor(kRed);

    hists2017[i]->Scale(hists2018[i]->Integral()/hists2017[i]->Integral());

    hists2018[i]->Draw();
    hists2017[i]->Draw("SAME");

    legend->Draw();

    c1->Print("2017vs2018comp.pdf");
  }

  c1->Print("2017vs2018comp.pdf]");

}
